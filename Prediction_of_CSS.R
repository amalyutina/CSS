#load data ( primary target matrix, SEA target matrix, fingerprint matrix ) to compare different regression methods  ( based on 4 statistics )
#for a selected cell line

library( caret )
library( dplyr )
library( tidyr )
library( plyr )
library( grid )
library( tibble )
library( data.table )
library( reshape2 )
library( readr )
library( stringr )
library( glmnet )


main<- function(  identifier, data, response, design, method, cv_p, num_seeds ) 
{ #now I do not inlcude checks for input parameters, because it is a test still, but in the future it will be added 
  
  collect_stats<- function( feature_matrix, response, method, cv_p, seeds, design_name ) 
  { #function to calculate statistics
    
    cl_stat <- list()
    cl_coef <- list()
    cl_imp <- list()
    res_stat <- list()
    
    #calculate separately for every seed
    for ( i in 1:length( seeds ) )  {
      
      cat("Seed: ", i, '\n')
      
      #split the data
      set.seed( seeds[i] )
      intrain <- createDataPartition( y = response, p = cv_p, list = FALSE )
      
      training <- feature_matrix[intrain,]
      testing <- feature_matrix[-intrain,]
      
      
      set.seed( seeds[i] )
      trctrl <- trainControl( method = "cv", number = 10 )
      
      
      #select the method
      if ( method == 'glmnet' ) {
        
        # elasticnet
       
        set.seed( seeds[i] )
        my_fit = train( training, as.vector( response[intrain] ), method = method,
                        trControl=trctrl )

       
      }
      
      if ( method == 'ranger' ) {
        
        set.seed( seeds[i] )
        my_fit <- train( training,as.vector( response[intrain] ), method = method, trControl=trctrl )
        
      }
      
      
      if ( method == 'lm' ||  method == 'svmRadial' ) {
        
        set.seed( seeds[i] )
        my_fit <- train( training, as.vector( response[intrain] ), method = method, trControl=trctrl )
        
      }
      
      #ridge is sensitive to zero variance columns
      if ( method == 'ridge' ) {
        
        zer_var <- nearZeroVar( training, saveMetrics = TRUE) 
        ind <- which( zer_var$zeroVar == 'TRUE' )
        
        if ( length( ind ) > 0 ) {
          
          cat("Removing zero variance columns: ", rownames( zer_var )[ind], '\n')
          
          set.seed( seeds[i] )
          new_training <- training[, -nearZeroVar( training ) ]
          testing <- testing[, -nearZeroVar( training ) ]
          my_fit <- train( new_training, as.vector( response[intrain] ), method = method, trControl=trctrl )
          
          
        } else {
          
          set.seed( seeds[i] )
          my_fit <- train( training, as.vector( response[intrain] ), method = method, trControl=trctrl )
          
        }
      }
      
      
      
      #build_predictions
      my_Predict <- predict( my_fit, newdata = testing )
      
      my_cor <- cor( my_Predict, as.vector( response[-intrain] ) )
      my_mean <- mean( abs( my_Predict - as.vector( response[-intrain] ) ) )
      my_rmse <- sqrt( mean( ( my_Predict - as.vector( response[-intrain] ) )^2 ) )
      my_r2 <- 1 - ( ( sum( ( my_Predict - as.vector( response[-intrain] ) )^2  ) )/( sum( ( as.vector( response[-intrain] )- mean(as.vector( response[-intrain]) ) )^2  ) ) ) 
      
      
      #MODEL
      tmp_coeffs <- coef(my_fit$finalModel, s = my_fit$bestTune$lambda)
    
      cl_stat[[i]] <- data.frame( cor = my_cor, MAE = my_mean, R2 = my_r2 , RMSE = my_rmse, seed = seeds[i], method = method, design_num = j, design = design_name  ) 
      
      #save coefficients
      cl_coef[[i]] <- data.frame( Feature = as.character( tmp_coeffs@Dimnames[[1]][tmp_coeffs@i + 1]), coef = tmp_coeffs@x, seed = seeds[i], method = method, design_num = j, design = design_name )
      #save importance
      cl_imp[[i]] <- data.frame( Feature = as.character( tmp_coeffs@Dimnames[[1]]), importance = c( 0, as.numeric( varImp(my_fit)$importance$Overall  )  ), seed = seeds[i], method = method, design_num = j, design = design_name )   
      
    }
    
    cat('\n')
    res_stat <- ldply( cl_stat, data.frame )
    res_coef <- ldply( cl_coef, data.frame )
    res_imp <-  ldply( cl_imp, data.frame )
    r <- list()
    
    r[[1]] <- res_stat
    r[[2]] <- res_coef
    r[[3]] <- res_imp
    r
  }
  
  #print information about the user selections
  print( identifier )
  cat( "method: ", method, "\n" )
  cat( "# of seeds: ", num_seeds, "\n" )
  cat('\n')
  
  #save results for every design option
  res_list <- list()  
  res_coef <- list()
  res_imp <- list()
  
  #loop over different design options  
  for ( j in 1:nrow( design ) )  
  {
    cat("Design: ", j , '\n')
    
    #find matrix or matrix combo to use
    ind <- which( design[j,] == 1 ) 
    
    design_name <- paste( names( data )[ind], collapse = "+")
    ####
    cat("Design name: ", design_name , '\n')
    ####
    design_data <- do.call( cbind, data[ind] )
    
    indi <- which( colSums( design_data ) == 0 | colSums( design_data ) == ncol( design_data ) )
    
    if ( length( indi ) >0 ) {
      design_data <- design_data[, -indi ] 
    }
    
    #specify the seeds
    seeds <- seq( 50, 10000, length.out = num_seeds )
    
    #collect statistics for every design specified by user
    f <- collect_stats( design_data, response, method, cv_p, seeds, design_name )
    
    res_list[[j]] <- f[[1]]
    res_coef[[j]] <- f[[2]]
    res_imp[[j]] <- f[[3]]
  }
  
  res <- ldply( res_list, data.frame )
  res_list <- list( identifier = identifier, res_statisics = res, method = method, coef = ldply( res_coef, data.frame ), imp = ldply( res_imp, data.frame )  )
}
plot_summary <- function( data_list) 
{ #function that plots the statistics; 
  
  multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) 
  { #function to plot several figures in one plot 
    
    main_title <- c( paste( 'method :', method ) )
    
    for ( i in 1:length( data_list[[1]] ) ) {
      
      main_title <- paste( main_title, ';', paste( names( data_list[[1]] )[i], ':', data_list[[1]][i] ) )
      
      
    }
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
      # Make the panel
      # ncol: Number of columns of plots
      # nrow: Number of rows needed, calculated from # of cols
      layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                       ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
      print(plots[[1]])
      
    } else {
      # Set up the page
      grid.newpage()
      pushViewport( viewport( layout = grid.layout( nrow(layout)+1, ncol(layout), heights = unit(c(0.5, 5), "null") ) ) ) 
      
      # Make each plot, in the correct location
      for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
        
        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row+1,
                                        layout.pos.col = matchidx$col))
        
      }
      
    }
    
    grid.text( main_title, vp = viewport(layout.pos.row = 1, layout.pos.col = 1:4 ) )
  }
  
  method <- as.character( data_list[[3]] )
  plot_data <- data_list[[2]]
  
  icpl <-  ggplot2::ggplot( plot_data, aes( seed, cor, colour = design ) ) +geom_point( size = 1.5, alpha = 0.7) + theme_bw() +
    labs(  x = "seeds ") + theme(plot.title = element_text(hjust = 0.5)) +theme(legend.text=element_text( size=4 ) ) +
    theme( axis.text.x = element_text( angle = 0, hjust = 1 ), axis.title.y=element_blank() )+ theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  icpl1 <- icpl + ggtitle( 'Correlation' ) 
  
  
  icpl <-  ggplot2::ggplot( plot_data, aes( seed, MAE, colour = design ) ) +geom_point( size = 1.5, alpha = 0.7) + theme_bw() +
    labs(  x = "seeds ") + theme(plot.title = element_text(hjust = 0.5)) +
    theme( axis.text.x = element_text( angle = 0, hjust = 1 ),  axis.title.y=element_blank()  )+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position="none")
  
  icpl4 <- icpl + ggtitle( 'MAE' )
  
  
  
  icpl <-  ggplot2::ggplot( plot_data, aes( seed, RMSE, colour = design ) ) +geom_point( size = 1.5, alpha = 0.7) + theme_bw() +
    labs(  x = "seeds ") + theme(plot.title = element_text(hjust = 0.5)) +
    theme( axis.text.x = element_text( angle = 0, hjust = 1 ), axis.title.y=element_blank()  )+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position="none")
  
  icpl2 <- icpl + ggtitle( 'RMSE' )
  
  
  icpl <-  ggplot2::ggplot( plot_data, aes( seed, R2, colour = design ) ) +geom_point( size = 1.5, alpha = 0.7) + theme_bw() +
    labs(  x = "seeds ") + theme(plot.title = element_text(hjust = 0.5)) +
    theme( axis.text.x = element_text( angle = 0, hjust = 1 ), axis.title.y=element_blank()  )+ theme(axis.text.x = element_text(angle = 90, hjust = 1))+ theme(legend.position="none")
  
  icpl3 <- icpl + ggtitle( 'R2' )
  
  
  multiplot(icpl3, icpl2, icpl4, icpl1, cols=4)
  
}


ii = 10#select the cell line number

# download data
load( 'finger_data.RData' ) # fpt - fingerprint matrix for all single and drug combos
load( 'primary_data.RData' ) # prim_data -  primary target data for all single and drug combos
load( 'sea_data.RData' ) # sea_data - sea target data for all single and drug combos
load( 'CSS.RData' ) # load CSS or dCSS



# cell line names 

cl_names <-  unique( gsub("^(.*?)_.*", "\\1", rownames( CSS ) ) )

# put matrices in a list object
list_of_matrices <- list( primary_targets = pri_mat, fingerprint = finger_mat,  sea_targets = sea_mat )

# select only one cell line for the test 
sel_cell_line <- as.numeric( ii )

#use only the rows that correspond to the selected cell line in the matrices and CSS response vector
ind <- which( stringr::str_detect( rownames( list_of_matrices[[1]] ), cl_names[sel_cell_line] ) == TRUE ) # find indexes that correspond to KMS-11 matrix ( all matrices have same rownames )
list_of_matrices <- lapply( list_of_matrices, function( x )  x <- x[ind, ]) 
CSS <- as.vector( CSS[ind,1] )


#make a design matrix 

design_matrix <- matrix( rbind( c( 1,0,0 ), c( 0,1,0 ), c( 0,0,1 ), c( 0,1,1), c(1,1,0) ), nrow = 5 )


#names for plots
identifier  <- list( cell_line = cl_names[sel_cell_line], response = 'CSS' )

#call main function and collect statistics
res <- main( identifier, list_of_matrices, CSS, design_matrix, 'glmnet', 0.7, 5 )
#plot the results 
plot_summary( res )

#save( res, file = paste( paste( 'res',cl_names[sel_cell_line], sep = '_'),"RData",sep=".") )



