#!/usr/bin/env Rscript
#
# Visualize supermarket map as network
#
################################################################################
library( "igraph" )
library( "readxl" )


###
# Get general transition matrix
#
# 40 nodes
##
plot_supermarket <- function( outdir )
{
    # input matrix    
    input_matrix <- 'doc/transition_matrix.xlsx' 
    
    # read
    raw <- readxl::read_xlsx( input_matrix )
    raw$`...1` <- NULL
    
    # to df
    tmat <- as.data.frame( raw )
    
    for( i in 1:ncol( tmat ) )
        tmat[ , i ] <- as.numeric( tmat[ , i ] )   
    
    # to matrix
    tmat <- as.matrix( tmat )

    # relative contribution of weights
    p_backward <- 2 # lower triangle
    p_forward <- 10 # upper triangle
    p_diag <- 1     # remain in place
    
    tmat[ lower.tri( tmat ) ] <- tmat[ lower.tri( tmat ) ] * p_backward 
    tmat[ upper.tri( tmat, diag = FALSE ) ] <- tmat[ upper.tri( tmat, diag = FALSE ) ] * p_forward
    diag( tmat ) <- p_diag
    
    # normalize [row have to sum to 1!]
    correction_value <- 1 / rowSums( tmat, na.rm = TRUE ) 
    pmat <- tmat * correction_value 

    # set NA to 0
    pmat[ is.na( pmat ) ] <- 0
    
    # set diag to 0
    diag( pmat ) <- 0

    # absorption in prob. at each cycle
    absorption_per_cycle <- 1 - rowSums( pmat, na.rm = TRUE )
    
    # build the graph object
    g <- graph_from_adjacency_matrix( as.matrix( pmat ), weighted = TRUE )
    
    # get grid-kind of layout
    set.seed( 123 )    
    l <- layout_with_lgl( g )
    
    # reposition 12 and 13
    l[ 12, ] <- c( -50, -26.36 )
    l[ 13, ] <- c( -48.2, -31.8 )
    
    # colors
    V( g )$color <- c( rep( '#56bcc3', 11 ), 
                       rep( '#bd81f9', 8 ), 
                       rep( '#e07b71', 7 ),
                       rep( '#87ac34', 14 )  )
    
    # save to disk
    outfile <- paste0( outdir, '/supermarket_network.png' )
    png( outfile, height = 3200, width = 3200, res = 500 )
               
    # plot
    plot( g, layout = l, vertex.size = 10, edge.arrow.size = 0.4, vertex.label.cex = 0.5, 
          vertex.label.color = 'white', vertex.label.font = 2 )

    dev.off()
}

################################################################################



# output dir
outdir <- 'out.04.igraph'
dir.create( outdir, showWarnings = FALSE )

# plot
plot_supermarket( outdir )

