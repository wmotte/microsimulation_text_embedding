#
# This code forms the basis for the microsimulation model of the article: 
#
# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22.
# 
# See GitHub for more information or code updates
# https://github.com/DARTH-git/Microsimulation-tutorial

library( 'ggplot2' )


##################################### Functions ###########################################

###
# Number of ticks
##
number_ticks <- function( n )
{
    function( limits )
        pretty( limits, n + 1 )
}

###
# Custom theme
##
theme_custom <- function( base_size = 16, base_family = "", 
                          base_line_size = base_size / 22, 
                          base_rect_size = base_size / 22 ) 
{
    theme_bw( base_size = base_size, base_family = base_family, 
              base_line_size = base_line_size, 
              base_rect_size = base_rect_size ) %+replace% 
        theme( legend.position = "top",
               axis.ticks = element_blank(), 
               legend.background = element_blank(), 
               legend.key = element_blank(), 
               panel.border = element_blank(), 
               strip.background = element_blank(),
               strip.text.x = element_text( face = "bold" ), 
               strip.text.y = element_text( face = "bold" ),
               complete = TRUE )
}

###
# Get general transition matrix
#
# 40 nodes
##
get_general_transition_matrix <- function()
{
    library( 'readxl' )
    
    input_matrix <- 'doc/transition_matrix.xlsx' 
    
    raw <- readxl::read_xlsx( input_matrix )
    raw$`...1` <- NULL
    
    tmat <- as.data.frame( raw )
    
    for( i in 1:ncol( tmat ) )
        tmat[ , i ] <- as.numeric( tmat[ , i ] )   
    
    tmat <- as.matrix( tmat )
    
    #library( 'igraph' )
    
    # build the graph object
    #network <- graph_from_adjacency_matrix( as.matrix( tmat ) )
    
    # plot it
    #plot( network, vertex.size = 8, vertex.label = NA )
    
    
    # TODO
    #state_names <- paste0( "H", 1:5 )
    
    # transition table [all potential transitions in the model (i.e., all arrows, in this case 8 arrows]
    #tmat <- rbind( c( NA,  1,  1,  1,  0 ),
    #               c( -1, NA,  1,  0,  0 ),
    #               c( -1, -1, NA,  1,  1 ),
    #               c( -1, -1, -1, NA,  1 ),
    #               c( NA, NA, NA, NA,  1 ) )
    
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
    
    rowSums( pmat, na.rm = TRUE )
    
    # set NA to 0
    pmat[ is.na( pmat ) ] <- 0
    
    return( pmat )
}

# The MicroSim function for the simple microsimulation of the 'Sick-Sicker' 
#
# Arguments:  
# v.M_1:   vector of initial states for individuals 
# n.i:     number of individuals
# n.t:     total number of cycles to run the model
# v.n:     vector of health state names
#
# TR.out:  should the output include a microsimulation trace? (default is TRUE)
# TS.out:  should the output include a matrix of transitions between states? (default is TRUE)
# seed:    starting seed number for random number generator (default is 1)
#
# Makes use of:
# Probs:   function for the estimation of transition probabilities
##
MicroSim <- function( general_transition_matrix, v.M_1, n.i, n.t, v.n, TR.out = TRUE, TS.out = TRUE, seed = 1 )
{

    # create the matrix capturing the state name/costs/health outcomes for all individuals at each time point 
    m.M <- matrix( nrow = n.i, ncol = n.t + 1, 
                                  dimnames = list( paste( "ind", 1:n.i, sep = " " ), 
                                                   paste( "cycle", 0:n.t, sep = " " ) ) )  
    
    # indicate the initial health state 
    m.M[ , 1 ] <- v.M_1                       
    
    # loop over subject
    for ( i in 1:n.i )
    {
        # set seed
        set.seed( seed + i )

        # loop over time
        for( t in 1:n.t ) 
        {
            M_it <- m.M[ i, t ]
            # get probs
            
            # get transition probability from matching row index
            idx <- colnames( general_transition_matrix ) %in% M_it
            v.p <- general_transition_matrix[ idx, ]
            
            # sample the next health state and store that state in matrix m.M 
            m.M[ i, t + 1 ] <- sample( v.n, prob = v.p, size = 1 )  
            
        } # close the loop for the time points 
        
        # print progress
        if( i / 100 == round( i / 100, 0 ) ) 
            cat( '\r', paste( i / n.i * 100, "% done", sep = " " ) )
        
    } # close the loop for the individuals 
    
    if (TS.out == TRUE) {  # create a  matrix of transitions across states
        TS <- paste(m.M, cbind(m.M[, -1], NA), sep = "->") # transitions from one state to the other
        TS <- matrix(TS, nrow = n.i)
        rownames(TS) <- paste("Ind",   1:n.i, sep = " ")   # name the rows 
        colnames(TS) <- paste("Cycle", 0:n.t, sep = " ")   # name the columns 
    } else {
        TS <- NULL
    }
    
    if (TR.out == TRUE) { # create a trace from the individual trajectories
        TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE))))
        TR <- TR / n.i                                       # create a distribution trace
        rownames(TR) <- paste("Cycle", 0:n.t, sep = " ")     # name the rows 
        colnames(TR) <- v.n                                  # name the columns 
    } else {
        TR <- NULL
    }
    
    # store the results from the simulation in a list  
    results <- list( m.M = m.M, TS = TS, TR = TR ) 
    
    return(results)
    
}  # end of the MicroSim function  

##################################### Run the simulation ##################################

# set seed
set.seed( 123 )

# output
outdir <- 'out.00.sim'
dir.create( outdir, showWarnings = FALSE )

# get transition matrix (from Excel sheet)
general_transition_matrix <- get_general_transition_matrix()

# Model input
n.i   <- 1500                  # number of simulated individuals
n.t   <- 250                   # time horizon in cycles
v.n   <- paste0( "H", 1:nrow( general_transition_matrix ) )    # model state names
v.M_1 <- rep( "H1", n.i )      # everyone begins in the healthy state 

# get micro-simulations
df <- MicroSim( general_transition_matrix, v.M_1, n.i, n.t, v.n ) # run for no treatment

# convert to long [Cycle x State x value]
data <- as.data.frame( df$TR )
data$time <- as.factor( 1:nrow( data ) )
rownames( data ) <- NULL

# get long
data_long <- reshape2::melt( data, id.vars = 'time' )
data_long$time <- as.numeric( data_long$time )

# select subset of network nodes to visualize
selected_nodes <- c( "H1", paste0( "H", seq( 3, 40, 5 ) ), "H40" )

# subset of nodes
sdata <- data_long[ data_long$variable %in% selected_nodes, ]

# save sdata
readr::write_csv( sdata, file = gzfile( paste0( outdir, '/network_flow.csv.gz' ) ), quote = 'all' )

# plot
p <- ggplot( data = sdata, aes( x = as.numeric( time ), y = value * 100, group = variable, colour = variable ) ) + 
    geom_line() + facet_wrap( ~variable, ncol = 2, scales = 'free' ) + theme_custom() + 
    xlab( 'time' ) + ylab( '%' ) + theme( legend.position = 'none' )

# save to disk
ggsave( plot = p, dpi = 300, height = 8, width = 8, file = paste0( outdir, '/network_flow.png' ) )



### sample products ###









