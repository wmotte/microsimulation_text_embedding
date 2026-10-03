#!/usr/bin/env Rscript
#
# W.M. Otte (w.m.otte@umcutrecht.nl)
#
#
# NOTE: *** 100k random walks are unique! ***
#
# Microsimulation of a supermarket model with shopping visitors: 
#
# Krijkamp EM, Alarid-Escudero F, Enns EA, Jalal HJ, Hunink MGM, Pechlivanoglou P. 
# Microsimulation modeling for health decision sciences using R: A tutorial. 
# Med Decis Making. 2018;38(3):400–22.
# 
# See GitHub for more information or code updates
# https://github.com/DARTH-git/Microsimulation-tutorial

library( 'ggplot2' )

################################################################################
##################################### Functions ################################
################################################################################

###
# Random sample from Bernoulli distribution
##
rbernoulli <- function( n, p = 0.5 ) 
{
    return( stats::runif( n ) > ( 1 - p ) )
}

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

# The MicroSim function for the microsimulation
#
# Arguments:  
# general_transition_metrix: probabilities for state transitions
# v.M_1:   vector of initial states for individuals 
# n.i:     number of individuals
# n.t:     total number of cycles to run the model
# v.n:     vector of state names
#
# TR.out:  should the output include a microsimulation trace?
# TS.out:  should the output include a matrix of transitions between states?
# seed:    starting seed number for random number generator (default is 1)
#
##
MicroSim <- function( general_transition_matrix, v.M_1, n.i, n.t, v.n, TR.out = TRUE, TS.out = FALSE, seed = 1 )
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
    
    return( results )
    
}  # end of the MicroSim function  

###
# Process data for plotting
##
process_plotting_data <- function( df )
{
    
    # convert to long [Cycle x State x value]
    data <- as.data.frame( df$TR )
    data$time <- as.factor( 1:nrow( data ) )
    rownames( data ) <- NULL
    
    # get long
    data_long <- reshape2::melt( data, id.vars = 'time' )
    data_long$time <- as.numeric( data_long$time )
    
    return( data_long )
}

###
# Return random sets of 'I', 'II' or 'III'
##
get_product_version <- function( n )
{
    result <- sample( x = c( 'I', 'II', 'III' ), n, replace = TRUE )  
    return( result )
}

###
# Place random products (1,2 or 3)
##
place_products <- function( vchain )
{
    for( type in c( 'I', 'II', 'III' ) )
    {
        # convert text to products
        vchain[ vchain %in% paste0( 'H1', '_', type ) ] <- paste0( 'Apples', '_', type )
        vchain[ vchain %in% paste0( 'H2', '_', type ) ] <- paste0( 'Bananas', '_', type )
        vchain[ vchain %in% paste0( 'H3', '_', type ) ] <- paste0( 'Grapefruit', '_', type )
        vchain[ vchain %in% paste0( 'H4', '_', type ) ] <- paste0( 'Grapes', '_', type )
        vchain[ vchain %in% paste0( 'H5', '_', type ) ] <- paste0( 'Kiwi', '_', type )
        vchain[ vchain %in% paste0( 'H6', '_', type ) ] <- paste0( 'Lime', '_', type )
        vchain[ vchain %in% paste0( 'H7', '_', type ) ] <- paste0( 'Mangoes', '_', type )
        vchain[ vchain %in% paste0( 'H8', '_', type ) ] <- paste0( 'Oranges', '_', type )
        vchain[ vchain %in% paste0( 'H9', '_', type ) ] <- paste0( 'Pineapples', '_', type )
        vchain[ vchain %in% paste0( 'H10', '_', type ) ] <- paste0( 'Strawberries', '_', type )
        vchain[ vchain %in% paste0( 'H11', '_', type ) ] <- paste0( 'Watermelons', '_', type )
        
        # vegetables
        vchain[ vchain %in% paste0( 'H12', '_', type ) ] <- paste0( 'Cauliflower', '_', type )
        vchain[ vchain %in% paste0( 'H13', '_', type ) ] <- paste0( 'Cucumbers', '_', type )
        vchain[ vchain %in% paste0( 'H14', '_', type ) ] <- paste0( 'Eggplant', '_', type )
        vchain[ vchain %in% paste0( 'H15', '_', type ) ] <- paste0( 'Onion', '_', type )
        vchain[ vchain %in% paste0( 'H16', '_', type ) ] <- paste0( 'Spinach', '_', type )
        vchain[ vchain %in% paste0( 'H17', '_', type ) ] <- paste0( 'Tomatoes', '_', type )
        vchain[ vchain %in% paste0( 'H18', '_', type ) ] <- paste0( 'Peppers', '_', type )
        vchain[ vchain %in% paste0( 'H19', '_', type ) ] <- paste0( 'Zucchini', '_', type )
        
        # beverages
        vchain[ vchain %in% paste0( 'H20', '_', type ) ] <- paste0( 'BitterLemon', '_', type )
        vchain[ vchain %in% paste0( 'H21', '_', type ) ] <- paste0( 'Cassis', '_', type )
        vchain[ vchain %in% paste0( 'H22', '_', type ) ] <- paste0( 'Coffee', '_', type )
        vchain[ vchain %in% paste0( 'H23', '_', type ) ] <- paste0( 'Cola', '_', type )
        vchain[ vchain %in% paste0( 'H24', '_', type ) ] <- paste0( 'Lemonade', '_', type )
        vchain[ vchain %in% paste0( 'H25', '_', type ) ] <- paste0( 'Sodawater', '_', type )
        vchain[ vchain %in% paste0( 'H26', '_', type ) ] <- paste0( 'Tea', '_', type )
        
        # extras
        vchain[ vchain %in% paste0( 'H27', '_', type ) ] <- paste0( 'CaramelBars', '_', type )
        vchain[ vchain %in% paste0( 'H28', '_', type ) ] <- paste0( 'ChewingGum', '_', type )
        vchain[ vchain %in% paste0( 'H29', '_', type ) ] <- paste0( 'ChocolateBar', '_', type )
        vchain[ vchain %in% paste0( 'H30', '_', type ) ] <- paste0( 'Crackers', '_', type )
        vchain[ vchain %in% paste0( 'H31', '_', type ) ] <- paste0( 'GummyBears', '_', type )
        vchain[ vchain %in% paste0( 'H32', '_', type ) ] <- paste0( 'MixedNuts', '_', type )
        vchain[ vchain %in% paste0( 'H33', '_', type ) ] <- paste0( 'Popcorn', '_', type )
        vchain[ vchain %in% paste0( 'H34', '_', type ) ] <- paste0( 'PotatoChips', '_', type )
        vchain[ vchain %in% paste0( 'H35', '_', type ) ] <- paste0( 'Snickers', '_', type )
        vchain[ vchain %in% paste0( 'H36', '_', type ) ] <- paste0( 'Sweets', '_', type )
        vchain[ vchain %in% paste0( 'H37', '_', type ) ] <- paste0( 'Snacks', '_', type )
        vchain[ vchain %in% paste0( 'H38', '_', type ) ] <- paste0( 'Twix', '_', type )
        vchain[ vchain %in% paste0( 'H39', '_', type ) ] <- paste0( 'Pringles', '_', type )        
        vchain[ vchain %in% paste0( 'H40', '_', type ) ] <- paste0( 'ShoppingBag', '_', type )
        
        # end of line
        vchain[ vchain %in% paste0( 'H41', '_', type ) ] <- paste0( 'Exit', '_', type )
    }
    
    return( vchain )
}

################################################################################
############################## End of functions ################################
################################################################################

# set seed
set.seed( 4321 )

# output
outdir <- 'out.00.sim'
dir.create( outdir, showWarnings = FALSE )

# get transition matrix (from Excel sheet)
general_transition_matrix <- get_general_transition_matrix()

# write transition matrix to disk (S1 Table; H41 is the absorbing exit state, not a shelf)
tm_out <- data.frame( from = paste0( "H", 1:nrow( general_transition_matrix ) ), round( general_transition_matrix, 4 ) )
colnames( tm_out ) <- c( 'from', paste0( "H", 1:ncol( general_transition_matrix ) ) )
readr::write_tsv( tm_out, file = paste0( outdir, '/transition_matrix.tsv' ) )

# Model input
n.i   <- 105000                # number of simulated individuals
n.t   <- 75                    # time horizon in cycles
v.n   <- paste0( "H", 1:nrow( general_transition_matrix ) )    # model state names
v.M_1 <- rep( "H1", n.i )      # everyone begins in the healthy state 

# get micro-simulations
df <- MicroSim( general_transition_matrix, v.M_1, n.i, n.t, v.n ) # run for no treatment

# save simulation to disk
save( df, file = paste0( outdir, "/simulation_df.RData" ) )

# selected long-format data
data_long <- process_plotting_data( df )

# select subset of network nodes to visualize
selected_nodes <- c( "H1", paste0( "H", seq( 3, 40, 5 ) ), "H41" )

# subset of nodes
sdata <- data_long[ data_long$variable %in% selected_nodes, ]

# save sdata
readr::write_csv( sdata, file = gzfile( paste0( outdir, '/network_flow.csv.gz' ) ), quote = 'all' )

# plot (H41 is the absorbing exit state, not a shelf: label it as such)
sdata$panel <- factor( ifelse( sdata$variable == 'H41', 'Exit (H41, absorbing state)', as.character( sdata$variable ) ),
                       levels = c( setdiff( selected_nodes, 'H41' ), 'Exit (H41, absorbing state)' ) )
p <- ggplot( data = sdata, aes( x = as.numeric( time ), y = value * 100, group = panel, colour = panel ) ) + 
    geom_line() + facet_wrap( ~panel, ncol = 2, scales = 'free' ) + theme_custom() + 
    scale_x_continuous( breaks = number_ticks( 6 ) ) +
    scale_y_continuous( breaks = number_ticks( 3 ) ) +
    xlab( 'Time step' ) + ylab( 'Shoppers at this location (%)' ) + theme( legend.position = 'none' )

# save to disk
ggsave( plot = p, dpi = 300, height = 8, width = 8, file = paste0( outdir, '/network_flow.png' ) )


######################################
######### sample products ############
######################################

# three product versions per node

# <n_subjects> x <n_cycles>
dim( mat <- df$m.M )

# get cycle length for each subject [76 means: still in supermarket]
cycle_length <- 76 - matrixStats::rowCounts( mat, value = 'H41' )

# get some numbers of length of stay of subjects in supermarket
stats_cycles <- data.frame( min_states_before_exit = min( cycle_length ),
            mean_states_before_exit = round( mean( cycle_length[ cycle_length != 76 ] ), 0 ),
            max_states_before_exit = max( cycle_length ) )

# write to disk
readr::write_tsv( stats_cycles, file = paste0( outdir, '/cycle_length__stats.tsv' ), quote = 'all' )

# write to disk
readr::write_tsv( data.frame( cycle_length ), file = paste0( outdir, '/cycle_length.tsv' ), quote = 'all' )

# plot histogram with percentages
p_cycle_l <- 
    ggplot( data = data.frame( x = cycle_length ), aes( x = x ) ) + 
    geom_histogram( aes( y = (..count..) / sum(..count..) ),
    colour = 'gray30', bins = 12, fill = 'orange' ) + theme_custom() + 
    scale_x_continuous( breaks = number_ticks( 12 ) ) +
    scale_y_continuous( breaks = number_ticks( 10 ), labels = scales::percent ) +
    xlab( 'Locations visited before exit (76 = still inside after step 75)' ) + ylab( 'Shoppers (%)' ) + theme( legend.position = 'none' )

# save to disk
ggsave( plot = p_cycle_l, dpi = 300, height = 8, width = 8, file = paste0( outdir, '/states_before_exit.png' ) )

# probability of product taking
prob <- 0.15

# string container (list instead of rbind: same random draws and output, but linear time)
all <- vector( 'list', nrow( mat ) )

# keep track of the source subject of each list
all_id <- rep( NA, nrow( mat ) )

i <- 1

# loop over subjects
for( i in 1:nrow( mat ) )
{
    if( ( i %% 1000 ) == 0 )
        print( i )
    
    # get subject route
    vsub <- mat[ i, ]
    
    # remove EXIT nodes
    vsub <- vsub[ vsub != 'H41' ]
    
    # get product taken
    vproducts <- vsub[ rbernoulli( length( vsub ), p = prob ) ]

    # only continue if at least six products are available
    if( length( vproducts ) > 5 )
    {
        # add type (I, II or III) to product
        vvs <- get_product_version( length( vproducts ) )
        vproducts <- paste0( vproducts, '_', vvs )
        
        vchain <- vproducts
        
        # place random products [i.e., H1 -> Apples_I, Apples_II, Apples_III ]
        vchain <- place_products( vchain )
        
        # collapse
        vsubject <- paste( vchain, collapse = " " )

        all[[ i ]] <- vsubject
        all_id[ i ] <- i
    }
}
all_id <- all_id[ !sapply( all, is.null ) ]
all <- unlist( all )

# remove words < 5
nwords <- stringr::str_count( all, ' ' ) + 1
tmp <- all[ nwords >= 5 ]
tmp_id <- all_id[ nwords >= 5 ]

# remove duplicated random walks
is_dup <- duplicated( tmp )
tmp <- tmp[ !is_dup ]
tmp_id <- tmp_id[ !is_dup ]

# select complete number 
final_set <- tmp[ 1:80000 ]
final_id <- tmp_id[ 1:80000 ]

# record what happened to all simulated subjects
filter_counts <- data.frame( 
    simulated = nrow( mat ),
    still_inside_at_step_75 = sum( mat[ , ncol( mat ) ] != 'H41' ),
    discarded_lt_6_products = nrow( mat ) - length( all ),
    duplicated_lists = sum( is_dup ),
    unique_lists = length( tmp ),
    not_used_surplus = length( tmp ) - 80000,
    in_corpus = length( final_set ),
    in_corpus_still_inside_at_step_75 = sum( mat[ final_id, ncol( mat ) ] != 'H41' ),
    tokens = sum( stringr::str_count( final_set, ' ' ) + 1 ) )
print( t( filter_counts ) )
readr::write_tsv( filter_counts, file = paste0( outdir, '/filter_counts.tsv' ) )

# example sentences (Box 1)
readr::write_lines( final_set[ 1:10 ], file = paste0( outdir, '/example_sentences.txt' ) )

# get summary
nwords_final <- stringr::str_count( final_set, ' ' ) + 1
length( final_set )
summary( nwords_final )

# plot density
p_sen <- 
    ggplot( data = data.frame( x = nwords_final ), aes( x = x ) ) + 
    geom_histogram( colour = 'gray30', bins = max( nwords ) - 5, fill = 'orange' ) + theme_custom() + 
    scale_x_continuous( breaks = number_ticks( max( nwords ) - 5 ) ) +
    scale_y_continuous( breaks = number_ticks( 10 ) ) +
    xlab( 'Products in list (N)' ) + ylab( 'Number of lists' ) + theme( legend.position = 'none' )

# save to disk
ggsave( plot = p_sen, dpi = 300, height = 8, width = 8, file = paste0( outdir, '/product_list_length.png' ) )

# write to plain text file
readr::write_lines( final_set, file = gzfile( paste0( outdir, '/plain_text.txt.gz' ) ) )


