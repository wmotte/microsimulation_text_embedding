#!/usr/bin/env Rscript
#
# Process glove embedding
#
################################################################################
library( "umap" )
library( "ggplot2" )
library( "ggrepel" )
library( "ggcorrplot" )

################################################################################
# FUNCTIONS
################################################################################

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
# Number of ticks
##
number_ticks <- function( n )
{
    function( limits )
        pretty( limits, n + 1 )
}

###
# Get section (i.e., groups)
##
get_sections <- function()
{
    # get group d.f.
    sections <- rbind(
        data.frame( label = c('Apples', 'Bananas', 'Grapefruit', 'Grapes', 'Kiwi', 'Lime', 'Mangoes', 
                              'Oranges', 'Pineapples', 'Strawberries', 'Watermelons' ), section = 'fruits' ),
        
        data.frame( label = c( 'Cauliflower', 'Cucumbers', 'Eggplant', 'Onion', 'Spinach', 
                               'Tomatoes', 'Peppers', 'Zucchini' ), section = 'vegetables' ),
        
        data.frame( label = c( 'BitterLemon', 'Cassis', 'Coffee', 'Cola', 'Lemonade', 
                               'Sodawater', 'Tea' ), section = 'beverages' ),
        
        data.frame( label = c( 'CaramelBars', 'ChewingGum', 'ChocolateBar', 'Crackers', 'GummyBears', 
                               'MixedNuts', 'Popcorn', 'PotatoChips', 'Snickers', 'Sweets', 'Snacks', 
                               'Twix', 'Pringles', 'ShoppingBag' ), section = 'extras' ) )
    
    return( sections )
}

###
# Get top-10 neighbors
##
calculate_neighbors <- function( embedding )
{
    # container
    all <- NULL
    
    vnames <- rownames( as.data.frame( embedding ) )
    vname <- vnames[ 1 ]
    
    for( vname in vnames )
    {
        # get single word vector
        single <- embedding[ vname, , drop = FALSE ]
        
        # pairwise similarities
        cos_sim <- text2vec::sim2( x = embedding, y = single, method = "cosine", norm = "l2" )
        
        # get highest correspond (n=10)
        corr <- sort( cos_sim[ , 1 ], decreasing = TRUE )[ 2:11 ]
        
        # into d.f.
        data <- data.frame( product = vname, closest_by = names( corr ), similarity = round( corr, 3 ) )
        
        # merge into container
        all <- rbind( all, data )
        
    }
    
    rownames( all ) <- NULL 
    return( all )
}

###
# Make heatmap
##
make_heatmap <- function( input_matrix, sname )
{
    # get max value to normalize matrix
    max_value <- max( abs( input_matrix ) )
    
    # normalize and transpose
    mat <- t( input_matrix / max_value )
    
    # get sections
    ss <- get_sections()
    
    # order as sections
    mat <- mat[ , rev( ss$label ) ]
    
    # long format
    dfmat <- reshape2::melt( mat, na.rm = TRUE )
    
    # matrix plot
    pmat <- ggplot( data = dfmat, mapping = aes_string( x = "Var1", y = "Var2", fill = "value" ) ) +
        geom_tile( color = 'white' ) + 
        scale_x_continuous( breaks = c( 1, 5, 10, 15, 20, 25, 30, 35, 40, 45, 50 ) ) +
        # color scheme, diverging with white in the middle
        scale_fill_gradient2( midpoint = 0, low = "#6D9EC1", mid = "white", high = "#E46726", space = "Lab" ) +
        xlab( 'embedding (dim)' ) + ylab( 'item' ) +
        theme_custom( 12 ) + theme( legend.position = 'none' )
    
    # save to disk
    ggsave( plot = pmat, dpi = 300, height = 8, width = 11, file = paste0( outdir, '/mat_', sname, '.png' ) )
}

################################################################################
# END FUNCTIONS
################################################################################

# load vocab, tcm, glove, wv_main, wv_context
load( "out.01.emb/saved_glove.RData" )

# output dir
outdir <- 'out.02.process'
dir.create( outdir, showWarnings = FALSE )

# plot matrix
#image( wv_main )    # 120 x 50
#image( wv_context ) # 50 x 120
#image( embedding )  # 120 x 50

# neighbors
nn <- calculate_neighbors( embedding )

# write to file
readr::write_csv( nn, file = paste0( outdir, '/nearest_neighbors.csv' ), quote = 'all' )

# loop over types
type <- 'I'

for( type in c( 'I', 'II', 'III' ) )
{
    identifier <- paste0( "_", type, "$" )
    
    # only select products 'I' (i.e., 40 x 50 matrix) -> wv_main
    wv_main_small <- wv_main[ grep( identifier, rownames( wv_main ) ), ]
    rownames( wv_main_small ) <- gsub( identifier, "", rownames( wv_main_small ) )
    
    # only select products 'I' (i.e., 40 x 50 matrix) -> wv_context
    twv_context <- t( wv_context )
    wv_context_small <- twv_context[ grep( identifier, rownames( twv_context ) ), ]
    rownames( wv_context_small ) <- gsub( identifier, "", rownames( wv_context_small ) )
    
    # only select products 'I' (i.e., 40 x 50 matrix) -> average of main and context
    embedding_small <- embedding[ grep( identifier, rownames( embedding ) ), ]
    rownames( embedding_small ) <- gsub( identifier, "", rownames( embedding_small ) )

    # write heat maps of all dimensions as matrices
    make_heatmap( wv_main_small, paste0( 'wv_main_type_', type ) )
    make_heatmap( wv_context_small, paste0( 'wv_context_type_', type ) )
    make_heatmap( embedding_small, paste0( 'wv_embedding_type_', type ) )
    
    # dimension reduction
    glove_umap <- umap( embedding_small, n_components = 2, metric = 'cosine', min_dist = 0.2 )

    # do the same for the GloVe embeddings
    df <- as.data.frame( glove_umap$layout, stringsAsFactors = FALSE )

    # Add the labels of the words to the data frame
    df$label <- rownames( embedding_small )
    colnames( df ) <- c( "UMAP1", "UMAP2", "label" )
    
    # merge with sections (groups)
    df <- merge( df, get_sections() )
    
    df$label <- as.factor( df$label )
    df$section <- as.factor( df$section )

    # get map
    p <- ggplot( df, aes( x = UMAP1, y = UMAP2 ) ) +
         geom_point( aes( fill = section ), shape = 21, colour = 'gray30', size = 3 ) +
         geom_label_repel( aes( label = label, fill = section ), color = 'white', segment.colour="gray30", size = 2.5, alpha = 0.8 ) +
         xlab( 'dimension I' ) +
         ylab( 'dimension II' ) +
         scale_x_continuous( breaks = number_ticks( 6 ) ) +
         scale_y_continuous( breaks = number_ticks( 6 ) ) +
        
         theme_bw( base_size = 12 ) %+replace% 
         theme( legend.position = "top",
               axis.ticks = element_blank(), 
               legend.background = element_blank(), 
               legend.key = element_blank(), 
               panel.border = element_blank(), 
               strip.background = element_blank(),
               strip.text.x = element_text( face = "bold" ), 
               strip.text.y = element_text( face = "bold" ),
               complete = FALSE )
    
    # save to disk
    ggsave( plot = p, dpi = 300, height = 7, width = 7, file = paste0( outdir, '/map_embedding_type_', type, '.png' ) )
}




