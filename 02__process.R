#!/usr/bin/env Rscript
#
# Process glove embedding
#
################################################################################
library( "umap" )
library( "ggplot2" )
library( "ggrepel" )

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

################################################################################
# END FUNCTIONS
################################################################################

# load vocab, tcm, glove, wv_main, wv_context
load( "out.01.emb/saved_glove.RData" )

# output dir
outdir <- 'out.02.process'
dir.create( outdir, showWarnings = FALSE )


# TODO plot matrix
image( wv_main )
image( wv_context )
image( embedding )

# neighbors
nn <- calculate_neighbors( embedding )

# write to file
readr::write_csv( nn, file = paste0( outdir, '/nearest_neighbors.csv' ), quote = 'all' )

# loop over types
type <- 'I'

for( type in c( 'I', 'II', 'III' ) )
{
    identifier <- paste0( "_", type, "$" )
    
    # only select products 'I' (i.e., 40 x 50 matrix)
    embedding_small <- embedding[ grep( identifier, rownames( embedding ) ), ]
    
    # clean names from suffix
    rownames( embedding_small ) <- gsub( identifier, "", rownames( embedding_small ) )

    # dimension reduction
    glove_umap <- umap( embedding_small, n_components = 2, metric = 'cosine', min_dist = 0.2 )
                    #n_components = 2, metric = "cosine", 
                    #n_neighbors = 5, min_dist = 0.1, spread = 15 )

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


sname <- 'average'

# get max value to normalize matrix
max_value <- max( abs( embedding_small ) )

# TODO: sort
#tmp <- as.data.frame( embedding_small )#$group <- 'piet'
#tmp$group <- rownames( tmp )
# df[order(df[,1],df[,2],decreasing=TRUE),]

# save to disk
ggsave( plot = p, dpi = 300, height = 7, width = 7, file = paste0( outdir, '/map_embedding_type_', type, '.png' ) )


# normalize and transpose
mat <- t( embedding_small / max_value )


####### matrices #######

library( "ggcorrplot" )



p_cor <- ggcorrplot( mat, outline.col = "white", ggtheme = ggplot2::theme_gray,
            colors = c("#6D9EC1", "white", "#E46726"), show.legend = FALSE, lab_size = 10 )

# save to disk
ggsave( plot = p_cor, dpi = 300, height = 11, width = 11, file = paste0( outdir, '/mat_', sname, '__type_', type, '.png' ) )
