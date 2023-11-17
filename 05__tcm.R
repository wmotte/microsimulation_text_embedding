#!/usr/bin/env Rscript
#
# Visualize tcm
#
################################################################################
library( "ggplot2" )
#library( "ggrepel" )
#library( "ggcorrplot" )

library( 'superheat' )

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



################################################################################

# output dir
outdir <- 'out.05.tcm'
dir.create( outdir, showWarnings = FALSE )

# load tcm
load( 'out.01.emb/saved_glove.RData' )
embedding <- glove <- vocab <- wv_context <- wv_main <- NULL


for( type in c( 'I', 'II', 'III' ) )
{
    
    identifier <- paste0( '_', type, '$' )

    mat <- tcm[ grep( identifier, rownames( tcm ) ), grep( identifier, colnames( tcm ) ) ]
    rownames( mat ) <- gsub( identifier, '', rownames( mat ) )
    colnames( mat ) <- gsub( identifier, '', colnames( mat ) )
    
    # get separate matrix for numbers [because of NA text issues]
    mmat <- round( t( as.matrix( mat ) ), 0 )
    mmat_text <- mmat
    mmat[ mmat == 0 ] <- NA

    # save to disk
    outfile <- paste0( outdir, '/tcm_', type, '.png' )
    png( outfile, height = 5000, width = 5000, res = 300 )
    
    # heat plot
    superheat( mmat,
           # place dendrograms on columns and rows 
           #row.dendrogram = F, col.dendrogram = F,
           
           heat.na.col = "white",
           
           X.text = mmat_text,
           X.text.size = 3.5,
           X.text.col = 'white',
           X.text.angle = 45,
           # make gridlines white for enhanced prettiness
           grid.hline.col = "gray30",
           grid.vline.col = "gray30",
           
           # rotate bottom label text
           bottom.label.text.angle = -90,
           
           left.label.text.size = 5,
           bottom.label.text.size = 5,
           legend = FALSE )

    dev.off()
}
