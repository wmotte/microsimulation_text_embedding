#!/usr/bin/env Rscript
#
# Cluster similarity matrices
#
################################################################################
library( "ggplot2" )

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
outdir <- 'out.06.cluster'
dir.create( outdir, showWarnings = FALSE )

# read: https://cran.r-project.org/web/packages/kmed/vignettes/kmedoid.html


type <- 'I'


# loop over other types to plot matrix
for( type in c( 'I', 'II', 'III' ) )
{

    # load 'simil' df
    infile <- paste0( 'out.03.superheat/wv_embedding_type_', type, '__similarity_matrix.Rdata' )
    load( infile )
    m <- as.matrix( simil )
    
    # rescale similarity between 1 and 100%
    m <- scales::rescale( m, c( 1, 100 ) )
    
    # convert to distance matrix
    d <- as.data.frame( 1 / m )
    
    ######## Clustering ########
    
    library( 'cluster' )
    
    # Calculate silhouette width for many k using PAM
    sil_width <- c( NA )

    # cluster k-medoids         
    pam_fit <- cluster::pam( d, diss = TRUE, k = 4 )

    # get the 4 gold standard cluster sections
    sec <- get_sections()
    
    # combine with gold standard
    tmp <- data.frame( cluster = pam_fit$clustering, label = names( pam_fit$clustering ) )
    rownames( tmp ) <- NULL
    comb <- merge( sec, tmp )
    
    # get crosstab [with medoids on the rows and gold-standards as cols]
    tt <- table( comb$cluster, comb$section )
    rownames( tt ) <- pam_fit$medoids
    
    # write to disk
    outfile_tt <- paste0( outdir, '/clustering_table_', type, '.tsv' )
    write.table( tt, outfile_tt, quote = TRUE )
    
    # write silhouette widths
    outfile_med <- paste0( outdir, '/clustering_details_', type, '.tsv' )
    xx <- data.frame( pam_fit$silinfo$widths )
    xx$label <- rownames( xx )
    rownames( xx ) <- NULL
    xx$neighbor <- NULL
    xx$sil_width <- round( xx$sil_width, 4 )
    write.table( xx[ , c( 'label', 'cluster', 'sil_width' ) ], outfile_med, quote = TRUE )
    
    # save to disk
    outfile <- paste0( outdir, '/distance_', type, '.png' )
    png( outfile, height = 5000, width = 5000, res = 300 )
    
    # heat plot
    superheat( d,
           
           heat.na.col = "white",
           
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
