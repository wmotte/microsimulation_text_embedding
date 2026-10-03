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

source( 'functions.R' )
library( 'cluster' )

# output dir
outdir <- 'out.06.cluster'
dir.create( outdir, showWarnings = FALSE )

# read: https://cran.r-project.org/web/packages/kmed/vignettes/kmedoid.html

# load embedding (120 x 50)
load( 'out.01.emb/saved_glove.RData' )

# departments (ground truth labels) and graph distance between shelves
sh <- get_shelves()
dg <- get_graph_distance()

###
# PAM on cosine distance (1 - cosine similarity), with silhouette for k = 2..8
##
run_pam <- function( x, k_range = 2:8, k_fixed = 4 )
{
    d <- as.dist( 1 - cosine_matrix( x ) )
    sil <- sapply( k_range, function( k ) cluster::pam( d, diss = TRUE, k = k )$silinfo$avg.width )
    fit <- cluster::pam( d, diss = TRUE, k = k_fixed )
    list( d = d, sil = data.frame( k = k_range, avg_silhouette = round( sil, 4 ) ), fit = fit )
}

# reference partition of the network itself: PAM on shortest-path distance between shelves (k = 4)
graph_fit <- cluster::pam( as.dist( dg ), diss = TRUE, k = 4 )
graph_cluster <- graph_fit$clustering

summary_all <- NULL
type <- 'I'

# per variant (40 products) and all 120 products together
for( type in c( 'I', 'II', 'III', 'all' ) )
{
    if( type == 'all' ) {
        x <- embedding
    } else {
        x <- embedding[ grep( paste0( '_', type, '$' ), rownames( embedding ) ), ]
    }
    labs <- sub( '_(I|II|III)$', '', rownames( x ) )
    section <- sh$section[ match( labs, sh$label ) ]
    
    res <- run_pam( x )
    fit <- res$fit
    
    # crosstab [medoids on the rows and departments as cols]
    tt <- table( fit$clustering, section )
    rownames( tt ) <- fit$medoids
    write.table( tt, paste0( outdir, '/clustering_table_', type, '.tsv' ), quote = TRUE )
    
    # medoids (Table 6 uses the medoids of the clustering of all 120 products)
    if( type == 'all' )
        writeLines( fit$medoids, paste0( outdir, '/medoids_all_k4.txt' ) )
    
    # silhouette widths per item
    xx <- data.frame( fit$silinfo$widths )
    xx$label <- rownames( xx )
    rownames( xx ) <- NULL
    xx$neighbor <- NULL
    xx$sil_width <- round( xx$sil_width, 4 )
    write.table( xx[ , c( 'label', 'cluster', 'sil_width' ) ], paste0( outdir, '/clustering_details_', type, '.tsv' ), quote = TRUE )
    
    # silhouette over k
    write.table( res$sil, paste0( outdir, '/silhouette_by_k_', type, '.tsv' ), quote = FALSE, sep = '\t', row.names = FALSE )
    
    # agreement with departments and with the graph partition
    gclust <- graph_cluster[ match( labs, rownames( dg ) ) ]
    summary_all <- rbind( summary_all, data.frame( 
        set = type, n_items = nrow( x ),
        best_k_silhouette = res$sil$k[ which.max( res$sil$avg_silhouette ) ],
        avg_silhouette_k4 = round( fit$silinfo$avg.width, 3 ),
        ari_departments = round( mclust::adjustedRandIndex( fit$clustering, section ), 3 ),
        ari_graph_partition = round( mclust::adjustedRandIndex( fit$clustering, gclust ), 3 ),
        purity_departments = round( sum( apply( table( fit$clustering, section ), 1, max ) ) / nrow( x ), 3 ),
        medoids = paste( fit$medoids, collapse = ' ' ) ) )
    
    # distance plot (type I, II, III)
    if( type != 'all' )
    {
        d <- as.matrix( res$d )
        rownames( d ) <- colnames( d ) <- labs
        outfile <- paste0( outdir, '/distance_', type, '.png' )
        png( outfile, height = 5000, width = 5000, res = 300 )
        superheat( d,
               heat.na.col = "white",
               grid.hline.col = "gray30",
               grid.vline.col = "gray30",
               bottom.label.text.angle = -90,
               left.label.text.size = 5,
               bottom.label.text.size = 5,
               legend = FALSE )
        dev.off()
    }
}

# graph partition and how it relates to departments
gp <- data.frame( label = names( graph_cluster ), graph_cluster = graph_cluster, 
                  section = sh$section[ match( names( graph_cluster ), sh$label ) ] )
write.table( table( gp$graph_cluster, gp$section ), paste0( outdir, '/graph_partition_vs_departments.tsv' ), quote = TRUE )
summary_all$ari_graph_partition_vs_departments <- round( mclust::adjustedRandIndex( gp$graph_cluster, gp$section ), 3 )

print( summary_all )
readr::write_tsv( summary_all, paste0( outdir, '/clustering_summary.tsv' ) )
