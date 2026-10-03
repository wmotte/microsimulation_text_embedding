#!/usr/bin/env Rscript
#
# Quantitative recovery of the supermarket network from the GloVe embedding (PLOS ONE revision R1).
#
# - Spearman correlation between embedding (cosine) distance and graph (shortest-path) distance, all pairs
# - Mantel permutation test
# - Neighbour retrieval: are the nearest products in the embedding the neighbouring shelves in the network?
# - The same for the 2D UMAP maps (plus trustworthiness and Procrustes fit to the Fig 1 layout)
# - Role of product frequency
#
# Input:  out.01.emb/saved_glove.RData, out.02.process/umap_layout_type_*.tsv
# Output: out.07.recovery/*
#
################################################################################
library( 'ggplot2' )
library( 'ggrepel' )
library( 'patchwork' )
source( 'functions.R' )

set.seed( 777 )

# output dir
outdir <- 'out.07.recovery'
dir.create( outdir, showWarnings = FALSE )

# load embedding, tcm, vocab
load( 'out.01.emb/saved_glove.RData' )

sh <- get_shelves()
dg <- get_graph_distance()

# section colours (as in 02__process.R)
custom_colors <- c( beverages = '#fdae61', extras = '#fee08b', fruits = '#a6d96a', vegetables = '#d9ef8b' )

################################################################################
# 1. Embedding distance vs graph distance
################################################################################

res <- NULL
prec_items <- NULL

for( type in c( 'I', 'II', 'III', 'all' ) )
{
    if( type == 'all' ) {
        x <- embedding
    } else {
        x <- embedding[ grep( paste0( '_', type, '$' ), rownames( embedding ) ), ]
    }

    m <- recovery_metrics( x = x, dg = dg, mantel = TRUE )
    res <- rbind( res, data.frame( representation = 'GloVe embedding (50D)', set = type, m ) )

    # per-item neighbour precision
    d_hat <- 1 - cosine_matrix( x )
    d_true <- expand_graph_distance( rownames( d_hat ), dg )
    p <- neighbour_precision( d_hat, d_true )
    prec_items <- rbind( prec_items, data.frame( set = type, item = names( p ), precision = round( p, 3 ) ) )
}

################################################################################
# 2. UMAP maps vs graph distance and vs the Fig 1 layout
################################################################################

fig1 <- get_fig1_layout()
umap_maps <- list()

for( type in c( 'I', 'II', 'III' ) )
{
    lay <- as.matrix( read.table( paste0( 'out.02.process/umap_layout_type_', type, '.tsv' ), row.names = 1 ) )
    rownames( lay ) <- paste0( rownames( lay ), '_', type )

    d_low <- as.matrix( dist( lay ) )
    m <- recovery_metrics( d_hat = d_low, dg = dg, mantel = TRUE )

    # trustworthiness w.r.t. graph distance (k = 5)
    d_true <- expand_graph_distance( rownames( d_low ), dg )
    m$trustworthiness_k5 <- trustworthiness( d_low, d_true, k = 5 )

    # Procrustes fit to the Fig 1 drawing (rotation, reflection, scaling)
    labs <- sub( '_(I|II|III)$', '', rownames( lay ) )
    target <- fig1[ sh$shelf[ match( labs, sh$label ) ], ]
    pt <- vegan::protest( X = target, Y = lay, permutations = 999 )
    m$procrustes_r <- pt$t0
    m$procrustes_p <- pt$signif

    res <- rbind( res, data.frame( representation = 'UMAP map (2D)', set = type, m[ , setdiff( names( m ), c( 'trustworthiness_k5', 'procrustes_r', 'procrustes_p' ) ) ] ) )
    if( type == 'I' ) umap_extra <- m[ , c( 'trustworthiness_k5', 'procrustes_r', 'procrustes_p' ) ]

    umap_maps[[ type ]] <- list( lay = lay, pt = pt, m = m, labs = labs )

    write.table( data.frame( set = type, m ), paste0( outdir, '/umap_metrics_type_', type, '.tsv' ), sep = '\t', quote = FALSE, row.names = FALSE )
}

res$spearman_rho <- round( res$spearman_rho, 3 )
res$neighbour_precision <- round( res$neighbour_precision, 3 )
res$neighbour_precision_chance <- round( res$neighbour_precision_chance, 3 )
print( res )
readr::write_tsv( res, paste0( outdir, '/recovery_metrics.tsv' ) )
readr::write_tsv( prec_items, paste0( outdir, '/neighbour_precision_per_item.tsv' ) )

################################################################################
# 3. Fig 8: ground truth (Fig 1 layout) next to the Procrustes-aligned UMAP map (type I)
################################################################################

g <- get_shelf_graph()
el <- igraph::as_edgelist( g )
el_lab <- cbind( sh$label[ match( el[ , 1 ], sh$shelf ) ], sh$label[ match( el[ , 2 ], sh$shelf ) ] )

make_panel <- function( xy, title, subtitle )
{
    nodes <- data.frame( label = rownames( xy ), x = xy[ , 1 ], y = xy[ , 2 ] )
    nodes$section <- sh$section[ match( nodes$label, sh$label ) ]
    edges <- data.frame( x = xy[ el_lab[ , 1 ], 1 ], y = xy[ el_lab[ , 1 ], 2 ],
                         xend = xy[ el_lab[ , 2 ], 1 ], yend = xy[ el_lab[ , 2 ], 2 ] )
    ggplot() +
        geom_segment( data = edges, aes( x = x, y = y, xend = xend, yend = yend ), colour = 'gray60', linewidth = 0.4 ) +
        geom_point( data = nodes, aes( x = x, y = y, fill = section ), shape = 21, colour = 'gray20', size = 3.5 ) +
        geom_text_repel( data = nodes, aes( x = x, y = y, label = label ), size = 2.6, colour = 'gray10',
                         segment.colour = 'gray70', max.overlaps = 50, seed = 1 ) +
        scale_fill_manual( values = custom_colors ) +
        coord_equal() +
        labs( title = title, subtitle = subtitle, fill = 'department' ) +
        theme_void( base_size = 11 ) +
        theme( legend.position = 'bottom', plot.title = element_text( face = 'bold' ) )
}

# ground truth, scaled to unit size
xy_true <- fig1[ sh$shelf, ]
rownames( xy_true ) <- sh$label
xy_true <- scale( xy_true, scale = FALSE ) / sqrt( sum( scale( xy_true, scale = FALSE )^2 ) )

# aligned UMAP (type I), same scaling
um <- umap_maps[[ 'I' ]]
xy_umap <- fitted( um$pt )
rownames( xy_umap ) <- um$labs
xy_umap <- xy_umap[ sh$label, ]
xy_umap <- scale( xy_umap, scale = FALSE ) / sqrt( sum( scale( xy_umap, scale = FALSE )^2 ) )

mI <- res[ res$representation == 'UMAP map (2D)' & res$set == 'I', ]
p_true <- make_panel( xy_true, 'A  Supermarket network (ground truth)', 'Shelves and walkable connections as in Fig 1' )
p_umap <- make_panel( xy_umap, 'B  UMAP map of the embedding (type I)',
                      'Rotated and scaled to A; grey lines are the true connections' )

p8 <- ( p_true | p_umap ) + plot_layout( guides = 'collect' ) & theme( legend.position = 'bottom' )
ggsave( plot = p8, dpi = 600, height = 5.2, width = 10.5, bg = 'white', file = paste0( outdir, '/fig8_truth_vs_umap.png' ) )

################################################################################
# 4. Cosine similarity as a function of graph distance (all pairs of different products)
################################################################################

S <- cosine_matrix( embedding )
D <- expand_graph_distance( rownames( S ), dg )
ut <- upper.tri( S )
pairs <- data.frame( i = rownames( S )[ row( S )[ ut ] ], j = colnames( S )[ col( S )[ ut ] ],
                     cosine = S[ ut ], graph_distance = D[ ut ] )
pairs$same_product <- sub( '_(I|II|III)$', '', pairs$i ) == sub( '_(I|II|III)$', '', pairs$j )

# frequencies
freq <- setNames( vocab$term_count, vocab$term )
pairs$f_i <- freq[ pairs$i ]
pairs$f_j <- freq[ pairs$j ]
pairs$log_f_min <- log10( pmin( pairs$f_i, pairs$f_j ) )
pairs$log_f_max <- log10( pmax( pairs$f_i, pairs$f_j ) )

readr::write_tsv( pairs, gzfile( paste0( outdir, '/all_pairs_cosine_graph_distance.tsv.gz' ) ) )

# summary by distance
by_dist <- aggregate( cosine ~ graph_distance, data = pairs, FUN = function( v ) c( n = length( v ), median = median( v ), q25 = quantile( v, 0.25 ), q75 = quantile( v, 0.75 ), min = min( v ), max = max( v ) ) )
by_dist <- do.call( data.frame, by_dist )
readr::write_tsv( by_dist, paste0( outdir, '/cosine_by_graph_distance.tsv' ) )

p_sim <- ggplot( pairs, aes( x = factor( graph_distance ), y = cosine ) ) +
    geom_boxplot( outlier.size = 0.4, fill = 'gray90', colour = 'gray30' ) +
    xlab( 'Distance in the supermarket network (number of steps)' ) +
    ylab( 'Cosine similarity of the two products' ) +
    theme_bw( base_size = 12 ) + theme( panel.grid.minor = element_blank() )
ggsave( plot = p_sim, dpi = 600, height = 4.5, width = 7, bg = 'white', file = paste0( outdir, '/cosine_by_graph_distance.png' ) )

################################################################################
# 5. Frequency: does product frequency affect similarity beyond graph distance?
################################################################################

dp <- pairs[ !pairs$same_product, ]

# rank-based partial correlation of x and y given z
partial_spearman <- function( x, y, z )
{
    rx <- resid( lm( rank( x ) ~ rank( z ) ) )
    ry <- resid( lm( rank( y ) ~ rank( z ) ) )
    cor( rx, ry )
}

fr <- data.frame(
    n_pairs = nrow( dp ),
    rho_cos_distance = cor( dp$cosine, dp$graph_distance, method = 'spearman' ),
    rho_cos_logfmin = cor( dp$cosine, dp$log_f_min, method = 'spearman' ),
    partial_rho_cos_distance_given_fmin = partial_spearman( dp$cosine, dp$graph_distance, dp$log_f_min ),
    partial_rho_cos_fmin_given_distance = partial_spearman( dp$cosine, dp$log_f_min, dp$graph_distance ),
    rho_distance_logfmin = cor( dp$graph_distance, dp$log_f_min, method = 'spearman' ) )

# linear model on standardised variables (graph distance as numeric)
fit <- lm( scale( cosine ) ~ scale( graph_distance ) + scale( log_f_min ) + scale( log_f_max ), data = dp )
co <- summary( fit )$coefficients
fr$beta_distance <- co[ 2, 1 ]
fr$beta_log_f_min <- co[ 3, 1 ]
fr$beta_log_f_max <- co[ 4, 1 ]
fr$r2_full <- summary( fit )$r.squared
fr$r2_distance_only <- summary( lm( cosine ~ graph_distance, data = dp ) )$r.squared
fr <- data.frame( lapply( fr, function( v ) if( is.numeric( v ) ) signif( v, 3 ) else v ) )
print( t( fr ) )
readr::write_tsv( fr, paste0( outdir, '/frequency_effects.tsv' ) )

# product frequency (summed over variants) vs shelf number: does the numbering create a frequency gradient?
pf <- aggregate( freq, by = list( label = sub( '_(I|II|III)$', '', names( freq ) ) ), FUN = sum )
colnames( pf ) <- c( 'label', 'frequency' )
pf$shelf_number <- as.numeric( sub( 'H', '', sh$shelf[ match( pf$label, sh$label ) ] ) )
pf$section <- sh$section[ match( pf$label, sh$label ) ]
pf <- pf[ order( pf$shelf_number ), ]
readr::write_tsv( pf, paste0( outdir, '/product_frequency_by_shelf.tsv' ) )
readr::write_tsv( data.frame( rho_frequency_shelf_number = signif( cor( pf$frequency, pf$shelf_number, method = 'spearman' ), 3 ),
                              min_frequency = min( freq ), max_frequency = max( freq ), median_frequency = median( freq ) ),
                  paste0( outdir, '/frequency_gradient.tsv' ) )

p_freq <- ggplot( pf, aes( x = shelf_number, y = frequency, fill = section ) ) +
    geom_col( colour = 'gray30', linewidth = 0.2 ) +
    scale_fill_manual( values = custom_colors ) +
    xlab( 'Shelf number (H1-H40)' ) + ylab( 'Times selected (three variants summed)' ) +
    theme_bw( base_size = 12 ) + theme( legend.position = 'top', panel.grid.minor = element_blank() )
ggsave( plot = p_freq, dpi = 600, height = 4.5, width = 8, bg = 'white', file = paste0( outdir, '/product_frequency_by_shelf.png' ) )

################################################################################
# 6. Vegetables and beverages: graph structure behind the merged cluster
################################################################################

# mean graph distance between departments (shelf level)
secs <- unique( sh$section )
md <- outer( secs, secs, Vectorize( function( a, b ) mean( dg[ sh$label[ sh$section == a ], sh$label[ sh$section == b ] ] ) ) )
dimnames( md ) <- list( secs, secs )

# number of edges between departments
cross <- table( factor( sh$section[ match( el_lab[ , 1 ], sh$label ) ], levels = secs ),
                factor( sh$section[ match( el_lab[ , 2 ], sh$label ) ], levels = secs ) )
cross <- cross + t( cross ) - diag( diag( cross ) )

sink( paste0( outdir, '/department_structure.txt' ) )
cat( 'Mean shortest-path distance between departments (shelves)\n' )
print( round( md, 2 ) )
cat( '\nNumber of walkable connections between departments\n' )
print( cross )
sink()

writeLines( capture.output( sessioninfo::session_info() ), paste0( outdir, '/session_info.txt' ) )
