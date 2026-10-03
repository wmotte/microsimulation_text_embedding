#!/usr/bin/env Rscript
#
# Summaries (mean and SD over 10 seeds) and Fig 10 from the output of 08__robustness.R.
#
# Output: out.08.robustness/robustness_summary.tsv, umap_summary.tsv, corpus_info_summary.tsv, fig10_robustness.png
#
################################################################################
library( 'ggplot2' )
library( 'patchwork' )

outdir <- 'out.08.robustness'
res <- as.data.frame( readr::read_tsv( paste0( outdir, '/robustness_long.tsv' ), show_col_types = FALSE ) )
umap_res <- as.data.frame( readr::read_tsv( paste0( outdir, '/umap_long.tsv' ), show_col_types = FALSE ) )
corpus_info <- as.data.frame( readr::read_tsv( paste0( outdir, '/corpus_info.tsv' ), show_col_types = FALSE ) )

################################################################################
# Summaries (mean and SD over seeds)
################################################################################

msd <- function( v ) sprintf( '%.3f (%.3f)', mean( v ), sd( v ) )

agg <- res
for( v in c( 'window', 'rank', 'n_lists' ) ) agg[[ v ]][ is.na( agg[[ v ]] ) ] <- -1
summ <- aggregate( cbind( spearman_rho, neighbour_precision ) ~ condition + representation + set + window + rank + n_lists,
                   data = agg, FUN = function( v ) c( mean = mean( v ), sd = sd( v ) ) )
summ <- do.call( data.frame, summ )
summ[ summ == -1 ] <- NA
num <- sapply( summ, is.numeric )
summ[ num ] <- lapply( summ[ num ], round, 3 )
readr::write_tsv( summ, paste0( outdir, '/robustness_summary.tsv' ) )

chance <- unique( res$neighbour_precision_chance[ res$set == 'all' ] )

um <- data.frame( spearman_rho = msd( umap_res$spearman_rho ), neighbour_precision = msd( umap_res$neighbour_precision ),
                  trustworthiness_k5 = msd( umap_res$trustworthiness_k5 ) )
readr::write_tsv( um, paste0( outdir, '/umap_summary.tsv' ) )

ci <- data.frame( n_simulated = msd( corpus_info$n_simulated ), fraction_inside_at_end = msd( corpus_info$fraction_inside_at_end ),
                  tokens = msd( corpus_info$tokens ) )
readr::write_tsv( ci, paste0( outdir, '/corpus_info_summary.tsv' ) )

################################################################################
# Figure (Fig 10): baselines/controls and sweeps, all 120 products
################################################################################

ra <- res[ res$set == 'all', ]
long <- rbind( data.frame( ra, metric = 'Spearman rho (embedding vs graph distance)', value = ra$spearman_rho ),
               data.frame( ra, metric = 'Neighbour retrieval (precision)', value = ra$neighbour_precision ) )
long$metric <- factor( long$metric, levels = c( 'Spearman rho (embedding vs graph distance)', 'Neighbour retrieval (precision)' ) )

base <- long[ long$condition %in% c( 'main', 'control' ), ]
lev <- c( 'GloVe (rank 50)', 'word2vec skip-gram (dim 50)', 'SVD of PPMI (rank 50)', 'PPMI', 'Raw co-occurrence counts',
          'GloVe, tokens shuffled within lists', 'GloVe, tokens shuffled across corpus' )
base$representation <- factor( base$representation, levels = rev( lev ) )

pA <- ggplot( base, aes( x = value, y = representation ) ) +
    geom_jitter( height = 0.15, width = 0, size = 1, alpha = 0.6, colour = 'gray30' ) +
    stat_summary( fun = mean, geom = 'point', shape = 23, size = 2.8, fill = 'orange' ) +
    facet_wrap( ~metric, scales = 'free_x' ) +
    geom_vline( data = data.frame( metric = factor( 'Neighbour retrieval (precision)', levels = levels( long$metric ) ), v = chance ), aes( xintercept = v ), linetype = 2, colour = 'gray50' ) +
    xlab( NULL ) + ylab( NULL ) + ggtitle( 'A  Representations and controls (10 seeds)' ) +
    theme_bw( base_size = 10 ) + theme( plot.title = element_text( face = 'bold' ), panel.grid.minor = element_blank() )

mean_sd <- function( v ) data.frame( y = mean( v ), ymin = mean( v ) - sd( v ), ymax = mean( v ) + sd( v ) )

sweep_plot <- function( cond, xvar, xlab, title, logx = FALSE )
{
    d <- long[ long$condition == cond | ( long$condition == 'main' & long$representation == 'GloVe (rank 50)' ), ]
    d <- d[ !is.na( d[[ xvar ]] ), ]
    p <- ggplot( d, aes( x = .data[[ xvar ]], y = value ) ) +
        stat_summary( fun.data = mean_sd, geom = 'pointrange', size = 0.3 ) +
        stat_summary( fun = mean, geom = 'line' ) +
        facet_wrap( ~metric, scales = 'free_y', ncol = 2 ) +
        xlab( xlab ) + ylab( NULL ) + ggtitle( title ) +
        theme_bw( base_size = 10 ) + theme( plot.title = element_text( face = 'bold' ), panel.grid.minor = element_blank() )
    if( logx ) p <- p + scale_x_log10()
    p
}

pB <- sweep_plot( 'size', 'n_lists', 'Number of shopping lists (log scale)', 'B  Corpus size', logx = TRUE )
pC <- sweep_plot( 'window', 'window', 'Context window (words)', 'C  Context window' )
pD <- sweep_plot( 'rank', 'rank', 'Embedding dimension (log scale)', 'D  Embedding dimension', logx = TRUE )

p10 <- pA / pB / pC / pD + plot_layout( heights = c( 1.3, 1, 1, 1 ) )
ggsave( plot = p10, dpi = 300, height = 8.75, width = 7.5, bg = 'white', file = paste0( outdir, '/fig10_robustness.png' ) )

