#!/usr/bin/env Rscript
#
# Robustness, baselines and controls (PLOS ONE revision R1).
#
# For 10 independent seeds (new corpus + new model fit for each seed):
#   - GloVe (as in 01__emb.R)
#   - baselines on the same co-occurrence matrix: raw weighted counts, PPMI, SVD of PPMI (rank 50)
#   - word2vec skip-gram (dimension 50, window 5)
#   - controls: tokens shuffled across the whole corpus (destroys all structure),
#               tokens shuffled within each list (keeps list membership, destroys order)
# Sweeps (GloVe, 10 seeds each): corpus size, context window, embedding dimension.
#
# Corpora are generated with simulate_corpus_fast() (functions.R): same model as 00__sim.R, vectorised.
#
# Output: out.08.robustness/robustness_long.tsv, umap_long.tsv, corpus_info.tsv
# Summaries and Fig 10 are made by 09__fig_robustness.R
#
################################################################################
source( 'functions.R' )

# output dir
outdir <- 'out.08.robustness'
dir.create( outdir, showWarnings = FALSE )

seeds <- 1:10
dg <- get_graph_distance()

###
# Metrics for one representation: all 120 products and the 40 type I products
##
score <- function( x, seed, condition, representation, extra = list() )
{
    lr_attr <- attr( x, 'learning_rate' )
    x <- x[ order( rownames( x ) ), , drop = FALSE ]
    attr( x, 'learning_rate' ) <- lr_attr
    xI <- x[ grep( '_I$', rownames( x ) ), , drop = FALSE ]
    out <- rbind( data.frame( set = 'all', recovery_metrics( x = x, dg = dg ) ),
                  data.frame( set = 'I', recovery_metrics( x = xI, dg = dg ) ) )
    lr <- attr( x, 'learning_rate' )
    out <- data.frame( seed = seed, condition = condition, representation = representation, out,
                       glove_learning_rate = ifelse( is.null( lr ), NA, lr ) )
    sweep <- list( window = NA, rank = NA, n_lists = NA )
    sweep[ names( extra ) ] <- extra
    for( nm in names( sweep ) ) out[[ nm ]] <- sweep[[ nm ]]
    return( out )
}

###
# UMAP of the type I products (as in 02__process.R) and its agreement with the graph
##
umap_score <- function( x, seed )
{
    xI <- x[ grep( '_I$', rownames( x ) ), , drop = FALSE ]
    lay <- umap::umap( xI, n_components = 2, metric = 'cosine', min_dist = 0.2, n_neighbors = 15, random_state = seed )$layout
    rownames( lay ) <- rownames( xI )
    d_low <- as.matrix( dist( lay ) )
    m <- recovery_metrics( d_hat = d_low, dg = dg )
    m$trustworthiness_k5 <- trustworthiness( d_low, expand_graph_distance( rownames( d_low ), dg ), k = 5 )
    data.frame( seed = seed, m )
}

shuffle_tokens <- function( lists, within = FALSE )
{
    toks <- strsplit( lists, ' ' )
    if( within ) {
        toks <- lapply( toks, function( v ) v[ sample.int( length( v ) ) ] )
    } else {
        len <- lengths( toks )
        all_t <- sample( unlist( toks ) )
        toks <- split( all_t, rep( seq_along( len ), len ) )
    }
    vapply( toks, paste, character( 1 ), collapse = ' ' )
}

res <- NULL
umap_res <- NULL
corpus_info <- NULL

for( s in seeds )
{
    cat( '\n==== seed', s, format( Sys.time(), '%H:%M:%S' ), '====\n' )

    lists <- simulate_corpus_fast( seed = 1000 + s, n_lists = 80000 )
    corpus_info <- rbind( corpus_info, data.frame( seed = s, n_simulated = attr( lists, 'n_simulated' ),
                                                   fraction_inside_at_end = attr( lists, 'fraction_inside_at_end' ),
                                                   tokens = sum( stringr::str_count( lists, ' ' ) + 1 ) ) )

    # main: GloVe, window 5, rank 50
    tc <- make_tcm( lists, window = 5 )
    X <- symmetric_tcm( tc$tcm )
    emb <- fit_glove( tc$tcm, rank = 50, seed = s )
    res <- rbind( res, score( emb, s, 'main', 'GloVe (rank 50)' ) )
    umap_res <- rbind( umap_res, umap_score( emb, s ) )

    # baselines on the same co-occurrence matrix
    res <- rbind( res, score( X, s, 'main', 'Raw co-occurrence counts' ) )
    res <- rbind( res, score( ppmi( X ), s, 'main', 'PPMI' ) )
    res <- rbind( res, score( svd_ppmi( X, k = 50 ), s, 'main', 'SVD of PPMI (rank 50)' ) )
    res <- rbind( res, score( fit_word2vec( lists, dim = 50, window = 5, seed = s ), s, 'main', 'word2vec skip-gram (dim 50)' ) )

    # controls
    set.seed( 2000 + s )
    tc_sh <- make_tcm( shuffle_tokens( lists, within = FALSE ), window = 5 )
    res <- rbind( res, score( fit_glove( tc_sh$tcm, rank = 50, seed = s ), s, 'control', 'GloVe, tokens shuffled across corpus' ) )
    tc_sw <- make_tcm( shuffle_tokens( lists, within = TRUE ), window = 5 )
    res <- rbind( res, score( fit_glove( tc_sw$tcm, rank = 50, seed = s ), s, 'control', 'GloVe, tokens shuffled within lists' ) )

    # sweep: context window
    for( w in c( 1, 2, 10 ) )
    {
        tcw <- make_tcm( lists, window = w )
        res <- rbind( res, score( fit_glove( tcw$tcm, rank = 50, seed = s ), s, 'window', 'GloVe (rank 50)', list( window = w ) ) )
    }

    # sweep: embedding dimension
    for( r in c( 2, 5, 10, 25, 100 ) )
        res <- rbind( res, score( fit_glove( tc$tcm, rank = r, seed = s, relax = TRUE ), s, 'rank', paste0( 'GloVe (rank ', r, ')' ), list( rank = r ) ) )

    # sweep: corpus size (first n lists)
    for( n in c( 1000, 5000, 10000, 20000, 40000 ) )
    {
        tcn <- make_tcm( lists[ 1:n ], window = 5 )
        res <- rbind( res, score( fit_glove( tcn$tcm, rank = 50, seed = s ), s, 'size', 'GloVe (rank 50)', list( n_lists = n ) ) )
    }

    # write intermediate results
    readr::write_tsv( res, paste0( outdir, '/robustness_long.tsv' ) )
}

# fill sweep columns for the main condition (window 5, rank 50, 80,000 lists)
main_glove <- res$condition == 'main' & res$representation == 'GloVe (rank 50)'
res$window[ main_glove ] <- 5
res$rank[ main_glove ] <- 50
res$n_lists[ main_glove ] <- 80000

readr::write_tsv( res, paste0( outdir, '/robustness_long.tsv' ) )
readr::write_tsv( umap_res, paste0( outdir, '/umap_long.tsv' ) )
readr::write_tsv( corpus_info, paste0( outdir, '/corpus_info.tsv' ) )

writeLines( capture.output( sessioninfo::session_info() ), paste0( outdir, '/session_info.txt' ) )
cat( '\nDone', format( Sys.time(), '%H:%M:%S' ), '\n' )
