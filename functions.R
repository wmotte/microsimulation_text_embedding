#!/usr/bin/env Rscript
#
# Shared functions for the recovery and robustness analyses (07__recovery.R, 08__robustness.R)
# and the clustering (06__cluster.R).
#
# Added in the PLOS ONE revision (R1), October 2026.
#
################################################################################

###
# Shelf -> product -> section (department) table [H1 ... H40]
##
get_shelves <- function()
{
    data.frame(
        shelf = paste0( 'H', 1:40 ),
        label = c( 'Apples', 'Bananas', 'Grapefruit', 'Grapes', 'Kiwi', 'Lime', 'Mangoes', 'Oranges',
                   'Pineapples', 'Strawberries', 'Watermelons',
                   'Cauliflower', 'Cucumbers', 'Eggplant', 'Onion', 'Spinach', 'Tomatoes', 'Peppers', 'Zucchini',
                   'BitterLemon', 'Cassis', 'Coffee', 'Cola', 'Lemonade', 'Sodawater', 'Tea',
                   'CaramelBars', 'ChewingGum', 'ChocolateBar', 'Crackers', 'GummyBears', 'MixedNuts', 'Popcorn',
                   'PotatoChips', 'Snickers', 'Sweets', 'Snacks', 'Twix', 'Pringles', 'ShoppingBag' ),
        section = c( rep( 'fruits', 11 ), rep( 'vegetables', 8 ), rep( 'beverages', 7 ), rep( 'extras', 14 ) ),
        stringsAsFactors = FALSE )
}

###
# Coordinates of the shelves as drawn in Fig 1 of the manuscript (pixels, y downwards).
# Used only to display the ground truth next to the UMAP map (Fig 8); all statistics use graph distances.
##
get_fig1_layout <- function()
{
    xy <- rbind(
        H1 = c( 67, 210 ), H2 = c( 248, 310 ), H3 = c( 333, 476 ), H4 = c( 406, 647 ), H5 = c( 373, 188 ),
        H6 = c( 455, 355 ), H7 = c( 528, 529 ), H8 = c( 497, 71 ), H9 = c( 579, 236 ), H10 = c( 646, 409 ),
        H11 = c( 798, 398 ), H12 = c( 886, 520 ), H13 = c( 994, 645 ), H14 = c( 1081, 752 ), H15 = c( 1161, 904 ),
        H16 = c( 1226, 1037 ), H17 = c( 994, 380 ), H18 = c( 1091, 492 ), H19 = c( 1199, 623 ), H20 = c( 1285, 772 ),
        H21 = c( 1361, 928 ), H22 = c( 1080, 218 ), H23 = c( 1196, 351 ), H24 = c( 1310, 489 ), H25 = c( 1406, 641 ),
        H26 = c( 1486, 809 ), H27 = c( 1559, 982 ), H28 = c( 1708, 919 ), H29 = c( 1817, 1058 ), H30 = c( 1726, 1214 ),
        H31 = c( 1644, 1379 ), H32 = c( 1614, 1570 ), H33 = c( 1457, 1572 ), H34 = c( 1322, 1474 ), H35 = c( 1361, 1286 ),
        H36 = c( 1465, 1138 ), H37 = c( 1734, 1697 ), H38 = c( 1835, 1552 ), H39 = c( 1639, 1850 ), H40 = c( 1848, 1830 ) )
    xy[ , 2 ] <- -xy[ , 2 ]
    colnames( xy ) <- c( 'x', 'y' )
    return( xy )
}

###
# Transition matrix (41 x 41; H41 = absorbing exit state), same weights as in 00__sim.R
##
get_transition_matrix <- function( p_backward = 2, p_forward = 10, p_diag = 1 )
{
    raw <- suppressMessages( readxl::read_xlsx( 'doc/transition_matrix.xlsx' ) )
    raw$`...1` <- NULL
    tmat <- suppressWarnings( as.matrix( as.data.frame( lapply( as.data.frame( raw ), as.numeric ) ) ) )  # empty cells -> NA

    tmat[ lower.tri( tmat ) ] <- tmat[ lower.tri( tmat ) ] * p_backward
    tmat[ upper.tri( tmat, diag = FALSE ) ] <- tmat[ upper.tri( tmat, diag = FALSE ) ] * p_forward
    diag( tmat ) <- p_diag
    pmat <- tmat / rowSums( tmat, na.rm = TRUE )
    pmat[ is.na( pmat ) ] <- 0

    # absorbing exit
    pmat[ 41, ] <- 0
    pmat[ 41, 41 ] <- 1

    dimnames( pmat ) <- list( paste0( 'H', 1:41 ), paste0( 'H', 1:41 ) )
    return( pmat )
}

###
# Undirected shelf graph (40 shelves, exit state removed)
##
get_shelf_graph <- function()
{
    pmat <- get_transition_matrix()[ 1:40, 1:40 ]
    adj <- ( pmat > 0 ) * 1
    diag( adj ) <- 0
    adj <- ( ( adj + t( adj ) ) > 0 ) * 1
    g <- igraph::graph_from_adjacency_matrix( adj, mode = 'undirected' )
    return( g )
}

###
# Shortest-path (hop) distance between shelves, named by product label (40 x 40)
##
get_graph_distance <- function()
{
    g <- get_shelf_graph()
    d <- igraph::distances( g )
    sh <- get_shelves()
    dimnames( d ) <- list( sh$label[ match( rownames( d ), sh$shelf ) ], sh$label[ match( colnames( d ), sh$shelf ) ] )
    return( d )
}

###
# Expand the 40 x 40 shelf distance to an arbitrary set of tokens 'Label_TYPE'
##
expand_graph_distance <- function( tokens, dg = get_graph_distance() )
{
    labs <- sub( '_(I|II|III)$', '', tokens )
    d <- dg[ labs, labs ]
    dimnames( d ) <- list( tokens, tokens )
    return( d )
}

###
# Cosine similarity between the rows of x
##
cosine_matrix <- function( x )
{
    xn <- x / sqrt( rowSums( x^2 ) )
    s <- tcrossprod( xn )
    s[ s > 1 ] <- 1
    s[ s < -1 ] <- -1
    return( s )
}

###
# Neighbour retrieval: for each item i with k_i true graph neighbours (distance 1),
# the fraction of its k_i nearest items (by d_hat) that are true neighbours.
# Items of the same product (other variants) are excluded from the candidate set.
##
neighbour_precision <- function( d_hat, d_true )
{
    n <- nrow( d_hat )
    labs <- sub( '_(I|II|III)$', '', rownames( d_hat ) )
    prec <- rep( NA, n )
    for( i in 1:n )
    {
        cand <- which( labs != labs[ i ] )
        truth <- cand[ d_true[ i, cand ] == 1 ]
        k <- length( truth )
        if( k == 0 ) next
        nn <- cand[ order( d_hat[ i, cand ] ) ][ 1:k ]
        prec[ i ] <- mean( nn %in% truth )
    }
    names( prec ) <- rownames( d_hat )
    return( prec )
}

###
# Chance level of neighbour_precision (expected precision of a random ranking)
##
neighbour_precision_chance <- function( d_true )
{
    labs <- sub( '_(I|II|III)$', '', rownames( d_true ) )
    ch <- sapply( 1:nrow( d_true ), function( i ) {
        cand <- which( labs != labs[ i ] )
        mean( d_true[ i, cand ] == 1 ) } )
    return( mean( ch ) )
}

###
# Trustworthiness (Venna & Kaski 2001) of a low-dimensional map with respect to the true graph distances.
# Ranks of the true distances use average ties.
##
trustworthiness <- function( d_low, d_true, k = 5 )
{
    n <- nrow( d_low )
    penalty <- 0
    for( i in 1:n )
    {
        others <- setdiff( 1:n, i )
        r_true <- rank( d_true[ i, others ], ties.method = 'average' )
        nn_low <- others[ order( d_low[ i, others ] ) ][ 1:k ]
        nn_true_rank <- r_true[ match( nn_low, others ) ]
        penalty <- penalty + sum( pmax( nn_true_rank - k, 0 ) )
    }
    return( 1 - 2 / ( n * k * ( 2 * n - 3 * k - 1 ) ) * penalty )
}

###
# All recovery metrics for one representation (rows = tokens)
#   x:      item x feature matrix (embedding or count/PPMI rows)
#   d_hat:  optional distance matrix (if given, x is ignored)
##
recovery_metrics <- function( x = NULL, d_hat = NULL, dg = get_graph_distance(), mantel = FALSE )
{
    if( is.null( d_hat ) )
        d_hat <- 1 - cosine_matrix( x )

    d_true <- expand_graph_distance( rownames( d_hat ), dg )
    ut <- upper.tri( d_hat )

    rho <- suppressWarnings( cor( d_hat[ ut ], d_true[ ut ], method = 'spearman' ) )
    prec <- neighbour_precision( d_hat, d_true )

    out <- data.frame( n_items = nrow( d_hat ), n_pairs = sum( ut ),
                       spearman_rho = rho,
                       neighbour_precision = mean( prec, na.rm = TRUE ),
                       neighbour_precision_chance = neighbour_precision_chance( d_true ) )

    if( mantel )
    {
        mt <- vegan::mantel( as.dist( d_hat ), as.dist( d_true ), method = 'spearman', permutations = 999 )
        out$mantel_p <- mt$signif
    }
    return( out )
}

###
# Fast (vectorised) simulation of shoppers.
# Same model as 00__sim.R (same transition matrix, 75 cycles, purchase probability 0.15, >= 6 products,
# three variants with equal probability, duplicates removed), but all shoppers are moved in parallel.
# The random number stream therefore differs from 00__sim.R; the distribution of the corpus is the same.
##
simulate_corpus_fast <- function( seed, n_lists = 80000, n_t = 75, prob = 0.15, min_products = 6, batch = 20000 )
{
    set.seed( seed )
    pmat <- get_transition_matrix()
    cum <- t( apply( pmat, 1, cumsum ) )
    cum[ , ncol( cum ) ] <- 1
    sh <- get_shelves()

    lists <- character( 0 )
    n_sim <- 0
    n_inside <- 0

    while( length( lists ) < n_lists )
    {
        # walks: batch x (n_t + 1) state indices; everybody starts at H1
        m <- matrix( 1L, nrow = batch, ncol = n_t + 1 )
        for( t in 1:n_t )
        {
            u <- runif( batch )
            cur <- m[ , t ]
            # next state: first column where cumulative probability exceeds u
            m[ , t + 1 ] <- 1L + rowSums( u > cum[ cur, , drop = FALSE ] )
        }
        m[ m > 41 ] <- 41L
        n_inside <- n_inside + sum( m[ , n_t + 1 ] != 41 )

        # purchases: every visited state (except the exit) with probability prob
        take <- matrix( runif( length( m ) ) < prob, nrow = batch ) & ( m != 41 )
        n_take <- rowSums( take )
        keep <- which( n_take >= min_products )

        variant <- c( 'I', 'II', 'III' )[ sample.int( 3, sum( take[ keep, ] ), replace = TRUE ) ]
        prod_idx <- t( m[ keep, , drop = FALSE ] )[ t( take[ keep, , drop = FALSE ] ) ]
        tokens <- paste0( sh$label[ prod_idx ], '_', variant )
        sentence_id <- rep( seq_along( keep ), n_take[ keep ] )
        new_lists <- vapply( split( tokens, sentence_id ), paste, character( 1 ), collapse = ' ' )

        lists <- c( lists, new_lists )
        lists <- lists[ !duplicated( lists ) ]
        n_sim <- n_sim + batch
    }

    out <- lists[ 1:n_lists ]
    attr( out, 'n_simulated' ) <- n_sim
    attr( out, 'fraction_inside_at_end' ) <- n_inside / n_sim
    return( out )
}

###
# Term co-occurrence matrix as in 01__emb.R (symmetric window, 1/d weighting, no counts across lists)
##
make_tcm <- function( lists, window = 5L )
{
    sep <- paste0( ' ', paste0( rep( '@', window * 2 ), collapse = ' ' ), ' ' )
    tokens <- text2vec::space_tokenizer( paste( lists, collapse = sep ) )
    it <- text2vec::itoken( tokens, progressbar = FALSE )
    vocab <- text2vec::create_vocabulary( it )
    vocab <- vocab[ vocab$term != '@', ]
    vectorizer <- text2vec::vocab_vectorizer( vocab )
    tcm <- text2vec::create_tcm( it, vectorizer, skip_grams_window = as.integer( window ) )
    return( list( tcm = tcm, vocab = vocab ) )
}

###
# Full symmetric co-occurrence matrix (dense) from the upper-triangular text2vec tcm
##
symmetric_tcm <- function( tcm )
{
    x <- as.matrix( tcm )
    x <- x + t( x ) - diag( diag( x ) )
    return( x )
}

###
# GloVe embedding as in 01__emb.R (main + context vectors), single thread for reproducibility
##
###
# rsparse's GloVe stops in the first epoch when the cost per non-zero cell exceeds 1. That heuristic aborts very
# low- and high-dimensional fits (rank 2, rank 100) on this corpus irrespective of the learning rate, although the
# fits converge normally without it. We relax the threshold for the dimension sweep only; the NaN check stays.
# (The main embedding in 01__emb.R never triggers the heuristic and is fitted with the unmodified package.)
##
relax_glove_cost_check <- function()
{
    for( gen in list( text2vec::GlobalVectors, rsparse::GloVe ) )
    {
        f <- gen$public_methods$fit_transform
        src <- deparse( f )
        if( !any( grepl( 'cost/n_nnz > 1)', src, fixed = TRUE ) ) ) next
        src <- sub( 'cost/n_nnz > 1)', 'cost/n_nnz > 1e+06)', src, fixed = TRUE )
        g <- eval( parse( text = src ) )
        environment( g ) <- environment( f )
        gen$set( 'public', 'fit_transform', g, overwrite = TRUE )
    }
    invisible( TRUE )
}

fit_glove <- function( tcm, rank = 50, x_max = 100, n_iter = 100, seed = 444, learning_rates = c( 0.15, 0.05, 0.01 ), relax = FALSE )
{
    if( relax ) relax_glove_cost_check()
    # text2vec stops when the cost explodes; in that case retry with a smaller learning rate (recorded as attribute)
    for( lr in learning_rates )
    {
        set.seed( seed )
        glove <- text2vec::GlobalVectors$new( rank = rank, x_max = x_max, learning_rate = lr, alpha = 0.75 )
        wv_main <- tryCatch( suppressMessages( glove$fit_transform( tcm, n_iter = n_iter, convergence_tol = 0.00001, n_threads = 1 ) ),
                             error = function( e ) NULL )
        if( !is.null( wv_main ) ) break
    }
    if( is.null( wv_main ) ) stop( 'GloVe did not converge' )
    emb <- wv_main + t( glove$components )
    attr( emb, 'learning_rate' ) <- lr
    return( emb )
}

###
# Positive pointwise mutual information of a symmetric count matrix
##
ppmi <- function( x )
{
    total <- sum( x )
    pr <- rowSums( x ) / total
    pmi <- log( ( x / total ) / outer( pr, pr ) )
    pmi[ !is.finite( pmi ) | pmi < 0 ] <- 0
    return( pmi )
}

###
# Truncated SVD of PPMI (rank k), U * sqrt(S) as in Levy, Goldberg & Dagan (2015)
##
svd_ppmi <- function( x, k = 50 )
{
    s <- svd( ppmi( x ), nu = k, nv = 0 )
    emb <- s$u[ , 1:k ] %*% diag( sqrt( s$d[ 1:k ] ) )
    rownames( emb ) <- rownames( x )
    return( emb )
}

###
# word2vec skip-gram (dimension 50, window 5) on the lists
# Note: the word2vec package has no seed argument and is not exactly reproducible, even with threads = 1.
# On the main corpus, three repeated fits gave rho 0.796-0.802 and neighbour precision 0.397-0.407,
# i.e. within the spread over the 10 simulation seeds. 'seed' only sets R's generator.
##
fit_word2vec <- function( lists, dim = 50, window = 5, seed = 444, iter = 20 )
{
    set.seed( seed )
    model <- word2vec::word2vec( x = lists, type = 'skip-gram', dim = dim, window = window,
                                 iter = iter, min_count = 1, threads = 1, split = c( ' \n', '\n' ),
                                 stopwords = character( 0 ), sample = 0, hs = FALSE, negative = 5 )
    emb <- as.matrix( model )
    emb <- emb[ rownames( emb ) != '</s>', , drop = FALSE ]
    return( emb )
}
