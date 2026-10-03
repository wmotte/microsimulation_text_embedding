#!/usr/bin/env Rscript
#
# https://text2vec.org/glove.html
#
# https://medium.com/cmotions/nlp-with-r-part-2-training-word-embedding-models-and-visualize-results-ae444043e234
#
################################################################################
library( "text2vec" )

# output dir
outdir <- 'out.01.emb'
dir.create( outdir, showWarnings = FALSE )

# input data
df <- readr::read_lines( 'out.00.sim/plain_text.txt.gz' )

# 804,422 words
sum( stringr::str_count( df, ' ' ) + 1 )

# input to GloVe is a single line, but we do not want word counts to influence at boundaries
# therefore, we concat with a dummy term in between
list_separator <- paste0( " ", paste0( rep( "@", 10 ), collapse = ' ' ), " " )

# get single string
single_df <- paste( df, collapse = list_separator )

# Create iterator over tokens
tokens <- space_tokenizer( single_df )

# Create vocabulary. Terms will be unigrams (simple words).
it <- itoken( tokens, progressbar = TRUE )
vocab <- create_vocabulary( it )

# prune to get rid of maximal "@"
vocab <- prune_vocabulary( vocab, term_count_max = vocab[ vocab$term == '@', 'term_count' ] - 1 )

# write to file
readr::write_tsv( vocab, file = paste0( outdir, '/vocab_summary.tsv' ), quote = 'all' )

# Use our filtered vocabulary
vectorizer <- vocab_vectorizer( vocab )

# use window of n context words
# term-co-occurrence matrix (TCM).
tcm <- create_tcm( it, vectorizer, skip_grams_window = 5L )

set.seed( 444 )

# 50 vector length (x_max = 100, alpha = 0.75, learning rate = 0.15, iter = 100, tol = 0.00001)
# single thread: multi-threaded (asynchronous) updates are not reproducible from a seed
glove <- GlobalVectors$new( rank = 50, x_max = 100, learning_rate = 0.15, alpha = 0.75 )
wv_main <- glove$fit_transform( tcm, n_iter = 100, convergence_tol = 0.00001, n_threads = 1 )

# get context matrix
wv_context <- glove$components

# 50 x 120
dim( wv_context )

# combine main embedding and context embedding (sum) into one matrix
embedding <- wv_main + t( wv_context )

# save files to disk
save( vocab, tcm, glove, wv_main, wv_context, embedding, file = paste0( outdir, "/saved_glove.RData" ) )

# write settings and software versions (for reproducibility)
settings <- data.frame( rank = 50, x_max = 100, alpha = 0.75, learning_rate = 0.15, n_iter = 100, 
                        convergence_tol = 0.00001, n_threads = 1, seed = 444, skip_grams_window = 5, 
                        weighting = '1/distance, symmetric, no counts across lists', 
                        final_cost = round( tail( glove$get_history()$cost_history, 1 ), 6 ), 
                        n_iter_run = length( glove$get_history()$cost_history ) )
readr::write_tsv( settings, file = paste0( outdir, '/glove_settings.tsv' ) )
writeLines( capture.output( sessioninfo::session_info() ), paste0( outdir, '/session_info.txt' ) )


 