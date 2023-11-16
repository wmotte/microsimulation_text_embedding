# https://text2vec.org/glove.html

# https://medium.com/cmotions/nlp-with-r-part-2-training-word-embedding-models-and-visualize-results-ae444043e234

library( "text2vec" )
library( "umap" )
library( "ggplot2" )

# input data
df <- readr::read_lines( 'out.00.sim/plain_text.txt.gz' )

# 491,185 words
sum( stringr::str_count( df, ' ' ) + 1 )

# input to GloVe is a single line, but we do not want word counts to influence at boundaries
# therefore, we concat with a dummy term in between and add that to stop words
list_separator <- paste0( " ", paste0( rep( "@", 10 ), collapse = ' ' ), " " )

# get single string
single_df <- paste( df, collapse = list_separator )

# Create iterator over tokens
tokens <- space_tokenizer( single_df )

# Create vocabulary. Terms will be unigrams (simple words).
it <- itoken( tokens, progressbar = TRUE )
vocab <- create_vocabulary( it, stopwords = c( "@" ) )

# prune minimal
vocab <- prune_vocabulary( vocab, term_count_min = 100 )

# Use our filtered vocabulary
vectorizer <- vocab_vectorizer( vocab )

# use window of n context words
# term-co-occurrence matrix (TCM).
tcm <- create_tcm( it, vectorizer, skip_grams_window = 10 )

set.seed( 444 )

# 50 vector length (x_max = 10, iter = 100, tol = 0.0001, threads = 4)
glove <- GlobalVectors$new( rank = 50, x_max = 100 )
wv_main <- glove$fit_transform( tcm, n_iter = 50, convergence_tol = 0.0001, n_threads = 4 )

# get context matrix
wv_context <- glove$components

# 50 x 120
dim( wv_context )

# combine main embedding and context embedding (sum) into one matrix
glove_embedding <- wv_main + t( wv_context )

# container
all <- NULL

vnames <- rownames( as.data.frame( glove_embedding ) )
vname <- vnames[ 1 ]

for( vname in vnames )
{
    # get single word vector
    single <- glove_embedding[ vname, , drop = FALSE ]
    
    # pairwise similarities
    cos_sim <- text2vec::sim2( x = glove_embedding, y = single, method = "cosine", norm = "l2" )
    
    # get highest correspond (n=5)
    corr <- sort( cos_sim[ , 1 ], decreasing = TRUE )[ 2:6 ]
    
    # into d.f.
    data <- data.frame( product = vname, closest_by = names( corr ), similarity = round( corr, 3 ) )

    # merge into container
    all <- rbind( all, data )
    
}
   
rownames( all ) <- NULL 
all


#######################

# glove dimension reduction
glove_umap <- umap( glove_embedding, 
                    n_components = 2, metric = "cosine", 
                    n_neighbors = 5, min_dist = 0.1, spread = 20 )

# dimensions of end result [120 x 2]
dim( glove_umap$layout )

# do the same for the GloVe embeddings
df_glove_umap <- as.data.frame( glove_umap$layout, stringsAsFactors = FALSE )

# Add the labels of the words to the dataframe
df_glove_umap$word <- rownames( glove_embedding )
colnames( df_glove_umap ) <- c( "UMAP1", "UMAP2", "word" )
df_glove_umap$technique <- 'GloVe'
cat( paste0('\n', 'Our GloVe embedding reduced to 2 dimensions:', '\n') )
str( df_glove_umap )

df_umap <- df_glove_umap


# Plot the UMAP dimensions for both Word2Vec and GloVe
ggplot( df_umap ) +
    geom_label( aes( x = UMAP1, y = UMAP2, label = word ) ) +
        geom_point( aes( x = UMAP1, y = UMAP2 ), colour = 'blue', size = 5 ) 
+

    #facet_wrap(~technique) +
    #labs(title = "Word embedding in 2D using UMAP") +
    #theme(plot.title = element_text(hjust = .5, size = 14))
