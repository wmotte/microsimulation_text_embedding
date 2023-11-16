# https://text2vec.org/glove.html

# https://medium.com/cmotions/nlp-with-r-part-2-training-word-embedding-models-and-visualize-results-ae444043e234

library( "text2vec" )
library( "umap" )


# single line
# wiki <- readLines( 'misc/text8', n = 1, warn = FALSE )


df <- readr::read_lines( 'out.00.sim/plain_text.txt' )

# 17,005,207 words
stringr::str_count( df, ' ' ) + 1


# Create iterator over tokens
tokens <- space_tokenizer( df )

# Create vocabulary. Terms will be unigrams (simple words).
it <- itoken( tokens, progressbar = TRUE )
vocab <- create_vocabulary( it )

vocab <- prune_vocabulary( vocab, term_count_min = 100 )

# Use our filtered vocabulary
vectorizer <- vocab_vectorizer( vocab )

# use window of 5 for context words
# term-co-occurence matrix (TCM).
tcm <- create_tcm( it, vectorizer, skip_grams_window = 3 )

set.seed( 1 )

# glove
#glove <- GlobalVectors$new( rank = 16, x_max = 1 )


library( 'rsparse' )

glove_model <- NULL
glove_model = GloVe$new( rank = 100, x_max = 1, learning_rate = 0.00000015, lambda = 0.01, shuffle = TRUE )
embeddings = glove_model$fit_transform( tcm, n_iter = 10, n_threads = 4 )

#this->vocab_size = as<size_t>(params["vocab_size"]);
#this->word_vec_size = as<size_t>(params["word_vec_size"]);
#this->x_max = as<uint32_t>(params["x_max"]);
#this->learning_rate = as<T>(params["learning_rate"]);
#this->alpha = as<T>(params["alpha"]);
#this->lambda = as<T>(params["lambda"]);

for( i in 10:100 )
{
    print( i )
    # maximum number of co-occurrences to use in the weighting function, we choose the entire token set divided by 50 => 0.8 => 1
    glove <- NULL
    glove <- GlobalVectors$new( rank = i, x_max = 1 )
    glove$initialize( rank = i, x_max = 1, learning_rate = 0.0000001, lambda = 0.00001, alpha = 1 )
    wv_main <- glove$fit_transform( tcm, n_iter = 10, convergence_tol = -1, n_threads = 1 )
    Sys.sleep( 5 )
}

wv_context <- glove$components
dim( wv_context )

# combine main embedding and context embeddings (sum) into one matrix
word_vectors <- wv_main + t( wv_context )

# container
all <- NULL

vnames <- rownames( as.data.frame( word_vectors ) )

for( vname in vnames )
{

    single <- word_vectors[ vname, , drop = FALSE ]
    cos_sim = sim2( x = word_vectors, y = single, method = "cosine", norm = "l2" )
    
    # get highest correspond (n=2)
    corr <- sort( cos_sim[ , 1 ], decreasing = TRUE )[2:4]
        
    data <- data.frame( vname = vname, closest_corresponence = names( corr ) )

    all <- rbind( all, data )
    
}
   
rownames( all ) <- NULL 
all


#######################



glove_embedding <- word_vectors

# GloVe dimension reduction
glove_umap <- umap( word_vectors, n_components = 2, metric = "cosine", n_neighbors = 15, min_dist = 0.1, spread = 2 )

# Dimensions of end result
dim(glove_umap$layout)


# Do the same for the GloVe embeddings
df_glove_umap <- as.data.frame( glove_umap$layout, stringsAsFactors = FALSE)

# Add the labels of the words to the dataframe
df_glove_umap$word <- rownames( glove_embedding )
colnames(df_glove_umap) <- c("UMAP1", "UMAP2", "word")
df_glove_umap$technique <- 'GloVe'
cat(paste0('\n', 'Our GloVe embedding reduced to 2 dimensions:', '\n'))
str(df_glove_umap)
ls()
df_umap <- df_glove_umap

library( 'ggplot2' )
# Plot the UMAP dimensions for both Word2Vec and GloVe
ggplot(df_umap) +
    geom_point( aes( x = UMAP1, y = UMAP2 ), colour = 'blue', size = 5) +
    geom_label( aes( x = UMAP1, y = UMAP2, label = word ) )
    #facet_wrap(~technique) +
    #labs(title = "Word embedding in 2D using UMAP") +
    #theme(plot.title = element_text(hjust = .5, size = 14))
