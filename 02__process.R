#!/usr/bin/env Rscript
#
# Process glove embedding
#
################################################################################
library( "umap" )
library( "ggplot2" )

################################################################################
# FUNCTIONS
################################################################################

###
# Get top-10 neighbors
##
calculate_neighbors <- function( embedding )
{
    # container
    all <- NULL
    
    vnames <- rownames( as.data.frame( embedding ) )
    vname <- vnames[ 1 ]
    
    for( vname in vnames )
    {
        # get single word vector
        single <- embedding[ vname, , drop = FALSE ]
        
        # pairwise similarities
        cos_sim <- text2vec::sim2( x = embedding, y = single, method = "cosine", norm = "l2" )
        
        # get highest correspond (n=10)
        corr <- sort( cos_sim[ , 1 ], decreasing = TRUE )[ 2:11 ]
        
        # into d.f.
        data <- data.frame( product = vname, closest_by = names( corr ), similarity = round( corr, 3 ) )
        
        # merge into container
        all <- rbind( all, data )
        
    }
    
    rownames( all ) <- NULL 
    return( all )
}

################################################################################
# END FUNCTIONS
################################################################################

# load vocab, tcm, glove, wv_main, wv_context
load( "out.01.emb/saved_glove.RData" )

# output dir
outdir <- 'out.02.process'
dir.create( outdir, showWarnings = FALSE )


# TODO plot matrix
image( wv_main )
image( wv_context )
image( embedding )

# neighbors
nn <- calculate_neighbors( embedding )

# write to file
readr::write_csv( nn, file = paste0( outdir, '/nearest_neighbors.csv' ), quote = 'all' )



# only select products 'I' (i.e., 40 x 50 matrix)
embedding_small <- embedding[ grep( "_I$", rownames( embedding ) ), ]

# clean names from suffix
rownames( embedding_small ) <- gsub( "_I$", "", rownames( embedding_small ) )

#######################

# glove dimension reduction
glove_umap <- umap( embedding_small, n_components = 2, spread = 1 )
                    #n_components = 2, metric = "cosine", 
                    #n_neighbors = 5, min_dist = 0.1, spread = 15 )

# dimensions of end result [40 x 2]
dim( glove_umap$layout )

# do the same for the GloVe embeddings
df_glove_umap <- as.data.frame( glove_umap$layout, stringsAsFactors = FALSE )

# Add the labels of the words to the dataframe
df_glove_umap$word <- rownames( embedding_small )
colnames( df_glove_umap ) <- c( "UMAP1", "UMAP2", "word" )



# Plot the UMAP dimensions for both Word2Vec and GloVe
#ggplot( df_glove_umap ) +
#    geom_label( aes( x = UMAP1, y = UMAP2, label = word ) ) +
#        geom_point( aes( x = UMAP1, y = UMAP2 ), colour = 'blue', size = 2, alpha = 0.6 ) 


set.seed(42)


p + geom_label_repel(aes(label = rownames(df),
                         fill = factor(cyl)), color = 'white',
                     size = 3.5) +
    theme(legend.position = "bottom")


    #facet_wrap(~technique) +
    #labs(title = "Word embedding in 2D using UMAP") +
    #theme(plot.title = element_text(hjust = .5, size = 14))

library( "ggrepel" )


df_glove_umap$group <- 'Vegetable'

p <- 
    ggplot( df_glove_umap, aes( x = UMAP1, y = UMAP2 ) ) +
       geom_point( color = 'red' ) +
       geom_label_repel( aes( label = word, fill = group ), color = 'white', segment.colour="gray30", size = 3.5 ) +
       theme( legend.position = "top" )

head( df )

head( df_glove_umap )


