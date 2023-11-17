#!/usr/bin/env Rscript
#
# Process glove embedding with superheat clustering
#
################################################################################
library( "superheat" )
library( "dplyr" ) # for %>%
library( "wordspace" ) # for norm-2 normalization
################################################################################
# FUNCTIONS
################################################################################

###
# Get section (i.e., groups)
##
get_sections <- function()
{
    # get group d.f.
    sections <- rbind(
        data.frame( label = c('Apples', 'Bananas', 'Grapefruit', 'Grapes', 'Kiwi', 'Lime', 'Mangoes', 
                              'Oranges', 'Pineapples', 'Strawberries', 'Watermelons' ), section = '4_fruits' ),
        
        data.frame( label = c( 'Cauliflower', 'Cucumbers', 'Eggplant', 'Onion', 'Spinach', 
                               'Tomatoes', 'Peppers', 'Zucchini' ), section = '3_vegetables' ),
        
        data.frame( label = c( 'BitterLemon', 'Cassis', 'Coffee', 'Cola', 'Lemonade', 
                               'Sodawater', 'Tea' ), section = '2_beverages' ),
        
        data.frame( label = c( 'CaramelBars', 'ChewingGum', 'ChocolateBar', 'Crackers', 'GummyBears', 
                               'MixedNuts', 'Popcorn', 'PotatoChips', 'Snickers', 'Sweets', 'Snacks', 
                               'Twix', 'Pringles', 'ShoppingBag' ), section = '1_extras' ) )
    
    # reverse and do stuff with factor ordering (for proper plotting)
    sections <- sections[ nrow( sections ):1, ]
    sections$section <- as.factor( sections$section )
    levels( sections$section ) <- stringr::str_split_fixed( levels( sections$section ), "_", 2 )[ , 2 ]
    
    return( sections )
}

###
# Return common order of cols
##
get_order_cols <- function( input_matrix )
{
    # get max value to normalize matrix
    max_value <- max( abs( input_matrix ) )
    
    # normalize and transpose
    mat <- t( input_matrix / max_value )
    
    # get sections and invert
    ss <- get_sections()
    
    # order rows as sections
    mat <- t( mat[ , ss$label ] )
    
    # cluster
    hclust_cols <- hclust( dist( t( mat ), method = "euclidean" ) )
    output <- hclust_cols$order
    
    return( output )
}

###
# Get cosine between two vectors
##
cosine_fun <- function( x, y )
{
    # calculate the cosine similarity between two vectors: x and y
    c <- sum(x*y) / (sqrt(sum(x * x)) * sqrt(sum(y * y)))
    return(c)
}

###
# Get cosine matrix
##
cosine_sim <- function( X )
{
    # calculate the pairwise cosine similarity between columns of the matrix X.
    # initialize similarity matrix
    m <- matrix(NA, 
                nrow = ncol(X),
                ncol = ncol(X),
                dimnames = list(colnames(X), colnames(X)))
    cos <- as.data.frame(m)
    
    # calculate the pairwise cosine similarity
    for(i in 1:ncol(X)) {
        for(j in i:ncol(X)) {
            co_rate_1 <- X[which(X[, i] & X[, j]), i]
            co_rate_2 <- X[which(X[, i] & X[, j]), j]  
            cos[i, j] <- cosine_fun(co_rate_1, co_rate_2)
            # fill in the opposite diagonal entry
            cos[j, i] <- cos[i, j]        
        }
    }
    return(cos)
}

###
# calculate the cosine silhouette width, which in cosine land is 
# (1) the lowest average dissimilarity of the data point to any other cluster, 
#  minus
# (2) the average dissimilarity of the data point to all other data points in 
#     the same cluster
# https://rlbarter.github.io/superheat-examples/Word2Vec/
##
cosine_silhouette <- function( cosine.matrix, membership ) 
{
    # Args:
    #   cosine.matrix: the cosine similarity matrix for the words
    #   membership: the named membership vector for the rows and columns. 
    #               The entries should be cluster centers and the vector 
    #               names should be the words.
    if (!is.factor(membership)) {
        stop("membership must be a factor")
    }
    # note that there are some floating point issues:
    # (some "1" entires are actually slightly larger than 1)
    cosine.dissim <- acos(round(cosine.matrix, 10)) / pi
    widths.list <- lapply(levels(membership), function(clust) {
        # filter rows of the similarity matrix to words in the current cluster
        # filter cols of the similarity matrix to words in the current cluster
        cosine.matrix.inside <- cosine.dissim[membership == clust, 
                                              membership == clust]
        # a: average dissimilarity of i with all other data in the same cluster
        a <- apply(cosine.matrix.inside, 1, mean)
        # filter rows of the similarity matrix to words in the current cluster
        # filter cols of the similarity matrix to words NOT in the current cluster
        other.clusters <- levels(membership)[levels(membership) != clust]
        cosine.matrix.outside <- sapply(other.clusters, function(other.clust) {
            cosine.dissim[membership == clust, membership == other.clust] %>%
                apply(1, mean) # average over clusters
        })
        # b is the lowest average dissimilarity of i to any other cluster of 
        # which i is not a member
        b <- apply(cosine.matrix.outside, 1, min)
        # silhouette width is b - a
        cosine.sil.width <- b - a
        data.frame(word = names(cosine.sil.width), width = cosine.sil.width)
    })
    widths.list <- do.call(rbind, widths.list)
    # join membership onto data.frame
    membership.df <- data.frame(word = names(membership), 
                                membership = membership)
    widths.list <- left_join(widths.list, membership.df, by = "word")
    
    return( widths.list )
}

###
# Make heatmap
##
make_heatmap <- function( input_matrix, sname, outdir, order_cols )
{
    # norm-2 normalization (i.e., dot product of identity vectors should be 1)
    mat <- t( wordspace::normalize.rows( input_matrix, method = "euclidean" ) )
    
    # get sections and invert
    ss <- get_sections()

    # order rows as sections
    mat <- t( mat[ , ss$label ] )

    # order cols after clustering of main matrix
    mat <- mat[ , order_cols ]
    
    # get similarity matrix
    simil <- cosine_sim( t( mat ) )
    
    ###### PLOT 1 ####
    
    # save to disk
    outfile <- paste0( outdir, '/', sname, '__similarity_matrix.png' )
    png( outfile, height = 3200, width = 3200, res = 300 )
    
    # plot similarity matrix
    superheat( simil, 
              # place dendrograms on columns and rows 
              row.dendrogram = F, col.dendrogram = F,
              
              # make gridlines white for enhanced prettiness
              grid.hline.col = "gray30",
              grid.vline.col = "gray30",
              
              # rotate bottom label text
              bottom.label.text.angle = -90,
              
              left.label.text.size = 4,
              bottom.label.text.size = 4,
              legend = FALSE )
    
    dev.off()
    
    # cosine-silhouette width
    #membership <- ss$section
    #names( membership ) <- ss$label
    #cs_width <- cosine_silhouette( mat, membership = as.factor( membership ) )$width   

    ######## PLOT 2 #######
        
    # save to disk
    outfile <- paste0( outdir, '/', sname, '__sections.png' )
    png( outfile, height = 2700, width = 4000, res = 300 )

    # sections [fruit, vegetables, beverages, extras]
    superheat( mat, 
               
               # vertical right bar with cosine silhouette width
               # https://rlbarter.github.io/superheat-examples/Word2Vec/
               #yr = cs_width,
               #yr.plot.type = "bar",
               #yr.bar.col = "grey30",
               #yr.axis.name = 'width',

               grid.hline.col = "gray90",
               grid.hline.size = 0.8,
               
               grid.vline.col = "gray30",
               bottom.label = "none",
               
               left.label.text.size = 5,
               left.label.text.alignment = "center",
               legend = FALSE,
               pretty.order.cols = FALSE,
               scale = TRUE,
               membership.rows = ss$section )
    
    dev.off()
    
    
    ##### PLOT 3 ######
    
    # save to disk
    outfile <- paste0( outdir, '/', sname, '__items.png' )
    png( outfile, height = 2700, width = 4000, res = 300 )
    
    # individual items
    superheat( mat, 
               grid.hline.col = "gray30",
               grid.vline.col = "gray30",
               
               bottom.label = "none",
              
               left.label.text.size = 4,
               left.label.text.alignment = "center",
               legend = FALSE,
               pretty.order.cols = FALSE, scale = TRUE )
    
    dev.off() 
    
}

################################################################################
# END FUNCTIONS
################################################################################

# load vocab, tcm, glove, wv_main, wv_context
load( "out.01.emb/saved_glove.RData" )

# output dir
outdir <- 'out.03.superheat'
dir.create( outdir, showWarnings = FALSE )

# plot matrix
#image( wv_main )    # 120 x 50
#image( wv_context ) # 50 x 120
#image( embedding )  # 120 x 50


# loop over types
type <- 'I'

for( type in c( 'I', 'II', 'III' ) )
{
    identifier <- paste0( "_", type, "$" )
    
    # only select products 'I' (i.e., 40 x 50 matrix) -> wv_main
    wv_main_small <- wv_main[ grep( identifier, rownames( wv_main ) ), ]
    rownames( wv_main_small ) <- gsub( identifier, "", rownames( wv_main_small ) )
    
    # only select products 'I' (i.e., 40 x 50 matrix) -> wv_context
    twv_context <- t( wv_context )
    wv_context_small <- twv_context[ grep( identifier, rownames( twv_context ) ), ]
    rownames( wv_context_small ) <- gsub( identifier, "", rownames( wv_context_small ) )
    
    # only select products 'I' (i.e., 40 x 50 matrix) -> average of main and context
    embedding_small <- embedding[ grep( identifier, rownames( embedding ) ), ]
    rownames( embedding_small ) <- gsub( identifier, "", rownames( embedding_small ) )

    # get hierarchical clustering order of main embedding (to use for other matrices as well)
    order_cols <- get_order_cols( embedding_small )
    
    # write heat maps of all dimensions as matrices
    make_heatmap( wv_main_small, paste0( 'wv_main_type_', type ), outdir, order_cols )
    make_heatmap( wv_context_small, paste0( 'wv_context_type_', type ), outdir, order_cols )
    make_heatmap( embedding_small, paste0( 'wv_embedding_type_', type ), outdir, order_cols )
    
}

