# R_Exercise-Graphs.R
#
# Purpose:  Exercises for working with graphs in R.
#
# Version: 1.4
#
# Date:    2016  02  05
# Author:  Boris Steipe
#
# V 1.4    add HotNet subnetwork discovery example
# V 1.3    Update from 2016 "Function" code; add heat-diffusion example
# V 1.2    Typos
# V 1.1    Added a biomaRt section
# V 1.0    First code
#
# TODO:
#  - storing extra annotations
#  - compare networks from different sources
#  - create a fusion network
#  - graph output formats
#  - demonstrate igraph / cytoscape interface
#  - More network visualization
#    https://rpubs.com/kateto/netviz
#    http://kateto.net/networks-r-igraph
#  - Biofabric and Hive-plots
#  http://stackoverflow.com/questions/22453273/how-to-visualize-a-large-network-in-r
#  http://datablick.com/2015/05/13/comb-the-hairball-with-biofabric-in-tableau-by-chris-demartini/
#  http://datablick.com/2015/04/13/circular-and-hive-plot-network-graphing-in-tableau-by-chris-demartini/
#
#  - Add GO analysis?
#
# == HOW TO WORK WITH THIS FILE ======================================
#
#  Go through this script line by line to read and understand the
#  code. Execute code by typing <cmd><enter>. When nothing is
#  selected, that will execute the current line and move the cursor to
#  the next line. You can also select more than one line, e.g. to
#  execute a block of code, or less than one line, e.g. to execute
#  only the core of a nested expression.
#
#  Edit code, as required, experiment with options, or just play.
#  Especially play.
#
#  DO NOT simply source() this whole file!
#
#  If there are portions you don't understand, use R's help system,
#  Google for an answer, or ask me. Don't continue if you don't
#  understand what's going on. That's not how it works ...
#
#
# ==============================================================================

# This tutorial covers basic concepts of graph theory and analysis in R. You
# should have typed init() to configure some utilities in the background.


# ==============================================================================
#        PART ONE: REVIEW
# ==============================================================================

# I assume you'll have read the Pavlopoulos review of graph theory concepts.
# Let's explore some of the ideas by starting with a small random graph."


# To begin let's write a little function that will create random "gene" names;
# there's no particular purpose to this other than to make our graphs look a
# little more like what we would find in a publication ...
makeRandomGenenames <- function(N) {
  nam <- character()
  while (length(nam) < N) {
    a <- paste(c(sample(LETTERS, 1), sample(letters, 2)),
               sep="", collapse="") # three letters
    n <- sample(1:9, 1)             # one number
    nam[length(nam) + 1] <- paste(a, n, sep="") # store in vector
    nam <- unique(nam)   # delete if this was a duplicate
  }
  return(nam)
}

N <- 20

set.seed(112358)
( Nnames <- makeRandomGenenames(N) )

# One way to represent graphs in a computer is as an "adjacency matrix". In this
# matrix, each row and each column represents a node, and the cell at the
# intersection of a row and column contains a value/TRUE if there is an edge,
# 0/FALSE otherwise. It's easy to see that an undirected graph has a symmetric
# adjacency matrix (i, j) == (j, i); and we can put values other than {1, 0}
# into a cell if we want to represent a weighted edge.

# At first, lets create a random graph: let's say a pair of nodes has
# probability p <- 0.1 to have an edge, and our graph is symmetric and has no
# self-edges. We use our Nnames as node labels, but I've written the function so
# that we could also just ask for any number of un-named nodes, we'll use that later.

makeRandomGraph <- function(nam, p = 0.1) {
  # nam: either a character vector of unique names, or a single
  #        number that will be converted into a vector of integers.
  # p:   probability that a random pair of nodes will have an edge.
  #
  # Value: an adjacency matrix
  #
  if (is.numeric(nam) && length(nam) == 1) { # if nam is  a single number ...
    nam <- as.character(1:nam)
  }
  N <- length(nam)
  G <- matrix(numeric(N * N), ncol = N)  # The adjacency matrix
  rownames(G) <- nam
  colnames(G) <- nam
  for (iRow in 1:(N-1)) { # Note how we make sure iRow != iCol
    for (iCol in (iRow+1):N) {
      if (runif(1) < p) {  # runif() creates uniform random numbers
        # between 0 and 1
        G[iRow, iCol] <- 1   # row, col !
        G[iCol, iRow] <- 1   # col, row !
      }
    }
  }
  return(G)
}

set.seed(112358)
( G <- makeRandomGraph(Nnames, p = 0.09) )


# Listing the matrix is not very informative - we should plot this graph. We'll
# go into more details of the igraph package a bit later, for now we just use it
# to plot:

if (!require(igraph)) {
  install.packages("igraph")
  library(igraph)
}

iG <- graph_from_adjacency_matrix(G)
iGxy <- layout_with_graphopt(iG, charge=0.001)   # calculate layout coordinates


# The igraph package adds its own function to the collection of plot()
# functions; R makes the selection which plot function to use based on the class
# of the object that we request to plot. This plot function has parameters
#  layout - the x,y coordinates of the nodes;
#  vertex.color - which I define to color by node-degree
#  vertex size - which I define to increase with node-degree
#  vertex.label - which I set to use our Nnames vector

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 800 + (150 * degree(iG)),
     vertex.label = as.character(degree(iG)/2),
     #     vertex.label = Nnames,
     edge.arrow.size = 0)
par(oPar)  # reset plot window


# ==============================================================================
#        PART TWO: COMPUTING WITH GRAPHS
# ==============================================================================

# To get a better feel for computing with graphs, let's write a simple
# algorithm: counting the number of connected components in the graph. A
# subgraph is connected if there is a path from any node to any other. As the
# plot above showed you, our graph G has three connected components.

# The algorithm employs a simple breadth-first search: we start at a node and
# initialize our first connected component with it. Then we add all its neigbors
# to the connected component. Then their neighbors, until there are none left to
# add. If there are remaining nodes that have not yet been processed, they must
# be part of another component. We initialize another component and repeat until
# all nodes have been processed.

# To write this, we need three data structures:
components <- list()             # A list to store the components of G
processed <- logical(ncol(G))    # A vector of nodes that have already been processed
names(processed) <- rownames(G)
queue <- character()             # A vector of nodes currently scheduled for processing

# Then we define a helper function
neighboursOf <- function(node, adjM) {
  # Returns the names of all neighbors of a node in an adjacency
  # matrix. One could easily substitute this function for a different graph
  # representation. This function assumes that elements of the adjacency
  # matrix evaluate as logical FALSE if they do NOT represent an edge (i.e.
  # they are 0, or 0.0) and TRUE if they do represent an edge (any non-zero
  # number). Thus this will work for weighted and unweighted graphs.
  #                                adjM[node, ]   contents of node's row of
  #                                               the matrix ...
  #                     as.logical(adjM[node, ])  ... turned into a logical
  #                                               vector ...
  #      colnames(adjM)[as.logical(adjM[node, ])] ... used to subset the
  #                                               return-value of colnames()
  return(colnames(adjM)[as.logical(adjM[node, ])])
}


for (node in rownames(G)) {           # Process each node.
  if (! processed[node]) {            # If we haven't processed it previously:
    iC <- length(components) + 1      #   initialize a new component
    components[[iC]] <- node          #   ... with this node.
    queue <- neighboursOf(node, G)    #   Append its neighbors to the queue
    processed[node] <- TRUE           #   Flag it as processed.
    while (length(queue) > 0) {       #   Now work through the queue:
      currentNode <- queue[1]         #     get the node from position one,
      queue <- queue[-1]              #     ... and remove it from the queue.
      if (! processed[currentNode]) { #     If we haven't yet processed it:
        components[[iC]] <- c(components[[iC]], currentNode) # ... add it to the
                                                             # current
                                                             # component,
        queue <- c(queue, neighboursOf(currentNode, G)) # ... add its neigbors
                                                        # to the queue
        processed[currentNode] <- TRUE      #  ... and flag it as processed.
      }
    }                                       #   Repeat until the queue is empty.
  }
}

# A brief note on using queues in R:

# R is a functional language. In clean R code, functions do not have side
# effects. In my experience this really makes a lot of difference for writing
# bug-free code that is easy to read. But sometimes this means we have to take
# extra steps. For example, some languages have functions to push() and pop()
# elements on/off stacks - or shift()/unshift(), or enqueue(), dequeue()
# elements from a queue. But a dequeue operation does two things: it retrieves a
# value AND it modifies the queue. In R, this normally requires two statements -
# that's the way we wrote it above. That's fine because we really DO make two
# changes to our data. But if you feel you MUST have a function that does this
# in a single expression, this can be done. See:
# https://gist.github.com/leeper/d6d085ac86d1e006167e  and
# https://gist.github.com/mpettis/b7bfeff282e3b052684f also
# https://www.researchgate.net/post/What_is_the_queue_data_structure_in_R Of
# course you could also hold a value and an array in a list that you modify with
# a single function. So there are ways to wrap this into one expression. But
# knowing that you CAN do this, doesn't mean you SHOULD.

# Now let's look at the results:


# size of each component
( sizes <- unlist(lapply(components, length)) )
length(sizes)                       # number of components
mean(sizes)                         # average size
sd(sizes)                           # sd of size
components                          # the actual components

# So far so good, and if you want to practice coding, you could perhaps write
# your own implementation of Dijkstra's algorithm to find shortest paths
# (https://en.wikipedia.org/wiki/Dijkstra's_algorithm). But of course most graph
# algorithms of interest are available in an R package. There are several graph
# packages in R, but currently the package of choice seems to be iGraph."


# ==============================================================================
#        PART TWO: DEGREE DISTRIBUTIONS
# ==============================================================================

# The simplest descriptor of a graph are the number of nodes, edges, and the
# degree-distribution. In our example, the number of nodes was given: N; the
# number of edges can easily be calculated from the adjacency matrix. In our
# matrix, we have entered 1 for every edge. Thus we simply sum over the matrix:
sum(G)

# Is that correct? Is that what you see in the plot?

# Yes and no: we entered every edge twice: once for a node [i,j], and again for
# the node [j, i]. Whether that is correct depends on what exactly we
# want to do with the matrix. If these were directed edges, we would need to
# keep track of them separately. Since we didn't intend them to be directed,
# we could divide the number of edges by 2. Why didn't we simply use an
# upper-triangular matrix? Because then we need to keep track of the ordering of
# edges if we want to know whether a particular edge exists or not. For example
# we could sort the nodes alphabetically, and make sure we always query a pair
# in alphabetical order. Then a triangular matrix would be efficient.

# What about the degree distribution? We can get that simply by summing over the
# rows (or the columns):"

rowSums(G)  # check this against the plot!

# Let's  plot the degree distribution in a histogram:
rs <- rowSums(G)
brk <- seq(min(rs)-0.5, max(rs)+0.5, by=1)  # define breaks for the histogram
hist(rs, breaks=brk, col="#A5CCF5",
     xlim = c(-1,8), xaxt = "n",
     main = "Node degrees", xlab = "Degree", ylab = "Number")  # plot histogram
axis(side = 1, at = 0:7)

# Note: I don't _have_ to define breaks, the hist() function usually does so
# quite well, automatically. But for this purpose I want the columns of the
# histogram to represent exactly one node-degree difference.

# A degree distribution is actually quite an important descriptor of graphs,
# since it is very sensitive to the generating mechanism. For biological
# networks, that is one of the key questions we are interested in: how was the
# network formed?

# Let's simulate a few graphs that are a bit bigger to get a better sense of
# their degree distributions:
#

# === random graph


set.seed(31415927)
G200 <- makeRandomGraph(200, p = 0.015)
iG200 <- graph_from_adjacency_matrix(G200)
iGxy <- layout_with_graphopt(iG200, charge=0.0001) # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG200,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iG200)+1))[degree(iG200)+1],
     vertex.size = 200 + (30 * degree(iG200)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

# This graph has thirteen singletons and one large, connected component. Many
# biological graphs look approximately like this.

# Calculate degree distributions
dg <- degree(iG200)/2   # here, we use the iGraph function degree()
                        # not rowsums() from base R.
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#A5CCF5",
     xlim = c(-1,11), xaxt = "n",
     main = "Node degrees", xlab = "Degree", ylab = "Number")  # plot histogram
axis(side = 1, at = 0:10)


# Note the characteristic peak of this distribution: this is not "scale-free".
# Here is a log-log plot of frequency vs. degree-rank:

( freqRank <- table(dg) )
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#A5CCF5",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a random network")

# === scale-free graph (Barabasi-Albert)

# What does one of those intriguing "scale-free" distributions look like? The
# iGraph package has a function to make random graphs according to the
# Barabasi-Albert model of scale-free graphs. It is: sample_pa(), where pa
# stands for "preferential attachment", one type of process that will yield
# scale-free distributions.


set.seed(31415927)
GBA <- sample_pa(200, power = 0.8)

iGxy <- layout_with_graphopt(GBA, charge=0.0001) # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(GBA,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(GBA)+1))[degree(GBA)+1],
     vertex.size = 200 + (30 * degree(GBA)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

# This is a very obviously different graph! Some biological networks have
# features that look like that - but in my experience the hub nodes are usually
# not that distinct. But then again, that really depends on the parameter
# "power". Feel encouraged to change "power" and get a sense for what difference
# this makes. Also: note that the graph has only a single component.

# What's the degree distribution of this graph?
( dg <- degree(GBA) )
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#A5D5CC",
     xlim = c(0,30), xaxt = "n",
     main = "Node degrees 200 nodes PA graph",
     xlab = "Degree", ylab = "Number")
axis(side = 1, at = seq(0, 30, by=5))

# Most nodes have a degree of 1, but one node has a degree of 28.

(freqRank <- table(dg))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#A5F5CC",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a preferential-attachment network")

# Sort-of linear, but many of the higher ranked nodes have a frequency of only
# one. That behaviour smooths out in larger graphs:
#
X <- sample_pa(100000, power = 0.8)  # 100,000 nodes
freqRank <- table(degree(X))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     xlab = "log(Rank)", ylab = "log(frequency)",
     pch = 21, bg = "#A5F5CC",
     main = "100,000 nodes in a random, scale-free network")
rm(X)

# === Random geometric graph

# Finally, let's simulate a random geometric graph and look at the degree
# distribution. Remember: these graphs have a high probability to have edges
# between nodes that are "close" together - an entirely biological notion.

# We'll randomly place our nodes in a box. Then we'll define the
# probability for two nodes to have an edge to be a function of their distance.

# Here is a function that makes such graphs. iGraph has sample_grg(), which
# connects nodes that are closer than a cutoff, the function I give you below is
# a bit more interesting since it creates edges according to a probability that
# is determined by a generalized logistic function of the distance. This
# sigmoidal function gives a smooth cutoff and creates more "natural" graphs.
# Otherwise, the function is very similar to the random graph function, except
# that we output the "coordinates" of the nodes together with the adjacency
# matrix. Lists FTW.
#
makeRandomGeometricGraph <- function(nam, B = 25, Q = 0.001, t = 0.6) {
  # nam: either a character vector of unique names, or a single
  #        number that will be converted into a vector of integers.
  # B, Q, t:   probability that a random pair (i, j) of nodes gets an
  #              edge determined by a generalized logistic function
  #              p <- 1 - 1/((1 + (Q * (exp(-B * (x-t)))))^(1 / 0.9)))
  #
  # Value: a list with the following components:
  #        G$mat : an adjacency matrix
  #        G$nam : labels for the nodes
  #        G$x   : x-coordinates for the nodes
  #        G$y   : y-coordinates for the nodes
  #
  nu <- 1  # probably not useful to change
  G <- list()

  if (is.numeric(nam) && length(nam) == 1) {
    nam <- as.character(1:nam)
  }
  G$nam <- nam
  N <- length(G$nam)
  G$mat <- matrix(numeric(N * N), ncol = N)  # The adjacency matrix
  rownames(G$mat) <- G$nam
  colnames(G$mat) <- G$nam
  G$x <- runif(N)
  G$y <- runif(N)
  for (iRow in 1:(N-1)) { # Same principles as in makeRandomGraph()
    for (iCol in (iRow+1):N) {
      # geometric distance ...
      d <- sqrt((G$x[iRow] - G$x[iCol])^2 +
                  (G$y[iRow] - G$y[iCol])^2)  # Pythagoras
      # distance dependent probability
      p <- 1 - 1/((1 + (Q * (exp(-B * (d-t)))))^(1 / nu))
      if (runif(1) < p) {
        G$mat[iRow, iCol] <- 1
        G$mat[iCol, iRow] <- 1
      }
    }
  }
  return(G)
}

# Getting the parameters of a generalized logistic right takes a bit of
# experimenting. If you are interested, you can try a few variations. Or you can
# look up the function at
# https://en.wikipedia.org/wiki/Generalised_logistic_function

# This function computes generalized logistics ...
# genLog <- function(x, B = 25, Q = 0.001, t = 0.5) {
#     # generalized logistic (sigmoid)
#     nu <- 1
#     return(1 - 1/((1 + (Q * (exp(-B * (x-t)))))^(1 / nu)))
# }
#
# ... and this code plots p-values over the distances we could encouter between
# our nodes: from 0 to sqrt(2) i.e. the diagonal of the unit sqaure in which we
# will place our nodes.
# x <- seq(0, sqrt(2), length.out = 50)
# plot(x, genLog(x), type="l", col="#AA0000", ylim = c(0, 1),
#      xlab = "d", ylab = "p(edge)")

# 200 node random geomteric graph
set.seed(112358)
GRG <- makeRandomGeometricGraph(200, t=0.4)


iGRG <- graph_from_adjacency_matrix(GRG$mat)
iGRGxy <- cbind(GRG$x, GRG$y) # use our node coordinates for layout

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])) * 1.1,
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iGRG)+1))[degree(iGRG)+1],
     vertex.size = 0.1 + (0.1 * degree(iGRG)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

# degree distribution:
( dg <- degree(iGRG)/2 )
brk <- seq(min(dg)-0.5, max(dg)+0.5, by=1)
hist(dg, breaks=brk, col="#FCD6E2",
     xlim = c(0, 25), xaxt = "n",
     main = "Node degrees: 200 nodes RG graph",
     xlab = "Degree", ylab = "Number")
axis(side = 1, at = c(0, min(dg):max(dg)))

# You'll find that this is kind of in-between the random, and the scale-free
# graph. We do have hubs, but they are not as extreme as in the scale-free case;
# and we have have no singletons, in contrast to the random graph.

( freqRank <- table(dg) )
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#FCD6E2",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "200 nodes in a random geometric network")



# ==============================================================================
#        PART FOUR: A CLOSER LOOK AT THE igraph PACKAGE
# ==============================================================================


# == BASICS ==========================================================

# The basic object of the igraph package is a graph object. Let's explore the
# first graph some more, the one we built with our random gene names:
summary(iG)

# This output means: this is an IGRAPH graph, with D = directed edges and N =
# named nodes, that has 20 nodes and 40 edges. For details, see
?print.igraph

mode(iG)
class(iG)

# This means an igraph graph object is a special list object; it is opaque in
# the sense that a user is never expected to modify its components directly, but
# through a variety of helper functions which the package provides. There are
# many ways to construct graphs - from adjacency matrices, as we have just done,
# from edge lists, or by producing random graphs according to a variety of
# recipes, called _games_ in this package.

# Two basic functions retrieve nodes "Vertices", and "Edges":
V(iG)
E(iG)

# As with many R objects, loading the package provides special functions that
# can be accessed via the same name as the basic R functions, for example:

print(iG)
plot(iG)

# ... where plot() allows the usual flexibility of fine-tuning the plot. We
# first layout the node coordinates with the Fruchtermann-Reingold algorithm - a
# force-directed layout that applies an ettractive potential along edges (which
# pulls nodes together) and a repulsive potential to nodes (so they don't
# overlap). Note the use of the degree() function to color and scale nodes and
# labels by degree and the use of the V() function to retrieve the vertex names.
# See ?plot.igraph for details."

iGxy <- layout_with_fr(iG)   # calculate layout coordinates

# Plot with some customizing parameters
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     vertex.color=heat.colors(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 9 + (2 * degree(iG)),
     vertex.label.cex = 0.5 + (0.05 * degree(iG)),
     edge.arrow.size = 0,
     edge.width = 2,
     vertex.label = toupper(V(iG)$name))
par(oPar)


# == Components

# The igraph function components() tells us whether there are components of the
# graph in which there is no path to other components.
components(iG)

# In the _membership_ vector, nodes are annotated with the index of the
# component they are part of. Sui7 is the only node of component 2, Cyj1 is in
# the third component etc. This is perhaps more clear if we sort by component
# index
sort(components(iG)$membership)

# Retrieving e.g. the members of the first component from the list can be done
# by subsetting:

                                 components(iG)$membership == 1  # logical ..
       components(iG)$membership[components(iG)$membership == 1]
names(components(iG)$membership)[components(iG)$membership == 1]



# == RANDOM GRAPHS AND GRAPH METRICS ===========================================


# Let's explore some of the more interesting, topological graph measures. We
# start by building a somewhat bigger graph. We aren't quite sure whether
# biological graphs are small-world, or random-geometric, or
# preferential-attachment ... but igraph has ways to simulate the basic ones
# (and we could easily simulate our own). Look at the following help pages:

?sample_gnm                      # see also sample_gnp for the Erdös-Rényi models
?sample_smallworld               # for the Watts & Strogatz model
?sample_pa                       # for the Barabasi-Albert model

# But note that there are many more sample_ functions. Check out the docs!

# Let's look at betweenness measures for our first graph: here: the nodes again
# colored by degree. Degree centrality states: nodes of higher degree are
# considered to be more central. And that's also the way the force-directed
# layout drawas them, obviously.

set.seed(112358)
iGxy <- layout_with_fr(iG)   # calculate layout coordinates
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 20 + (10 * degree(iG)),
     vertex.label = Nnames,
     edge.arrow.size = 0)
par(oPar)

# == Diameter

diameter(iG)  # The diameter of a graph is its maximum length shortest path.

# let's plot this path: here are the nodes ...
get_diameter(iG)

# ... and we can get the x, y coordinates from iGxy by subsetting with the node
# names. The we draw the diameter-path with a transparent, thick pink line:
lines(iGxy[get_diameter(iG),], lwd=10, col="#ff63a788")

# == Centralization scores

?centralize
# replot our graph, and color by log_betweenness:

bC <- centr_betw(iG)  # calculate betweenness centrality
nodeBetw <- bC$res
nodeBetw <- round(log(nodeBetw +1)) + 1

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     rescale = FALSE,
     xlim = c(min(iGxy[,1]), max(iGxy[,1])) * 1.1,
     ylim = c(min(iGxy[,2]), max(iGxy[,2])) * 1.1,
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 20 + (10 * degree(iG)),
     vertex.label = Nnames,
     edge.arrow.size = 0)
par(oPar)

# Note that the betweenness - the number of shortest paths that pass through a
# node, is in general higher for high-degree nodes - but not always: Eqr2 has
# higher betweenness than Itv7: this measure really depends on the detailed
# local topology of the graph."

# Can you use centr_eigen() and centr_degree() to calculate the respective
# values? That's something I would expect you to be able to do.
#
# Lets plot betweenness centrality for our random geometric graph:

bCiGRG <- centr_betw(iGRG)  # calculate betweenness centrality

nodeBetw <- bCiGRG$res
nodeBetw <- round((log(nodeBetw +1))^2.5) + 1

# colours and size proportional to betweenness

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])),
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])),
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 0.1 + (0.03 * nodeBetw),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

diameter(iGRG)
lines(iGRGxy[get_diameter(iGRG),], lwd=10, col="#ff335533")



# == CLUSTERING ================================================================

# Clustering finds "communities" in graphs - and depending what the edges
# represent, these could be complexes, pathways, biological systems or similar.
# There are many graph-clustering algorithms. One approach with many attractive
# properties is the Map Equation, developed by Martin Rosvall. See:
# http://www.ncbi.nlm.nih.gov/pubmed/18216267 and htttp://www.mapequation.org


iGRGclusters <- cluster_infomap(iGRG)
modularity(iGRGclusters) # ... measures how separated the different membership
# types are from each other
membership(iGRGclusters) # which nodes are in what cluster?
table(membership(iGRGclusters))  # how large are the clusters?

# The largest cluster has 48 members, the second largest has 25, etc.

# Lets plot our graph again, coloring the nodes of the first five communities by
# their cluster membership:

# first, make a vector with as many grey colors as we have communities ...
commColors <- rep("#f1eef6", max(membership(iGRGclusters)))
# ... then overwrite the first five with "real colors" - something like rust,
# lilac, pink, and mauve or so.
commColors[1:5] <- c("#980043", "#dd1c77", "#df65b0", "#c994c7", "#d4b9da")


oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iGRG,
     layout = iGRGxy,
     rescale = FALSE,
     xlim = c(min(iGRGxy[,1]), max(iGRGxy[,1])),
     ylim = c(min(iGRGxy[,2]), max(iGRGxy[,2])),
     vertex.color=commColors[membership(iGRGclusters)],
     vertex.size = 0.1 + (0.1 * degree(iGRG)),
     vertex.label = "",
     edge.arrow.size = 0)

par(oPar)



# ==============================================================================
#        PART FIVE: EXPLORE FUNCTIONAL EDGES IN THE HUMAN PROTEOME
# ==============================================================================

# In order for you to explore a real, biological network, I give you a
# dataframe of functional relationships of human proteins that I have downloaded
# from the STRING database. The full table has 8.5 million records, here is a
# subset of records with combined confidence scores > 980

# The selected set of edges with a confidence of > 980 is a dataframe with about
# 50,000 edges and 6,500 unique proteins. You can load the saved dataframe here
# (and also read more about what the numbers mean at
# http://www.ncbi.nlm.nih.gov/pubmed/15608232 ).

load("data/STRINGedges.RData")

head(STRINGedges)

# The protein identifiers are Ensemble IDs, but they are prefixed with the
# taxonomy ID  for which this subset of all STRING edges was created (here:
# 9606, homo sapiens). For clarity, we'll remove this prefix. gsub() would
# normally be the function of choice, but in this case, since the format is
# exactly defined ...
summary(nchar(c(STRINGedges$protein1, STRINGedges$protein2)))
#  ... we can go with the faster substr() approach:

STRINGedges$protein1 <- substr(STRINGedges$protein1, 6, 20)
STRINGedges$protein2 <- substr(STRINGedges$protein2, 6, 20)

# We use an igraph function to make a graph from this dataframe
?graph_from_data_frame

gSTR <- graph_from_data_frame(STRINGedges)

# CAUTION you DON'T want to plot a graph with 6,500 nodes and 50,000 edges -
# layout of such large graphs is possible, but requires specialized code. Google
# for <layout large graphs> if you are curious. Also, consider what one can
# really learn from plotting such a graph ...

# Of course simple computations on this graph are reasonably fast:

compSTR <- components(gSTR)
summary(compSTR) # our graph is fully connected!

dg <- degree(gSTR)
hist(log(dg), col="#FEE0AF")
# this actually does look rather scale-free

(freqRank <- table(dg))
plot(log10(as.numeric(names(freqRank)) + 1),
     log10(as.numeric(freqRank)), type = "b",
     pch = 21, bg = "#FEE0AF",
     xlab = "log(Rank)", ylab = "log(frequency)",
     main = "6,500 human proteins with high-confidence  from the  STRING network",
     cex.main = 0.9)

# This looks very scale-free indeed.
#
# Now explore some more:

# === CLIQUES   ========
# Let's find the largest cliques. Remember: a clique is a fully connected
# subgraph, i.e. a subgraph in which every node is connected to every other.
# Biological complexes often appear as cliques in interaction graphs.

clique_num(gSTR)
# The largest clique has 63 members.

largest_cliques(gSTR)[[1]]

# Pick one of the proteins and find out what this fully connected cluster of 63
# proteins is (you can simply Google for the ID). Is this expected?



# === BETWEENNESS CENTRALITY   =======================================

# Let's find the nodes with the 10 - highest betweenness centralities.
#
BC <- centr_betw(gSTR)

# remember: BC$res contains the results
head(BC$res)

BC$res[1]   # betweeness centrality of node 1 in the graph ...
# ... which one is node 1?
V(gSTR)[1]

# to get the ten-highest nodes, we simply label the elements BC with their
# index ...
names(BC$res) <- as.character(1:length(BC$res))

# ... and then we sort:
sBC <- sort(BC$res, decreasing = TRUE)
head(sBC)

# This ordered vector means: node 3,862 has the highest betweeness centrality,
# node 1,720 has the second highest.. etc.

BCsel <- as.numeric(names(sBC)[1:10])
BCsel
# We can use the first ten labels to subset the nodes in gSTR and fetch the
# IDs...
ENSPsel <- names(V(gSTR)[BCsel])


# Could you define in a short answer quiz what these IDs are? And what their
# biological significance is? I expect you to be able to.

#  Next, to find what these proteins are...

# We could now Google for all of these IDs to learn more about them. But really,
# googling for ID one after the other, that would be lame. Let's instead use
# the very, very useful biomaRt package to translate these Ensemble IDs into
# gene symbols.

# == biomaRt =========================================================

# IDs are just labels, but for _bio_informatics we need to learn more about the
# biological function of the genes or proteins that we retrieve via graph data
# mining. biomaRt is the tool of choice. It's a package distributed by the
# bioconductor project. This here is not a biomaRt tutorial (that's for another
# day), simply a few lines of sample code to get you started on the specific use
# case of retrieving descriptions for ensembl protein IDs.

if (!require(biomaRt)) {
  source("http://bioconductor.org/biocLite.R")
  biocLite("biomaRt")
  library("biomaRt")
}

# define which dataset to use ...
myMart <- useMart("ensembl", dataset="hsapiens_gene_ensembl")

# what filters are defined?
filters <- listFilters(myMart)
filters

# and what attributes can we filter for?
attributes <- listAttributes(myMart)
attributes

# Soooo many options - let's look for the correct name of filters that are
# useful for ENSP IDs ...
filters[grep("ENSP", filters$description), ]

# ... and the correct attribute names for gene symbols and descriptions ...
attributes[grep("symbol", attributes$description, ignore.case=TRUE), ]
attributes[grep("description", attributes$description, ignore.case=TRUE), ]


# ... so we can put this together: here is a syntax example:
getBM(filters = "ensembl_peptide_id",
      attributes = c("hgnc_symbol",
                     "wikigene_description",
                     "interpro_description",
                     "phenotype_description"),
      values = "ENSP00000000442",
      mart = myMart)

# A simple loop will now get us the information for our 10 most central genes
# from the human subset of STRING.

CPdefs <- list()  # Since we don't know how many matches one of our queries
# will return, we'll put the result dataframes into a list.

for (ID in ENSPsel) {
  CPdefs[[ID]] <- getBM(filters = "ensembl_peptide_id",
                        attributes = c("hgnc_symbol",
                                       "wikigene_description",
                                       "interpro_description",
                                       "phenotype_description"),
                        values = ID,
                        mart = myMart)
}

# So what are the proteins with the ten highest betweenness centralities?
#  ... are you surprised? (I am! Really.)

# Final task: Write a loop that will go through your list and
#    for each ID:
#    --  print the ID,
#    --  print the first row's symbol, and
#    --  print the first row's wikigene description.
#
# (Hint, you can structure your loop in the same way as the loop that
# created CPdefs. )

# Print the R code for your loop and its output for the ten genes onto a sheet
# of paper.



# ==============================================================================
#        PART SIX: ATTRIBUTE PROPAGATION
# ==============================================================================

# ==== preparations - if you start at this section
library(igraph)
load("data/STRINGedges.RData")
STRINGedges$protein1 <- substr(STRINGedges$protein1, 6, 20)
STRINGedges$protein2 <- substr(STRINGedges$protein2, 6, 20)
gSTR <- graph_from_data_frame(STRINGedges)
# ================================================


# In graph theory, a label is a categorical attribute of a vertex or edge,
# typically represented by a boolean, or integer. Of course, attributes in
# general don't have to be categorical, but can be anything, however when we
# speak of "attribute propagation", we are considering real-valued attributes,
# and how network topology can help us to infer attributes by using knowledge
# about attributes of neighboring or close vertices. Label propagation on the
# other hand considers whether we can assign categories to nodes for which
# category membership is unknown, based on their proximity to nodes of known
# categories. The two tasks intersect when we propagate attributes in a network
# to determine vertex labels based on a threshold value of the attribute.

# In this section we will build a subgraph from the STRING network of human,
# functionally interacting proteins, add labels and attributes from the IntOGen
# database of cancer driver genes, and then illustrate various approaches to
# label and attribute propagation.

# The basic question this example aims to illustrate is: which genes are a part
# of a subnetwork that is significantly enriched in cancer driver genes.
# However, the approach is general, and you can readily substitute many other
# properties of ineterst for the "driver genes" we use here.

# ==== A RANDOM, CONNECTED SUBGRAPH ============================================

# Our first task is to define a subnetwork that we can work with (the complete
# STRING network is too large for demonstration purposes). Since propagation
# requires edges, we need a fully connected subgraph. And since we have no
# particular ideas about the toplogy our subgraph should have, we'll just
# proceed randomly. Here is a function that will run a random walk on a graph G,
# collect all vertices it encounters, and return the subgraph of order N once it
# has encountered N distinct vertices. Note that this is not strictly speaking a
# random subgraph, because we are more likely to explore the neighbourhood of
# highly connected subgraphs, or cliques, more exhaustively than other regions
# of the network - the subgraph topology and degree distribution is therefore
# not random, but depends on the original network topology.

pickSub <- function(G, N, mySeed = 112358) {
  # Pick a random, connected subgraph of N nodes on G.
  # G    an igraph graph.
  # N    integer
  # Value:
  #    igraph graph: induced subgraph on G
  # Details:
  #    A random walk is performed on G from a random starting vertex and its
  #    trajectory is stored in tra, until the number of vertices in tra is N, or
  #    the number of steps exceeds nMAX. Note that the number of vertices in a
  #    trajectory is not the same as number of steps, since vertices can be
  #    encountered multiple times! The subgraph containing all vertices in tra
  #    and all of their mutual edges is returned.
  nMAX <- 10 * N
  set.seed(mySeed)
  tra <- sample(as.numeric(V(G)), 1)  # random starting vertex
  nSteps <- 1

  while (length(unique(tra)) < N && nSteps < nMAX) {
    nei <- as.numeric(neighbors(G, V(G)[tra[1]])) # neighbors of most recent
                                                  # vertex in tra
    if (length(nei) > 1) {    # if more than one neighbor ...
      nei <- sample(nei, 1)   # ... pick only one.
    }
    tra <- c(nei, tra)   # add the neighbor to the front of the trajectory
    nSteps <- nSteps + 1
  }
  return(induced_subgraph(G, V(G)[unique(tra)]))
}

N <- 100  # number of vertices in our random graph
G <- pickSub(gSTR, N, mySeed = 31416)
# ... other interesting graphs you could explore
# G <- pickSub(gSTR, N, mySeed = 11235813)
# G <- pickSub(gSTR, N, mySeed = 27)
# G <- pickSub(gSTR, N, mySeed = 272)
# G <- pickSub(gSTR, N, mySeed = 271828183)
# G <- pickSub(gSTR, N, mySeed = 31)

# In order to make sense of our network, we should translate the ENSP IDs to
# gene symbols. Unfortunately not all of STRINGs IDs can be found in biomart -
# the databases are not synchronized and when Ensembl updates identifiers that
# does not necessarily mean the STRING identifiers get updated as well. However,
# the STRING database has a mapping file of ENSP IDs to UniProt IDs, and the
# latter can be used with the UniProt ID mapping service to map to HGNC gene
# symbols. For details, install and work through the project at
# https://www.github.com/hyginn/R_exercise-IDmapping Here we simply use the
# resulting mapping vector:
load("data/ID2symMap.RData")
head(ID2symMap)

# cheating :-)     (Two ENSP IDs in G are not in ID2symMap. I looked them up.)
x <- c("UBC", "DVL1")
names(x) <- c("ENSP00000344818", "ENSP00000368169")
ID2symMap <- c(ID2symMap, x)

# Change the vertex $name attribute in G
nV <- names(V(G))            # current vertex names
sV <- ID2symMap[nV]          # map to HGNC symbols
na <- which(is.na(sV))       # index of unmapped names
sV[na] <- nV[na]             # use original ENSP IDs instead
G <- set_vertex_attr(G, "name", value = sV)

# Compute a 2-D layout of the graph, for plotting
Gxy <- layout_with_fr(G); a <- 8; b <- 1.5

# Other possible layouts ...
# Gxy <- layout_with_graphopt(G, niter=10000, charge = 0.1); a <- 1500; b <- 300
# Gxy <- layout_with_kk(G); a <- 15; b <- 0.4
# Gxy <- layout_with_mds(G); a <- 8; b <- 0.3

# Centre the graph, to simplify scaling the plot frame to be a bit larger than
# the graph.
Gxy[,1] <- Gxy[,1] - mean(Gxy[,1])
Gxy[,2] <- Gxy[,2] - mean(Gxy[,2])

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(G,
     layout = Gxy, # calculated layout coordinates
     rescale = FALSE,
     xlim = c(min(Gxy[,1]), max(Gxy[,1])) * 1.1,
     ylim = c(min(Gxy[,2]), max(Gxy[,2])) * 1.1,
     vertex.color=heat.colors(max(degree(G)+1))[degree(G)+1],
     vertex.size = a + (b * degree(G)),
     vertex.label = "",
     edge.arrow.size = 0)
text(Gxy, names(V(G)), cex = 0.5)
par(oPar)

# Note that our random walk has picked up 27 of the 31 protein components of the
# mediator complex:
names(V(G))[grep("^MED|CCNC|CDK8", names(V(G)))]
# Physical complexes are highly connected thus random walks have a high
# probability to find all members, if they encounter one member. That's the
# principle behind the map equation graph clustering algorithm.

# Note also that one of the highly connected nodes is UBC. Nodes such as this
# illustrate a conceptual limitation of all network propagation algorithms: the
# fact that two nodes are connected to a high-degree network "hub" such as
# ubiqutin does not _necessarily_ mean they are functionally related.

# Also note that the above statements about the biological interpretation of
# a network are impossible to make if you don't annotate vertices with the
# name/function of the proteins they represent. Gene symbols are rather minimal
# information in that regard, but ENSP IDs or UniProt IDs are useless.


# ==== ANNOTATING CANCER DRIVER GENES ==========================================

# Are there any cancer driver genes in our graph? The IntOGen database publishes
# a convenient list at http://www.intogen.org - I have placed a table into the
# data directory:

intogenDrivers <- read.delim("data/intogen-drivers-data.2017-02-04.tsv",
                             header = TRUE,
                             stringsAsFactors = FALSE)
head(intogenDrivers)

# Are there driver genes in G? (Analyze the expression!)
V(G)[which(names(V(G)) %in% intogenDrivers$SYMBOL)]

# Indeed.

# Lets add an attribute to the graph in which we label "driver" genes:
G <- set_vertex_attr(G, "driver",
                     value = (names(V(G)) %in% intogenDrivers$SYMBOL) )

# Confirm:
head(V(G)$driver, 10)
sum(V(G)$driver)

# Let's plot the driver genes in the graph. vcount() returns the number of
# vertices in an igraph object (the order of the graph).

# For convenience, we define a coloring function that returns colors for values
# from 0 to 1 - we'll need this a lot later. Note how we use apply() and a bound
# to make sure we don't exceed the computed index in our vector of color values.
# The function is written to work on single values, as well as vectors:

scoreCol <- function(val, N = 20) {
  # returns a color from a gradient for val between 0 and 1 discretized
  # in N steps
  val <- apply(cbind(val, 0.00001), 1, max) # bound just above 0
  val <- apply(cbind(val, 0.99999), 1, min) # bound just below 1
  vCol <- colorRampPalette(c("#FFFFDD",
                             "#BBBBC5",
                             "#DD0000"),
                           bias = 2.5)(N)
  return(vCol[floor(val * N) + 1])
}

# display ...
barplot(rep(1, 20), col = scoreCol(seq(0, 1, length.out = 20)),
        axes = FALSE)


# Next we define a convenient plotting function:
plotG <- function(G, Gxy, myVcol = "#FFFFDD", myEcol = "#AAAAAA") {
  oMar <- par("mar")
  par(mar = rep(0,4)) # Turn margins off
  plot(G,
       layout = Gxy, # calculated layout coordinates
       rescale = FALSE,
       xlim = c(min(Gxy[,1]), max(Gxy[,1])) * 1.1,
       ylim = c(min(Gxy[,2]), max(Gxy[,2])) * 1.1,
       vertex.color = myVcol,
       edge.color = myEcol,
       vertex.size = 40,
       vertex.label = "",
       edge.arrow.size = 0)
  text(Gxy, names(V(G)), col = "#00000088", cex = 0.5)
  par(mar = oMar)
}

# Now let's look at our driver genes:

plotG(G, Gxy, scoreCol(as.numeric(V(G)$driver)))

# ===== Label propagation by majority rule

# Perhaps the simplest method for label propagation is "majority rule". A label
# is propagated to an unlabelled node, if the majority of its neighbors have
# that label. A quick look at our network shows however that none of the nodes
# have a majority of driver gene neighbors.

# ===== Label propagation by random walk

# An alternative idea is to determine a label of an unlabelled node by doing a
# random walk until we encounter a labelled node, then counting how often we
# have encountered which label. We would say: the walk has been "absorbed" by a
# label. Let's seed our graph with a few random labels from two categories and
# try such absorptive random walks.

set.seed(11235)
# randomly assign six nodes each to a category 1 and 2
cat1 <- logical(vcount(G))
cat1[sample(1:vcount(G), 6)] <- TRUE
cat2 <- logical(vcount(G))
cat2[sample(1:vcount(G), 6)] <- TRUE

twoLabelCol <- rep("#FFFFDD", vcount(G))
twoLabelCol[cat1] <- "#CC0088"   # Label 1
twoLabelCol[cat2] <- "#0088CC"   # Label 2
plotG(G, Gxy, twoLabelCol)

nWalks = 21
for (i in 1:vcount(G)) {  # for each node in G
  pBar(i, vcount(G))
  if (! cat1[i] && ! cat2[i]) {  # node is not in category 1 or 2
    n1 <- 0
    n2 <- 0
    for (j in 1:nWalks) {
      v <- i   # initialize current node
      while(! cat1[v] && ! cat2[v]) { # continue as long as we don't hit a
                                      # labelled node
          nei <- as.numeric(neighbors(G, V(G)[v])) # neighbors of current node
          if (length(nei) > 1) {    # if more than one neighbor ...
            nei <- sample(nei, 1)   # ... pick only one.
          }
          v <- nei                  # update current node
      } # end while: we have hit a label
      if (cat1[v]) { # count it
        n1 <- n1 + 1
      } else {
        n2 <- n2 + 1
      }
    } # end for nWalks: decide
    if (n1 > n2) {
      twoLabelCol[i] <- "#CC008833" # transparent color for cat1
    } else {
      twoLabelCol[i] <- "#0088CC33"
    }
  }
}

# Let's see how our labels have propagated ...
plotG(G, Gxy, twoLabelCol)

# For a graph bisection (or multisection) problem, this would be a very
# reasonable outcome.

# But we don't have two labels - we have only one kind of label: our driver
# genes. So we need to figure out a different approach. For example, we could
# say that at every step in which we do not encounter a driver node, the random
# walk could be absorbed with probability p by the background. Here is an
# implementation - with a bit of a speedup: rather than compute the expensive
# neighborhood calculation every time, we precompute neighborhoods for each node
# in G once, and store them in a list:

catD <- V(G)$driver
driverLabelCol <- scoreCol(as.numeric(catD))
plotG(G, Gxy, driverLabelCol)

neiList <- list()
for (i in 1:vcount(G)) {
  pBar(i, vcount(G))
  neiList[[i]] <- as.numeric(neighbors(G, V(G)[i]))
}


set.seed(11235)
nWalks = 21
p <- 0.1
for (i in 1:vcount(G)) {  # for each node in G
  pBar(i, vcount(G))
  if (! catD[i]) {  # node is not a driver
    nD  <- 0
    nBg <- 0
    for (j in 1:nWalks) {
      v <- i   # initialize current node
      while(! catD[v] && runif(1) > p) { # continue as long as we don't hit a
                                         # driver node or are absorbed by
                                         # background
        nei <- neiList[[i]]       # neighbors of current node
        if (length(nei) > 1) {    # if more than one neighbor ...
          nei <- sample(nei, 1)   # ... pick only one.
        }
        v <- nei                  # update current node
      } # end while: we have been absorbed
      if (catD[v]) { # count it
        nD <- nD + 1
      } else {
        nBg <- nBg + 1
      }
    } # end for nWalks: decide
    if (nD > nBg) {
      driverLabelCol[i] <- "#DD000033" # transparent color for driver-related
    }
  }
}

# How have the labels propagated?
plotG(G, Gxy, driverLabelCol)

# Again - not unreasonable for a representation of a locally increased density
# of driver genes. If you look carefully, you can find nodes that are neighbors
# to driver genes, but _not_ labelled, because they are embedded in a high
# density of background.

# However, not all driver genes have the same importance and in order to
# propagate real-valued scores we need to look at attribute propagation.


# ==== ATTRIBUTE PROPAGATION ===================================================

# The common way to model attribute propagation is by diffusion, and diffusion
# can by modeled by a simple random walk. Assume we have an attribute valued as
# 1.0 on a particular vertex. Let's call this attribute "heat" to make it easier
# to have some intuition and a mental model about what we are doing. Diffusion,
# modelled as a random walk, means dividing that heat into a number of portions,
# walking for a defined number of steps (that number would represent the "time"
# of a physical diffusion process), and depositing the heat at the node where we
# end up. Let's illustrate this on our graph. Before we do this, we need to
# modify our neighbour list to include every node itself in its neighborhood,
# i.e. our network gets a self-edge at every node. Why? If we would not do this,
# in a walk of one step, _all_ of the heat would have to move to the node's
# neighbour(s).

for (i in 1:length(neiList)) {
  neiList[[i]] <- c(i, neiList[[i]])
}

# Here is a function to "heat" a vertex/vertices v of G, diffuse the heat over t timesteps, and return the result.

diffuseH <- function(v, t, heat, N = 100) {
  # neighList and G must exist in global environment.
  # Arguments:
  #     v:    indices of vertices to diffuse heat from
  #     t:    length of each random walk ("time")
  #     h:    initial heat vector, will be initialized to heat of 1
  #              for each v if missing.
  #     N:    N * t is the number of walks, smoothing the stochastic result
  #
  # Note: since we know that the lenght of the neighbor list of each vertex
  #       is always greater than one, we can simplify the code to select
  #       our next step from our previous walking code.
  # Value:
  #     heat vector for all vertices of G, scaled to set the maximum value
  #     to 1.0

  if (missing(heat)) {
    # initialize heat vector with heat of 1.0 on each element of v
    heat <- numeric(vcount(G))  # initialize
    heat[v] <- 1.0              # add heat
  }
  nWalks <- N * t
  for (thisV in v) {                   # do for each node in v
    dH <- heat[thisV] / nWalks         # size of heat portion
    for (i in 1:nWalks) {
      heat[thisV] <- heat[thisV] - dH  # pick up heat portion
      currV <- thisV                   # start at thisV
      for (j in 1:t) {                 # walk for t random steps
        currV <- sample(neiList[[currV]], 1)
      } # end walking
      heat[currV] <- heat[currV] + dH # deposit heat portion
    }
  }
  heat <- heat / max(heat)  # rescale for better visibility
  return(heat)
}


plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 1)))
plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 2)))
plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 3)))
plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 4)))
plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 5)))
plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 6)))
plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 7)))
plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 8)))
plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 9)))

# Notice that as the length of the walk increases, heat gets concentrated in
# densely connected regions - this is where the random walk easily gets
# "trapped".

plotG(G, Gxy, scoreCol(diffuseH(v = 3, t = 30)))


# Here is what happens when we have multiple nodes - our driver genes:
plotG(G, Gxy, scoreCol(diffuseH(v = which(V(G)$driver), t = 1)))
plotG(G, Gxy, scoreCol(diffuseH(v = which(V(G)$driver), t = 2)))
plotG(G, Gxy, scoreCol(diffuseH(v = which(V(G)$driver), t = 3)))



# ==== OncodriveFM scores as "heat"

# We know that not all driver gene mutations are equally significant - i.e. the
# amount of "heat" on each gene is not the same. For example, we can use the
# Oncodrive score in the IntOGen table as "heat":


head(intogenDrivers)
# "CS" mutations are coding sequence mutations, "PAM" are protein affecting
# mutations. Signal types are C (Clustered), F (Functional), and R (Recurrent
# mutations). For Oncodrive scores, see here:
# http://bg.upf.edu/group/projects/oncodrive-fm.php


# === Preparing attributes

# Let's look at the score distributions, and what scores we find in the genes in
# our network:

# get an ordering vector, that orders the intogen table by Oncodrive scores
ord <- order(intogenDrivers$ONCODRIVE_ROLE, decreasing = TRUE)

# get a color vector for all driver genes
myCol <- scoreCol(intogenDrivers$ONCODRIVE_ROLE[ord])

# plot all driver gene scores
plot(intogenDrivers$ONCODRIVE_ROLE[ord],
     type = "l", ylab = "OncodriveFM score",
     main = "Driver scores (all, and in G)")

# add a color scale to the y-axis
N <- 20
yBreaks <- seq(0, 1, length.out = N + 1)
for (i in 1:N) {
  rect(par("usr")[1],
       yBreaks[i],
       par("usr")[1] / 2,
       yBreaks[i + 1],
       col = scoreCol(mean(c(yBreaks[i], yBreaks[i+1]))),
       border = NA)
}

# collect the driver gene names
namDG <- V(G)$name[V(G)$driver]

# add the driver gene labels to the plot
for (i in 1:nrow(intogenDrivers)) {
  if ((intogenDrivers$SYMBOL[ord])[i] %in% namDG) {
    # gene symbol
    text(i,
         (intogenDrivers$ONCODRIVE_ROLE[ord])[i],
         labels = (intogenDrivers$SYMBOL[ord])[i],
         cex = 0.7)
    # horizontal line
    points(c(par("usr")[1], i),
           rep((intogenDrivers$ONCODRIVE_ROLE[ord])[i], 2),
           type = "l",
           col = myCol[i],
           lwd = 0.7)
  }
}

# This plot shows where the scores for driver genes in G fall in the global
# distribution of scores. Let's add the scores as attributes to our graph.

gScores <- rep(0.0, vcount(G))  # vector of zero scores

for (i in which(V(G)$driver)) {
  # overwrite gScores[i] with driver gene score
  dName <- V(G)$name[i]
  iInto <- which(intogenDrivers$SYMBOL == dName)
  gScores[i] <- intogenDrivers$ONCODRIVE_ROLE[iInto]
}
G <- set_vertex_attr(G, "oncodrive", value = gScores )

# With this, we can replot our network, coloring the nodes with the score-color
# scale we defined above.

plotG(G, Gxy, scoreCol(V(G)$oncodrive))


# === Propagating attributes by diffusion

# Now we can propagate these values ...
hO <- V(G)$oncodrive
vD <- which(V(G)$driver)
plotG(G, Gxy, scoreCol(diffuseH(v = vD, t = 1, heat = hO)))
plotG(G, Gxy, scoreCol(diffuseH(v = vD, t = 2, heat = hO)))
plotG(G, Gxy, scoreCol(diffuseH(v = vD, t = 3, heat = hO)))
plotG(G, Gxy, scoreCol(diffuseH(v = vD, t = 4, heat = hO)))
plotG(G, Gxy, scoreCol(diffuseH(v = vD, t = 5, heat = hO)))

plotG(G, Gxy, scoreCol(diffuseH(v = vD, t = 20, heat = hO)))

# What difference between unit-heat and score based heat?
plotG(G, Gxy, scoreCol(diffuseH(v = vD, t = 3)))             # unit heat
plotG(G, Gxy, scoreCol(diffuseH(v = vD, t = 3, heat = hO)))  # scores


# As you see, there are differences in the details of distribution, which
# presumably reflect the underlying biology better. Regarding the global
# distribution we see that diffusing for a few steps identifies related nodes,
# modulated by the topology of the network. But diffusing for too many steps
# concentrates most of the heat in the densely connected clusters and in the
# nodes with high betweeness centrality ("hubs"). For short times, the heat
# distribution is dominated by the initial heat at time zero. For long times,
# the heat distribution is dominated by the topology of the network. We need to
# keep this in mind when we use attribute diffusion to define subnetworks.



# ==============================================================================
#        PART SEVEN: ATTRIBUTE DIFFUSION, EDGE WEIGHTS AND SUBNETWORKS
# ==============================================================================

# Are there subnetworks of our graph that are important for the cancer
# phenotype?

# If we can find subnetworks in our graph that point to some functional
# relationship, AND find cancer driver genes to be enriched in such a
# subnetwork, then we could conclude that that entire set of genes somehow
# collaborates towards an important function. But we can't enumerate all
# subnetworks - there are far too many -, and even if we did, all statistical
# power for enrichment would be lost because of the huge multiple testing
# problem we would then be confronted with. Therefore, it is of interest to
# discover only subnetworks that have some biological significance. And for this
# purpose, we can use attribute diffusion. In particular: if we can use
# attribute diffusion to compute edge-weights, and then remove all edges below a
# suitable threshold, the network will naturally fall apart into components,
# which we can then further analyze as subnetworks. The components will then
# depend on both the network topology, and the "heat" scores which initally
# label driver genes.

# There are certainly many ways to do this - here we'll discuss the approach
# taken by the authors of HotNet and HotNet2.


# ==== THE HotNet APPROACH =====================================================

# In Hotnet:
#  - First an edge "influence" is determined based on the amount of "heat"
#    that flows along the edge if unit heat is placed on either of its vertices.
#  - Second, an edge weight is determined, by multiplying the influence with
#    the score of either vertex.
#  - Finally, edges whose weight falls below a threshold are removed.

# Intuitively, you can imagine that this process is a mix of the short- and long
# walk regimes we explored above. In particular, it is NOT true that at
# equilibrium all genes have the same heat, as Vandin et al. claim. (2012, 59).
# Rather the heat in the network redistributes into densely connected clusters
# and nodes with a high betweenness centrality.

# You will find that the literature computes such heat diffusion influences via
# a "heat kernel" which is the matrix exponential of the Laplacian matrix of the
# graph. This procedure assumes unweighted, undirected graphs. Therefore we will
# continue simply estimating the required values from random walks, which can
# trivially be adapted to weighted, directed graphs with all manners of
# additional enhancements.

# We could query and update edge attributes within iGraph, but for efficiency
# and explicitness we will work on the adjacency matrix in what follows. That's
# no problem since our demo-network is tiny with 100 vertices and 1,010 edges.
# Just as an aside, note that in general, for genome scale networks with OTO
# 10e04 vertices and 10e05 edges, a regular adjacency matrix becomes too large.
# However this is a _very_ sparse matrix, as the matrix grows with |V|^2 but |E|
# only grows with k*|V| where k is on the order of 10. Therefore we need to work
# with sparse matrix data structures - from the R "matrix" package, as supported
# by igraph, whose size scales with |E|.

# The adjacency matrix:
A <- as_adjacency_matrix(G, type = "both")
str(A)

# print one row ...
A[which(names(V(G)) == "EIF4A3"), ]

# Degree of this node:
sum(A[which(names(V(G)) == "EIF4A3"), ] == 1)

# (You can count to confirm ...)

# ==== Calculating "influence"

# If we take "influence" between two vertices (u, v) to be the amount of heat
# that arrives on a vertex v after a time t if unit heat is placed on vertex u.
# Note that we are only interested in (u, v) that have an edge between them. So we can simply run a number of random walks of length t for each vertex u, count the number of times we end on a neighbor of u, and divide this by the total number of walks.

# We need a neighbour list  without self-edges
for (i in 1:vcount(G)) {
  neiList[[i]] <- as.numeric(neighbors(G, V(G)[i]))
}

# Let's look at neighbors of vertex 62 (UBC):
A[62, ]
sum(A[62, ])

myE <- numeric(vcount(G))
thisV <- 62
t <- 7
nWalks <- 100
for (i in 1:nWalks) {
  currV <- thisV              # start at thisV
  for (j in 1:t) {            # walk for t random steps
    nei <- neiList[[currV]]   # neighbors of current node
    if (length(nei) > 1) {    # if more than one neighbour ...
      nei <- sample(nei, 1)   # ... pick only one.
    }
    currV <- nei                  # update current node
  } # end walking
  myE[currV] <- myE[currV] + 1
}

# Result:
head(myE, 10)

# Scale:
myE <- myE / nWalks
head(myE, 10)

# Keep only values of neighbour edges by multypling myE with the row of the
# adjacency matrix:
myE <- myE * A[thisV, ]

# Plot the result:

# Vertex colors:
vc <- rep("#FFFFDD", vcount(G))
vc[thisV] <- "#DD0000"

# make a vector of vertex pairs to select edges. (see: ?E)
nei <- which(as.logical(A[thisV, ]))
( myP <- c(rbind(thisV, nei)) ) # analyze this!
# add the reverse edges too
( myP <- c(myP, c(rbind(nei, thisV))) ) # analyze this too!

# baseline color ...
E(G)$influence <- 0.05
# Influences (scaled)
E(G, P = myP)$influence <- c(myE[nei], myE[nei]) / max(myE[nei])
# map to to color values
ec <- scoreCol(E(G)$influence)

plotG(G, Gxy, myVcol = vc, myEcol = ec)

# Note how the influences are weaker on nodes that have a high probability that
# paths passing through them will not return.

# To calculate all influences, we simly repeat what we have done above in a
# loop. We store the values in an adjacency matrix H.

H <- A

t <- 5
nWalks <- 100
for (i in 1:vcount(G)){
  pBar(i, vcount(G))
  thisV <- i
  H[thisV, ] <- numeric(vcount(G))
  for (i in 1:nWalks) {
    currV <- thisV              # start at thisV
    for (j in 1:t) {            # walk for t random steps
      nei <- neiList[[currV]]   # neighbors of current node
      if (length(nei) > 1) {    # if more than one neighbour ...
        nei <- sample(nei, 1)   # ... pick only one.
      }
      currV <- nei                  # update current node
    } # end walking
    H[thisV, currV] <- H[thisV,currV] + 1
  }
  H[thisV, ] <- H[thisV, ] / nWalks       # scale
  H[thisV, ] <- H[thisV, ] * A[thisV, ]   # remove counts for non-neighbours
}

# ... a subset of the results.
H[1:15, 1:15]

# Note that the matrix is not symmetric and H(u, v) is usually not the same as
# H(v, u). Vandin et al. (2012) make the adjacency matrix symmetric by selecting
# min(H(u, v), H(v, u)) for both edges, ... but that seems debatable to me.
# Let's reproduce this anyway.

for (u in 1:(vcount(G) - 1)) {
  for (v in (u + 1):vcount(G)) {
    if (H[u, v] != 0 || H[v, u] != 0) {
      h <- min(H[u, v], H[v, u])
      H[u, v] <- h
      H[v, u] <- h
    }
  }
}

H[1:15, 1:15]
max(H)

# Next we update the edge weights by taking scores into account. HotNet simply
# multiplies the values from H with the max of the scores of the two vertices
# that share the edge. We have stored the OncodriveFM stores in:
head(V(G)$oncodrive)

W <- H

for (u in 1:(vcount(G) - 1)) {
  for (v in (u + 1):vcount(G)) {
    if (W[u, v] != 0) {
      s <- max(V(G)$oncodrive[c(u, v)]) * H[u, v]
      W[u, v] <- s
      W[v, u] <- s
    }
  }
}

# for display purpose, we rescale W:
W <- W / max(W)

# Done. Need to convert this to a list of edge colors to plot. Just a bit of bookkeeping ...

# baseline color ...
E(G)$weight <- 0.03

# build a vector of vertices that share edges, and a vector of weights

vP <- numeric()
vW <- numeric()

for (u in 1:(vcount(G) - 1)) {
  for (v in (u + 1):vcount(G)) {
    if (W[u, v] != 0) {
      vP <- c(vP, u, v, v, u)        # vertex pairs
      vW <- c(vW, W[u, v], W[u, v])  # weights
    }
  }
}

# set edge attributes
E(G, P = vP)$weight <- vW

# map to to color values
ec <- scoreCol(unlist(E(G)$weight), N = 100)

plotG(G, Gxy, myVcol = vc, myEcol = ec)

# Easy to see that we could now delete all vertices that don't, say, share an edge with a weight >= 0.2 ...

G2 <- G                                      # make a copy
sel <- which(E(G2)$weight < 0.2)             # select edges below threshold ...
G2 <- delete_edges(G2, sel)                  # ... and delete them.
sel <- which(degree(G2) == 0)                # select disconnected vertices ...
G2 <- delete_vertices(G2, sel)               # ... and delete them
G2xy <- Gxy[-sel, ]                          # ... and delete their layout x, y
vc <- scoreCol(V(G)$oncodrive)               # update vertex colors
ec <- scoreCol(unlist(E(G)$weight), N = 100) # update edge colors

plotG(G2, G2xy, myVcol = vc, myEcol = ec)

# Done ... here are our subnetworks.
components(G2)


# ==== THE HotNet2 APPROACH ====================================================


# TBC




# ==============================================================================
#        PART EIGHT: OUTLOOK
# ==============================================================================

# There are many more functions for graph and network analysis that this
# tutorial did not cover. You should know about the following. Look up the
# function and write a short bit of example code that uses it:"

?motifs           # to find network motifs
?neighbors        # to find neighbors
?neighborhood     # to find nodes that are not further than a given
                  #    distance from a center


# Then you should know about the following packages. There is an extensive set
# of biological graph algorithms in bioconductor packages:
#
# http://bioconductor.org/packages/release/BiocViews.html#___GraphAndNetwork
#
# ...and definitely check out MCL:
#
# https://cran.r-project.org/web/packages/MCL/index.html



# [END]
