# R_Exercise-Graphs.R
#
# Purpose:  Exercises for working with graphs in R.
#
# Version: 1.0
#
# Date:    2016  02  28
# Author:  Boris Steipe
#
# V 1.0    First code
#
# TODO:
#
#
# == HOW TO WORK WITH THIS FILE ======================================
#
#  Go through this script line by line to read and understand the
#  code. Execute code by selecting it and typing <cmd><enter>. Edit it
#  as required, experiment and play. DO NOT simply source() this whole
#  file! If there are portions you don't understand, use R's help
#  system, google for an answer, or ask me. Don't continue if you
#  don't understand what's going on. That's not how it works ...
#
# ====================================================================

" This tutorial covers basic concepts of graph theory and analysis in R. You should have typed init() to configure some utilities in the background.

"


# == INITIALIZATION ==================================================

source("utilities.R") # load a few convenience functions


# ====================================================================
#        PART ONE: REVIEW
# ====================================================================
"I assume you'll have read the Pavlopoulos review of graph theory concepts. Let's explore some of the ideas by starting with a small random graph."

set.seed(112358)

N <- 20
# create a vector of N random node names
Nnames <- ABC.rSyll(N)
Nnames

# create an edge matrix with a given node-edge probability.

G <- matrix(numeric(N * N), ncol = N)
rownames(G) <- Nnames
colnames(G) <- Nnames
pEdge <- 0.06
for (iRow in 1:N-1) {
  for (iCol in (iRow+1):N) {
    if (runif(1) < pEdge) {
      G[iRow, iCol] <- 1
      G[iCol, iRow] <- 1
    }
  }
}
G

ABC.fig("randGraphSample.png") # We'll learn to plot graphs later ...

"The simplest descriptor of a graph are the number of nodes, edges, and the degree-distribution. In our example, the number of nodes was given: N; the number of edges can easily be calculated from the adjacency matrix. In our matrix, we have entered 1 for every edge. Thus we simply sum over the matrix:"
sum(G)

"Is that correct? Is that what you see in the plot?"

"Yes and no: we entered every edge twice: once for a node [i,j], and again for the combination [j, i]. Whether that is correct depends on what exactly we want to do with the matrix. If these were directed edges, we would need to keep track of them separately. Since we didn't intend them to be directed, we'll just divide the number of edges by 2. Why didn't we simply use an upper-triangular matrix? Because then we need to keep track of the ordering of edges if we want to know whether a particular edge exists or not. For example we could sort the nodes alphabetically, and make sure we always query a pair in alphabetical order. Then a triangular matrix would be efficient.

What about the degree distribution? We can get that simply by summing over the rows (or the columns):"

rowSums(G)

rs <- rowSums(G)
brk <- seq(min(rs)-0.5, max(rs)+0.5, by=1)
hist(rs, breaks=brk, col=ABC.pal(length(brk)))

"The degree distribution is actually quite an important descriptor of graphs, since it is very sensitive to the generating mechanism. For biological networks, that is one of the key questions we are interested in: how was the network formed?

I think you can easily implement measures such as total connectivity or graph density. Many more functions can be found in graph analysis packages, such as the iGraph package that we will use below.
"
# ====================================================================
#        PART TWO: COMPUTING WITH GRAPHS
# ====================================================================

"To get a better feel for computing with graphs, let's write a simple algorithm: counting the number of connected components in the graph. A subgraph is connected if there is a path from any node to any other. As the image above showed you, our graph G has four connected components.

The algorithm employs a simple breadth-first search: we start at a node and initialize our first connected component with it. Then we add all its neigbors to the connected component. Then their neighbors, until there are none left to add. If there are remaining nodes that have not yet been processed, they must be part of another component. We initialize another component and repeat until all nodes have been processed."

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
  # they are 0, or 0.0) and TRUE if they do represent an edge (any non-zero number).
  # Thus this will work for weighted and unweighted graphs.
  #                                adjM[node, ]   contents of node's row of the matrix,
  #                     as.logical(adjM[node, ])  ... turned into a logical vector,
  #      colnames(adjM)[as.logical(adjM[node, ])] ... used to subset the return-value of colnames()
  return(colnames(adjM)[as.logical(adjM[node, ])])
}


for (node in rownames(G)) {                 # Process each node.
  if (! processed[node]) {                  # If we haven't processed it previously:
    iC <- length(components) + 1            #   initialize a new component
    components[[iC]] <- node                #   ... with this node.
    queue <- neighboursOf(node, G)          #   Append its neighbors to the queue
    processed[node] <- TRUE                 #   Flag it as processed.
    while (length(queue) > 0) {             #   Now work through the queue:
      currentNode <- queue[1]               #     get the node from position one,
      queue <- queue[-1]                    #     ... and remove it from the queue.
      if (! processed[currentNode]) {       #     If we haven't processed it previously:
        components[[iC]] <- c(components[[iC]], currentNode) # ... add it to the current component,
        queue <- c(queue, neighboursOf(currentNode, G)) # ... add its neigbors to the queue
        processed[currentNode] <- TRUE      #        ... and flag it as processed.
      }
    }                                       #   Repeat until the queue is empty.
  }
}

"A brief note on using queues in R:

R is a functional language. In clean R code, functions do not have side effects. In my experience this really makes a lot of difference for writing bug-free code that is easy to read. But sometimes this means we have to take extra steps. For example, some languages have functions to push() and pop() elements on/off stacks - or shift()/unshift(), or enqueue(), dequeue() elements from a queue. But a dequeue operation does two things: it retrieves a value AND it modifies the queue. In R, this normally requires two statements - that's the way we wrote it above. That's fine because we really DO make two changes to our data. But if you feel you MUST have a function that does this in a single expression, this can be done. See:
   https://gist.github.com/leeper/d6d085ac86d1e006167e  and
   https://gist.github.com/mpettis/b7bfeff282e3b052684f also
   https://www.researchgate.net/post/What_is_the_queue_data_structure_in_R
Of course you could also hold a value and an array in a list that you modify with a single function. So there are ways to wrap this into one expression. But knowing that you CAN do this, doesn't mean you SHOULD.

Now let's look at the results:
"

unlist(lapply(components, length))  # size of each component
sizes <- unlist(lapply(components, length))
length(sizes)                       # number of components
mean(sizes)                         # average size
sd(sizes)                           # sd of size
components                          # the actual components

"So far so good, and if you want to practice coding, you could perhaps write your own implementation of Dijkstra's algorithm to find shortest paths (https://en.wikipedia.org/wiki/Dijkstra's_algorithm). But of course most graph algorithms of interest are available in an R package. There are several graph packages in R, but currently the package of choice seems to be iGraph."

# ====================================================================
#        PART THREE: THE igraph PACKAGE
# ====================================================================

if (!require(igraph)) {
  install.packages("igraph")
  library(igraph)
}

# == BASICS ==========================================================

"The basic object of the igraph package is a graph object. Let's turn our adjecency matrix into such an object and explore:"
iG <- graph_from_adjacency_matrix(G)
summary(iG)

"Have a look at what the output IGRAPH DN-- 20 34 -- means: study the help page"
?print_igraph

mode(iG)
class(iG)

"So an igraph graph object is a special list object; it is opaque in the sense that a user is never expected to modify its components directly, but through a variety of helper functions the package provides. There are many ways to construct graphs - from adjacency matrices, as we have just done, from edge lists, or by producing random graphs according to a variety of recipes, called _games_ in this package.

As with many R objects, loading the package provides special functions that can be accessed via the same name as the basic R functions, for example:
"
print(iG)
plot(iG)

"... where plot() allows the usual flexibility of fine-tuning the plot. We first layout the node coordinates with the Fruchtermann-Reingold algorithm - a force-directed layout that applies an ettractive potential along edges (which oulls nodes together) and a repulsive potential to nodes (so they don't overlap). Note the use of the degree() function to color and scale nodes and labels by degree and the use of the V() function to retrieve the vertex names. See ?plot.igraph for details."

iGxy <- layout_with_fr(iG)   # calculate layout coordinates

# Plot with some customizing parameters
oPar <- par(mar= rep(0,4)) # Turn margins off
plot(iG,
     layout = iGxy,
     vertex.color=ABC.pal(max(degree(iG)+1))[degree(iG)+1],
     vertex.size = 9 + (2 * degree(iG)),
     vertex.label.cex = 0.5 + (0.05 * degree(iG)),
     edge.arrow.size = 0,
     edge.width = 2,
     vertex.label = toupper(V(iG)$name))
par(oPar)


# == Components

"The igraph function components() works in the same way as the code we wrote above. But its output is a bit different."
components(iG)

"In the _membership_ vector, nodes are annotetd with the index of the component they are part of. feemp is in component 1, seaks is in the second component etc. This is perhaps more clear if we sort by component index:"
sort(components(iG)$membership)

"Getting e.g. the members of the first component simply can be done by filtering:"

                                 components(iG)$membership == 1  # logical ..
       components(iG)$membership[components(iG)$membership == 1] # extracting the elements
names(components(iG)$membership)[components(iG)$membership == 1] # names of the elements



# == RANDOM GRAPHS AND GRAPH METRICS =================================


"Let's explore some of the more interesting, topological graph measures from the Pavlopoulos paper. We start by building a soemwhat bigger graph. We aren't quite sure whether biological graphs are small-world, or random-geometric, or preferential-attachment ... but igraph has ways to simulate the basic ones (and we could easily simulate our own). Look at the following help pages:"

?sample_gnm                      # see also sample_gnp for the Erdös-Rényi models
?sample_smallworld               # for the Watts & Strogatz model
?sample_pa                       # for the Barabasi-Albert model

"But note that there are many more sample_ functions. Check out the docs!

sample_pa() is a stochastic algorithm for a scale-free graph based on the Barabasi-Albert preferential attachment model. Let's plot a moderately large graph, with a force-directed layout, colored by nmode-degree."


set.seed(27172)
gBA <- sample_pa(1000, power=1.05, out.dist = c(0, 1, 0.3, 0.09, 0.027, 0.0081))
gBAxy <- layout_with_graphopt(gBA, charge=0.001)   # calculate layout coordinates

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(gBA,
     layout = gBAxy,
     rescale = FALSE,
     xlim = c(min(gBAxy[,1]), max(gBAxy[,1])),
     ylim = c(min(gBAxy[,2]), max(gBAxy[,2])),
     vertex.color=ABC.pal(max(degree(gBA)+1))[degree(gBA)+1],
     vertex.size = 800 + (80 * degree(gBA)),
     vertex.label = "",
     edge.arrow.size = 0)
par(oPar)

"The large nodes arise from the fact that a node of high degree is more likely to attract even more additional edges. "


# == Diameter

diameter(gBA)  # The maximum length shortest path
# let's plot the diameter though:
lines(gBAxy[get_diameter(gBA),], lwd=1.5, col="#AA0000")
" ... a good reminder that as nice as the netwoirk plots are, they don't necessarily show us the true topological structure of a graph well."


# == Centralization scores

?centralize
centr_betw(gBA)  # calculate betweenness centrality

"replot, and color by log_betweenness"

nodeBetw <- centr_betw(gBA)$res
nodeBetw <- round(log(nodeBetw +1)) + 1

oPar <- par(mar= rep(0,4)) # Turn margins off
plot(gBA,
     layout = gBAxy,
     rescale = FALSE,
     xlim = c(min(gBAxy[,1]), max(gBAxy[,1])),
     ylim = c(min(gBAxy[,2]), max(gBAxy[,2])),
     vertex.color=heat.colors(max(nodeBetw))[nodeBetw],
     vertex.size = 800 + (80 * degree(gBA)),
     vertex.label = "",
     edge.arrow.size = 0)

par(oPar)

"Note that the betweenness - the number of shortest paths that pass through a node, is in general higher for high-degree nodes - but not always: one of the larger nodes has a very low betweenness: this measure really depends on the detailed local topology of the graph."


# == Degree distribution

"Let's calculate the cumulative degree distribution of our network. The frist element are zero-degree nodes, the second are 0n-degree nodes and so on. Our largest node has degree 141. We remove the first element from this vector, since our network has no zero-degree nodes."

degree_distribution(gBA, cumulative=TRUE)[-1]
dDeg <- degree_distribution(gBA, cumulative=TRUE)[-1]

"Pavlopoulos mentioned three popular ways to plot the degree distribution:

1: log(rank) vs. log(degree):"
plot(log(rank(degree(gBA))), log(degree(gBA)))

"This is not a power-law distribution!

2: histogram of log(degree):"
hist(log(degree(gBA)))

"3. degrees vs. cumulative degree distribution"
plot(log(degree(gBA)), log(dDeg[degree(gBA)]))

"This plot gives a nearly straight line for the preferential attachment model.

While one may attempt to interpret the parameters of these plots - slopes, exponents, intercepts ... they are really most useful as sensitive measures to compare a biological network with a random, synthetic network to genrate hypotheses regarding which TYPE of model could be underlying our observation.
"


# == CLUSTERING ======================================================

'Clustering finds "communities" in graphs - and depending what the edges represent, these could be complexes, pathways, biological systems or similar. There are many graph-clustering algorithms. One approach with many attractive properties is the Map Equation, developed by Martin Rosvall. See: http://www.ncbi.nlm.nih.gov/pubmed/18216267 and htttp://www.mapequation.org '

gBAclusters <- cluster_infomap(gBA)
modularity(gBAclusters)   # ... measures how separated the different membership types are from
                          # each other
membership(gBAclusters)
table(membership(gBAclusters))  # The sizes of the clusters

"Lets plot our graph again, coloring the nodes of the first five communities by their membership:"
commColors <- rep("#f1eef6", max(membership(gBAclusters)))  # as many light colors as we
                                                            # have communities
commColors[1:5] <- c("#980043", "#dd1c77", "#df65b0", "#c994c7", "#d4b9da")


oPar <- par(mar= rep(0,4)) # Turn margins off
plot(gBA,
     layout = gBAxy,
     rescale = FALSE,
     xlim = c(min(gBAxy[,1]), max(gBAxy[,1])),
     ylim = c(min(gBAxy[,2]), max(gBAxy[,2])),
     vertex.color=commColors[membership(gBAclusters)],
     vertex.size = 800 + (80 * degree(gBA)),
     vertex.label = "",
     edge.arrow.size = 0)

par(oPar)



# ====================================================================
#        PART FOUR: PLAY
# ====================================================================

"In order for you to explore some real, biological networks, I give you a dataframe of functional relationships of human proteins downloaded from the STRING database. The full table has 8.5 million records, here is a subset of records with combined scores > 980"

ABC.fig("STRINGscores.png") # STRING scores distribution

"The selected set of edges with a confidence of > 980 is a dataframe with about 50,000 edges and 6,500 unique proteins. You can load the saved dataframe here (and also read more about what the numbers mean at http://www.ncbi.nlm.nih.gov/pubmed/15608232 ).
"

load("STRINGedges.RData")
# make a graph from this dataframe
?graph_from_data_frame

gSTR <- graph_from_data_frame(STRINGedges)

compSTR <- components(gSTR)
summary(compSTR) # our graph is fully connected!


hist(log(degree(gSTR)))
dDegSTR <- degree_distribution(gSTR, cumulative=TRUE)[-1]
plot(log(degree(gSTR)), log(dDegSTR[degree(gSTR)]))

"This is very cool! I haven't seen this before. What does it mean?

Now explore some more:

- Find the largest cliques with largest_cliques(). Pick one of the proteins and find out what it is (you can simply google for it). Is this expected?

- Find the nodes with the 10 - highest betweenness centralities. Use centr_betw() to calculate the values, V() to get the names, and how many nodes there are. The N - 10 highest ranked nodes is what you are looking for. Get the list of IDs. Then find what these proteins are...

... are you surprised? (I am!)

Pick ten random proteins for comparison:
"
sample(V(gSTR), 10 )

"... these seem to be a very different bunch."







# ====================================================================
#        APPENDIX: OUTLOOK
# ====================================================================
"There are many more functions for graph and network analysis that this tutorial did not cover. You should know about the following. Look up the function and write a short bit of example code that uses it:"

?motifs           # to find network motifs
?neighbors        # to find neighbors
?neighborhood     # to find nodes that are not further than a given distance from a center


"Then you should know about the following packages. There is an extensive set of biological graph algorithms in bioconductor packages:

http://bioconductor.org/packages/release/BiocViews.html#___GraphAndNetwork

...and definitely check out MCL:

https://cran.r-project.org/web/packages/MCL/index.html

"




# [END]
