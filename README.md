# CODAS
Clustering of Online Data Streams into Arbitrary Shapes

Hyde, R.; Angelov, P., "A New Local Density Based Approach for Clustering Online Data in Arbitrary Shapes Clusters," CYBCONF2015
doi: 10.1109/CYBConf.2015.7175937
Downloadable from: http://ieeexplore.ieee.org/stamp/stamp.jsp?arnumber=7175937

This algorithm deals with the clustering of online, continuous data-streams. This is distinct from the CEDAS algorithm, https://rhyde67.github.io/CEDAS/ which deals with evolving data streams. There is often confusion between these two, distinct types of data. Online data refers to data that arrives consecutively over time. Algorithms to deal with these must react to each data sample and then discard or archive the data. The nature of these streams of data allow for the clusters to migrate or change size and possibly shape over time. However they do not 'evolve', i.e. the clusters do not die out, or be re-born over time.

To illustrate the difference between online and evolving data consider two cluster A and B at time 't' occupying two different regions in data space. If, after time t+1000 cluster A has moved to the data space previously occupied by cluster B, and cluster B has moved to A, then these clusters will have merged. This limits the usefulness of online clustering for time series data to 'short term' data streams, where 'short term' is relative to the time dependent changes in the data.

In contrast, 'evolving' clustering will allow for this and consider only 'recent' data, such that cluster A and B will be separate at time t+1000 as they occupy different data space.

Additional confusion is added by the use of the terminology 'evolving clustering' to refer to algporithms where clusters change their parameters e.g. radii and centre of time. In most examples of 'evolving clustering' there is no geuine evolution of the clusters, rather a single generational adaptation to variations in relatively stationary data clusters.

The simple test is to see if the clusters have some method of birth, decay and death. If not then it is 'online', if it does then it is 'evolving'.

##Files:

CODAS2_ver02.m - Matlab code for the CODAS algorithm.

CODAS2_Test.m - Matlab code to run the algorithm across the various test data sets. Suggested parameters are included.

.csv file contain sample datasets

circles.m - not my work, Matlab script for displaying circles when plotting the clustering results. Downloadable from Matlab's file exchange.

distinguishable_colors.m - not my work, Matlab script for creating easily distiguishable colours. Downloadable from Matlab's file exchange.
