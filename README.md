# QUBIC: QUalitative BIClustering algorithm #
## Algorithm ##

A qualitative biclustering algorithm for analyses of gene expression data.
Biclustering extends the traditional clustering techniques by attempting to find (all) subgroups of genes with similar expression patterns under to-be-identified subsets of experimental conditions when applied to gene expression data. Still the real power of this clustering strategy is yet to be fully realized due to the lack of effective and efficient algorithms for reliably solving the general biclustering problem. We report a QUalitative BIClustering algorithm (QUBIC) that can solve the biclustering problem in a more general form, compared to existing algorithms, through employing a combination of qualitative (or semi-quantitative) measures of gene expression data and a combinatorial optimization technique. One key unique feature of the QUBIC algorithm is that it can identify all statistically significant biclusters including biclusters with the so-called ‘scaling patterns’, a problem considered to be rather challenging; another key unique feature is that the algorithm solves such general biclustering problems very efficiently, capable of solving biclustering problems with tens of thousands of genes under up to thousands of conditions in a few minutes of the CPU time on a desktop computer. We have demonstrated a considerably improved biclustering performance by our algorithm compared to the existing algorithms on various benchmark sets and data sets of our own.

**Citing us**: Li G, Ma Q, Tang H, Paterson AH, Xu Y.
QUBIC: a qualitative biclustering algorithm for analyses of gene expression data.
*Nucleic Acids Research*. 2009;**37(15)**:e101. doi:10.1093/nar/gkp491.

## Data Source ##
The related data sets and results can be found at http://csbl.bmb.uga.edu/~maqin/bicluster/benchmark.html.


## Web-server ##
**QServer**: http://csbl.bmb.uga.edu/publications/materials/ffzhou/QServer/

**Citing us**: Zhou F, Ma Q, Li G, Xu Y.
QServer: A Biclustering Server for Prediction and Assessment of Co-Expressed Gene Clusters.
*PLoS ONE*. 2012;**7(3)**:e32660. doi: 10.1371/journal.pone.0032660


## Evaluation ##
Based on an evaluation paper (see reference below) published at Brief Bioinformatics by _Kemal Eren_ et al., QUBIC handle noise **much better** than the other biclustering programs and has the **highest** enriched bicluster ratio in real data sets.

Eren K, Deveci M, Küçüktunç O, Catalyürek UV, A comparative analysis of biclustering algorithms for gene expression data. **Brief Bioinform.**, 2013.

# Usage #
This software provides a biclustering program for microarray data. For a set of genes and a set of conditions, the program outputs a block-like structure which shows uniform pattern within the block, the block would contain only subsets of all given genes under subsets of all given conditions.

Certain parts of the code uses open-source data structure library codes, including:
  * fib <http://resnet.uoregon.edu/~gurney_j/jmpc/fib.html>, copyright information in fib.c
  * Mark A. Weiss's data structure codes <http://www.cs.fiu.edu/~weiss/>

# Installation #
Simply put "qubic.tar.gz" in any directory, type
```
$ tar zxvf qubic.tar.gz
```
enter the folder '`qubic`' and type '`make`' then the compiled codes are within the same directory as the source.

# Inputs and outputs #
The major program in the provided package is '`qubic`', it can parse two formats of files, discrete data and continuous data, and examples for each are provided in the same folder.

  * To See help and look at all available options.
```
$ ./qubic -h
```
  * Take a look at `toy_example` (discrete data) first. And try to run biclustering
```
$ ./qubic -i toy_example -d
```
`-d` is important here since it tells the program that this is discrete data.

  * Then look at a continuous data, example. Try to run
```
$ ./qubic -i example -f .25
```
This restricts no two output blocks overlap more than 0.25 of the size of each one. And the other parameters are default value.

  * For each input file, our program generates three output files, namely,'`.blocks`' file', '`.chars`'file' and '`.rules`' file'.

In '`.blocks`' file', you can see all the identified biclusters, especially, we use a blank line to separate the positively and the negatively (if any) correlated genes in each bicluster. As to '`.chars`' file', it provides the qualitative matrix of the microarray data for users with some details of how to discrete the data in '`.rules`' file'. You can find further details about how to represent a microarray data set with a qualitative matrix in our paper.

# Parameters of QUBIC #
Parameters of QUBIC
QUBIC has a number of parameters, namely,
  * the range `r` of possible ranks
  * the percentage `q` of the regulating conditions for each gene
  * the required consistency level `c` for a bicluster (default value is 0.95).
  * the desired number `o` of the output biclusters (default value is 100).
  * and the control parameter `f` for overlaps among to-be-identified biclusters.

For each of these parameters, we allow the user to adjust the default value to provide some flexibility.
  * The parameters `r` and `q` affect the granularity of the biclusters. A user is recommended to start with a small value of r (the default value is 1 so the corresponding data matrix consists of values '+1', '–1' and '0'), evaluate the results, and then use larger values (should not be larger than half of the number of the columns) to look for fine structures within the identified biclusters. The choice of `q`'s value depends on the specific application goals; that is if the goal is to find genes that are responsive to local regulators, we should use a relatively small `q`-value; otherwise we may want to consider larger `q`-values. The default value of `q` is 0.06 in QUBIC (this value is selected based on the optimal biclustering results on simulated data).
  * We have a parameter `f` to control the level of overlaps between to-be-identified biclusters (not discussed in the above algorithm); its default value is set to 1 to ensure that no two reported biclusters overlap more than `f`.
  * QUBIC also provides the option (`-d`) that a user can skip the step of using ranks to represent the actual gene expression values to go directly to the biclustering step on the provided matrix.

# Advanced Usage #

A new function that can expand identified biclusters in specific environment. Suppose you have two expression matrices **A** and **B**, where **B** is subset of **A**, you can extend the biclusters of **B** in the matrix **A** as following:

```
$./qubic -i A
$./qubic -i B
$./qubic -i A.chars -b B.blocks -s
```

The program will generate a '`B.blocks.expansion`' file, containing the enlarged biclusters in 'B.blocks'.

# Results access #
The related data sets and results can be found at http://csbl.bmb.uga.edu/~maqin/bicluster/benchmark.html.
