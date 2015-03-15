A qualitative biclustering algorithm for analyses of gene expression data


# Usage #
This software provides a biclustering program for microarray data. For a set of genes and a set of conditions, the program outputs a block-like structure which shows uniform pattern within the block, the block would contain only subsets of all given genes under subsets of all given conditions.

Certain parts of the code uses open-source data structure library codes, including:
  * fib <http://resnet.uoregon.edu/~gurney_j/jmpc/fib.html>, copyright information in fib.c
  * Mark A. Weiss's data structure codes <http://www.cs.fiu.edu/~weiss/>

# Installation #
Simply put "qubic.tar.gz" in any directory,type
> _$ tar zxvf qubic.tar.gz_
enter the folder '**qubic**' and type '**make**' then the compiled codes are within the same directory as the source.

# Inputs and outputs #
The major program in the provided package is '**qubic**', it can parse two formats of files, discrete data and continuous data, and examples for each are provided in the same folder.

  * To See help and look at all available options.

> _$ ./qubic -h_

  * Take a look at `toy_example` (discrete data) first. And try to run biclustering

> _$ ./qubic -i toy\_example -d_

> where -d is important here since it tells the program that this is discrete data.

  * Then look at a continuous data, example. Try to run

> _$ ./qubic -i example -f .25_

> This restricts no two output blocks overlap more than 0.25 of the size of each one. And the other parameters are default value.

  * For each input file, our program generates three output files, namely,'**.blocks' file**', '**.chars'file**' and '**.rules' file**'.

> In '**.blocks' file**', you can see all the identified biclusters, especially, we use a blank line to separate the positively and the negatively (if any) correlated genes in each bicluster. As to '**.chars' file**', it provides the qualitative matrix of the microarray data for users with some details of how to discrete the data in '**.rules' file**'. You can find further details about how to represent a microarray data set with a qualitative matrix in our paper.

# Parameters of QUBIC #
Parameters of QUBIC
QUBIC has a number of parameters, namely,
  * the range **_r_** of possible ranks
  * the percentage **_q_** of the regulating conditions for each gene
  * the required consistency level **_c_** for a bicluster (default value is 0.95).
  * the desired number **_o_** of the output biclusters (default value is 100).
  * and the control parameter **_f_** for overlaps among to-be-identified biclusters.

For each of these parameters, we allow the user to adjust the default value to provide some flexibility.
  * The parameters **_r_** and **_q_** affect the granularity of the biclusters. A user is recommended to start with a small value of r (the default value is 1 so the corresponding data matrix consists of values ‘+1’, ‘–1’ and ‘0’), evaluate the results, and then use larger values (should not be larger than half of the number of the columns) to look for fine structures within the identified biclusters. The choice of **_q_**'s value depends on the specific application goals; that is if the goal is to find genes that are responsive to local regulators, we should use a relatively small **_q_**-value; otherwise we may want to consider larger **_q_**-values. The default value of **_q_** is 0.06 in QUBIC (this value is selected based on the optimal biclustering results on simulated data).
  * We have a parameter **_f_** to control the level of overlaps between to-be-identified biclusters (not discussed in the above algorithm); its default value is set to 1 to ensure that no two reported biclusters overlap more than **_f_**.
  * QUBIC also provides the option (**_-d_**) that a user can skip the step of using ranks to represent the actual gene expression values to go directly to the biclustering step on the provided matrix.

# Advanced Usage #

  * A new function that can expand identified biclusters in specific environment. Suppose you have two expression matrices **A** and **B**, where **B** is subset of **A**, you can extend the biclusters of **B** in the matrix **A** as following,
    1. $./qubic -i A
    1. $./qubic -i B
    1. $./qubic -i A.chars -b B.blocks -s
> and the program will generate a '**B.blocks.expansion**' file, containing the enlarged biclusters in 'B.blocks'.

# Results access #
The related data sets and results can be found at http://csbl.bmb.uga.edu/~maqin/bicluster/benchmark.html.

# Contact #
Any questions, problems, bugs are welcome and should be dumped to
  * Qin Ma <maqin@uga.edu>
  * Haibao Tang <htang@jcvi.org>

Creation: Dec. 16, 2008