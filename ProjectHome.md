# **QUBIC: QUalitative BIClustering algorithm** #
## Algorithm ##

A qualitative biclustering algorithm for analyses of gene expression data.
Biclustering extends the traditional clustering techniques by attempting to find (all) subgroups of genes with similar expression patterns under to-be-identified subsets of experimental conditions when applied to gene expression data. Still the real power of this clustering strategy is yet to be fully realized due to the lack of effective and efficient algorithms for reliably solving the general biclustering problem. We report a QUalitative BIClustering algorithm (QUBIC) that can solve the biclustering problem in a more general form, compared to existing algorithms, through employing a combination of qualitative (or semi-quantitative) measures of gene expression data and a combinatorial optimization technique. One key unique feature of the QUBIC algorithm is that it can identify all statistically significant biclusters including biclusters with the so-called ‘scaling patterns’, a problem considered to be rather challenging; another key unique feature is that the algorithm solves such general biclustering problems very efficiently, capable of solving biclustering problems with tens of thousands of genes under up to thousands of conditions in a few minutes of the CPU time on a desktop computer. We have demonstrated a considerably improved biclustering performance by our algorithm compared to the existing algorithms on various benchmark sets and data sets of our own.

**Citing us**: Li G, Ma Q, Tang H, Paterson AH, Xu Y. QUBIC: a qualitative biclustering algorithm for analyses of gene expression data. Nucleic Acids Res. 2009;37(15):e101. doi: 10.1093/nar/gkp491.

## Data Source ##
The related data sets and results can be found at http://csbl.bmb.uga.edu/~maqin/bicluster/benchmark.html.


## Web-server ##
**QServer**: http://csbl.bmb.uga.edu/publications/materials/ffzhou/QServer/

**Citing us**: Fengfeng Zhou, Qin Ma (co-first), Guojun Li, Ying Xu, QServer: a biclustering server for prediction and validation of co-expressed gene clusters, PLoS ONE 7(3): e32660. doi:10.1371/journal.pone.0032660, 2012


## Evaluation ##
Based on an evaluation paper (see reference below) published at Brief Bioinformatics by _Kemal Eren_ et al., QUBIC handle noise **much better** than the other biclustering programs and has the **highest** enriched bicluster ratio in real data sets.

Eren K, Deveci M, Küçüktunç O, Catalyürek UV, A comparative analysis of biclustering algorithms for gene expression data. **Brief Bioinform.**, 2013.


## Contact us ##
Any questions, problems, bugs are welcome and should be dumped to Qin Ma <**_maqin@uga.edu_**>