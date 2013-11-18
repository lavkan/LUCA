LUCA
====

R code to find the probability of presence of a gene at the root of an evolutionary tree using its presence absence data from the current day species

1. Enter the phyletic vector of the gene: The phyletic vector is a vector whose entries are 0, 1 depending on the prsence and absence of the gene. In some models, it could be m, denoting multiple copies of the same gene is present.

2. Enter the tree

Output: Probability of 0, 1 or m of the gene at the root and all internal vertices of the tree. The probability is obtained using a maximum likelihood approach. Also as outputs, one can see the optimal value of likelihood, with the optimal parameters of the corresponding model.

