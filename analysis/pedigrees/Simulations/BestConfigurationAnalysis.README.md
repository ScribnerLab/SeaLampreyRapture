<img src="./Figure.PNG" width="1430" style="display: block; margin: auto;" />

There are there major steps when analyzing a COLONY BestConfig (BC) file
generated using simulated offspring genotypes where the user provides no
parental genotypes. When COLONY does not have candidate parental
genotypes it will infer a parent, and its associated genotypes. These
parents are given generic names (e.g., *1 or \#3). Given that the user
does not know which column represents each sex, each parent (P) column
is given a generic column name (P1, P2), and the inferred parent IDs are
changed slightly to reflect the header columns (e.g., *1 and \#3 are
changed to P1\_1 and P2\_3). It is also important to note that the
offspring ID column (off) contains information about the known parents
that produced each offspring, which is embedded in each offspring ID.
For instance, the user knows that Mom 1 (m1) and Dad 1 (d1) produced
offspring “o1\_m1\_d1”. Step 2 consists of evaluating all unique
pairwise comparisons among all offspring genotyped (i.e., dyads). For
each dyad (i.e. row) the user can determine the inferred relationship
(IR) and known relationship (KR) between the two individuals by counting
how many parents the two individuals share. That is, full-siblings (FS,
2 parents shared), half-siblings (HS, 1 parent shared), and unrelated
(UR, 0 parents shared) dyads can be determined using both the Inferred
and Known Parent IDs. Using information in step three the accuracy and
false negative rate can be calculated for each simulation. See
definitions and example in Figure for more details.
