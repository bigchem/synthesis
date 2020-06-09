# USPTO-full

The directory contains processed USPTO-full datasets from

Dai, H.; Li, C.; Coley, C.; Dai, B.; Song, L. In Retrosynthesis Prediction with Conditional Graph Logic Network, Advances in Neural Information Processing Systems 32, Wallach, H.; Larochelle, H.; Beygelzimer, A.; Alché-Buc, F.; Fox, E.; Garnett, R., Eds. Curran Associates, Inc: 2019; pp 8872-8882 

The files were processed to remove reactions which contained no products or just single ions as reactants (e.g., US08163899B2,>>[OH2:11], US06048982,CC(=O)OCCCCC[I:22]>>[I-:22], US07425593B2,>>[K:12], US08114877B2,CC[I:13]>>[I-:13]).
We also eliminated such reactions as well as those where reactants had less than five atoms in total, since these were unlikely to be correct reactions.

uspto_testT.csv - processed test set
uspto_testT100.csv.bz2_xaa - part a
uspto_testT100.csv.bz2_xab - part b of uspto_testT.csv augmeted 100 times  (join both files before decompressing)
uspto_testT100.csv.can.bz2 - canonised version of the same file

result_uspto_testT100.csv.bz2_xaa - part a
result_uspto_testT100.csv.bz2_xab - part b of result_uspto_testT100.csv (join both files before decompressing)
result_uspto_testT100.csv.can.bz2 - canonised version of the same file

uspto_trainT.csv - processed training set
uspto_valT.csv - processed validation set
uspto_train_valTx5Shuff.csv.bz2 - combined training + validation set augmented 5 times

These sets and results were calculated as described in our pre-print, which is submitted for the review.
