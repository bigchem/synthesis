# USPTO-full

The directory contains processed USPTO-full datasets from

Dai, H.; Li, C.; Coley, C.; Dai, B.; Song, L. In Retrosynthesis Prediction with Conditional Graph Logic Network, Advances in Neural Information Processing Systems 32, Wallach, H.; Larochelle, H.; Beygelzimer, A.; Alché-Buc, F.; Fox, E.; Garnett, R., Eds. Curran Associates, Inc: 2019; pp 8872-8882 

The files were processed to remove reactions which contained no products or just single ions as reactants (e.g., US08163899B2,>>[OH2:11], US06048982,CC(=O)OCCCCC[I:22]>>[I-:22], US07425593B2,>>[K:12], US08114877B2,CC[I:13]>>[I-:13]).
We also eliminated reactions where reactants had less than five atoms in total, since these were likely to be incorrect reactions.

All files are compressed using xz (to use it on Mac install xz in MacPort; xz is available on Ubuntu by default )

uspto_testT.csv.xz - processed test set

uspto_testT100.csv.xz  - augmeted 100 times test set: it is too large to be uploaded and thus it is splitted to:
uspto_testT100.csv.xz_a 
uspto_testT100.csv.xz_b 

uspto_testT100.csv.can.xz - canonised version of the same file

result_uspto_testT100.csv.xz - results predicted by the model for the test file: again, the file is too big and thus it is splitted to:
result_uspto_testT100.csv.xz_a 
result_uspto_testT100.csv.xz_b
result_uspto_testT100.csv.xz_c

result_uspto_testT100.csv.can.xz - canonised version of the result file

uspto_trainT.csv.xz - processed training set
uspto_valT.csv.xz - processed validation set
uspto_train_valTx5Shuff.csv.zz - combined training + validation set, which was augmented 5 times and used for model development.

These sets and results are described in our pre-print, which is submitted for a review.
