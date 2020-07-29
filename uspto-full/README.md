# USPTO-full

The directory contains processed USPTO-full datasets 


uspto_trainT.csv
uspto_valT.csv
uspto_testT.csv


from Dai, H.; Li, C.; Coley, C.; Dai, B.; Song, L. In Retrosynthesis Prediction with Conditional Graph Logic Network, Advances in Neural Information Processing Systems 32, Wallach, H.; Larochelle, H.; Beygelzimer, A.; Alché-Buc, F.; Fox, E.; Garnett, R., Eds. Curran Associates, Inc: 2019; pp 8872-8882 

The files were processed to remove reactions which contained no products or just single ions as reactants (e.g., US08163899B2,>>[OH2:11], US06048982,CC(=O)OCCCCC[I:22]>>[I-:22], US07425593B2,>>[K:12], US08114877B2,CC[I:13]>>[I-:13]). We also eliminated reactions where reactants had less than five atoms in total, since these were likely to be incorrect reactions.

All files are compressed using xz (to use it on Mac install xz in MacPort; xz is available on Ubuntu by default )

Results obtained with these files are described in a pre-print, which is submitted for a review.

Tetko, I. V.; Karpov, P.; Van Deursen, R.; Godin, G. Augmented Transformer for Direct and Single-Step Retrosynthesis Predictions Overperforms All Published Approaches. Nat. Commun. 2020, submitted.
