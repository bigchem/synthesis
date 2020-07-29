# SYNTHESIS

The directory contains data and results from article

Tetko, I. V.; Karpov, P.; Van Deursen, R.; Godin, G. Augmented Transformer for Direct and Single-Step Retrosynthesis Predictions Overperforms All Published Approaches. Nat. Commun. 2020, submitted.

The data used to develop models and predict test sets are provided in the respective directories

USPTO-50k
Schneider, N.; Lowe, D. M.; Sayle, R. A.; Tarselli, M. A.; Landrum, G. A., Big Data from Pharmaceutical Patents: A Computational Analysis of Medicinal Chemists' Bread and Butter. J. Med. Chem. 2016, 59 (9), 4385-402.

USPTO-MIT
Jin, W.; Coley, C.; Barzilay, R.; Jaakkola, T. Predicting Organic Reaction Outcomes with Weisfeiler-Lehman Network. In Advances in Neural Information Processing Systems 30; Guyon, I., Luxburg, U. V., Bengio, S., Wallach, H., Fergus, R., Vishwanathan, S., Garnett, R., Eds.; Curran Associates, Inc., 2017; pp 2607–2616.

USPTO-full
Dai, H.; Li, C.; Coley, C.; Dai, B.; Song, L. In Retrosynthesis Prediction with Conditional Graph Logic Network, Advances in Neural Information Processing Systems 32, Wallach, H.; Larochelle, H.; Beygelzimer, A.; Alché-Buc, F.; Fox, E.; Garnett, R., Eds. Curran Associates, Inc: 2019; pp 8872-8882 

For the latter set we remove reactions which contained no products or just single ions as reactants (e.g., US08163899B2,>>[OH2:11], US06048982,CC(=O)OCCCCC[I:22]>>[I-:22], US07425593B2,>>[K:12], US08114877B2,CC[I:13]>>[I-:13]). We also eliminated reactions where reactants had less than five atoms in total, since these were likely to be incorrect reactions.

The modified training and test set are
uspto_trainT.csv
uspto_testT.csv

The augmened versions of these files were used for  training (5x times) and testing (100x times) of the models.

All files are compressed using xz (to use it on Mac install, e.g. xz in MacPort; xz is available on some linux systems, e.g. Ubuntu by default )

These sets and results are described in our pre-print, which is submitted for a review now.

The USPTO-Full, USPTO-MIT results files were too large (>100MB in archive) to be uploaded to GitHub and thus we uploaded a subset of the test file with augmentation 20x instead of 100x. The processed canonical results files with x100, which were much better compressed, are included. The full test files can be received on the request. 

The accuracy of model predictions using 20x and 100x augmentations was very similar and within 0.1% for Top-1 predictions. 

The models will be made after the acceptance of the article.

Script compare.pl can be used to calculate statsitical performance using canonical files, e.g.

perl compare.pl patents_test100.csv.can result_patents_test100.csv.can 1

        will calculate top-1 performance

perl compare.pl patents_test100.csv.can result_patents_test100.csv.can 2

        will calculates top-2 perfrormance and so on. There are other processing options, which are listed in the script output. 

