#  Phylogenetic Nearest Neighbours analysis

This repository includes analysis of the phyloKNN model (https://github.com/alrichardbollans/phylokNN) against other conventional and less conventional approaches.

## Citation


## A Note on Tree Size
The analysis in the related publication is based on trees with 100 extant species. As noted in the manuscript, it is of course possible that the machine learning methods perform better using trees with more species records for training. We have tested some cases using trees with 1000 extant species and have found similar results [large_compiled_score_outputs](https://github.com/alrichardbollans/phylokNN_analysis/tree/master/evaluation/large_compiled_score_outputs). To limit excessive computation, we looked at only ultrametric trees in the MCAR scenario. We looked at only a single evolutionary simulation model for both binary and continuous cases, which was selected by looking in the main analysis where the performance of one of the ML methods was closest to the corHMM and phylopars models, which turned out to be BiSSE and BMT. This is a limited result which would require further investigation, especially with real traits, but helps to support the original hypothesis.

## Licence

This work is licensed under a
[Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License][cc-by-nc-sa].

[![CC BY-NC-SA 4.0][cc-by-nc-sa-image]][cc-by-nc-sa]

[cc-by-nc-sa]: http://creativecommons.org/licenses/by-nc-sa/4.0/

[cc-by-nc-sa-image]: https://licensebuttons.net/l/by-nc-sa/4.0/88x31.png

[cc-by-nc-sa-shield]: https://img.shields.io/badge/License-CC%20BY--NC--SA%204.0-lightgrey.svg
