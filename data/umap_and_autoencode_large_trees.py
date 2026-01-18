import os

import pandas as pd
from tqdm import tqdm

from data.helper_functions import input_data_path
from data.umapping import reduction_factor, unsupervised_umap_wrapper
from phyloAutoEncoder import autoencode_pairwise_distances


def main():
    for tag in tqdm(range(1, 10 + 1)):
        dir_path = os.path.join(input_data_path, 'simulations', 'ultrametric','BiSSE_large_trees', str(tag))
        distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)
        unsup_model, unsupervised = autoencode_pairwise_distances(distances, reduction_factor, dir_path,plot=False)
        unsupervised.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny.csv'))

        dir_path = os.path.join(input_data_path, 'simulations', 'ultrametric','large_trees', str(tag))
        distances = pd.read_csv(os.path.join(dir_path, 'tree_distances.csv'), index_col=0)
        unsup_model, unsupervised = autoencode_pairwise_distances(distances, reduction_factor, dir_path,plot=False)
        unsupervised.to_csv(os.path.join(dir_path, 'unsupervised_autoencoded_phylogeny.csv'))

        dir_path = os.path.join(input_data_path, 'simulations', 'ultrametric', 'BiSSE_large_trees', str(tag))
        unsupervised_umap_wrapper(dir_path)

        dir_path = os.path.join(input_data_path, 'simulations', 'ultrametric', 'large_trees', str(tag))
        unsupervised_umap_wrapper(dir_path)

if __name__ == '__main__':
    main()
