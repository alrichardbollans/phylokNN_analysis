import os

import pandas as pd

from evaluation.helper_functions import binary_model_order, bin_model_names


def main():
    # find ev scenario where second best ML model overall does best against traditional models
    closest_score = 10
    closest_file = ''
    closest_model = ''
    for file in os.listdir(os.path.join('ev_model_outputs', 'binary')):
        if file.endswith('.csv') and 'Clonality' not in file:
            scores = pd.read_csv(os.path.join('ev_model_outputs', 'binary', file), index_col=0)
            corhmm_score = scores['corHMM'].loc['mean']
            for model in ['logit_PEM']:
                if model != 'corHMM':
                    model_score = scores[model].loc['mean']
                    diff = model_score-corhmm_score
                    ratio = diff / corhmm_score
                    if ratio<closest_score:
                        closest_score = ratio
                        closest_file = file
                        closest_model = model


    #BiSSE
    print(closest_file)
    print(closest_score)
    print(closest_model)

    closest_score = 10
    closest_file = ''
    closest_model = ''
    for file in os.listdir(os.path.join('ev_model_outputs', 'continuous')):
        if file.endswith('.csv') and 'Seed' not in file:
            scores = pd.read_csv(os.path.join('ev_model_outputs', 'continuous', file), index_col=0)
            phylopars_score = scores['phylopars'].loc['mean']
            for model in ['linear_PEM']:
                if model != 'phylopars':
                    model_score = scores[model].loc['mean']
                    diff = model_score - phylopars_score
                    ratio = diff / phylopars_score
                    if ratio < closest_score:
                        closest_score = ratio
                        closest_file = file
                        closest_model = model

    # BMT
    print(closest_file)
    print(closest_score)
    print(closest_model)


if __name__ == '__main__':
    main()
