import os
import pathlib

import pandas as pd
from matplotlib import pyplot as plt
import seaborn as sns
from sklearn.metrics import brier_score_loss, mean_absolute_error, average_precision_score

from data.helper_functions import number_of_simulation_iterations, simulation_types
from evaluation.compile_prediction_scores import read_all_results, evaluate_model_outputs_for_specific_case, output_results_from_df
from evaluation.helper_functions import bin_model_names, cont_model_names, rename_models_and_ev_models, \
    binary_model_order, continuous_model_order
from imputation.helper_functions import get_input_data_paths, get_prediction_data_paths, missingness_types, get_bin_or_cont_from_ev_model
from imputation.run_large_predictions import large_ev_models


def collate_simulation_outputs(ev_model: str,
                               range_to_eval: int = number_of_simulation_iterations, scorer=None, out_dir=None, signal=True):
    full_df = pd.DataFrame()
    for tag in range(1, range_to_eval + 1):
        for case in ['ultrametric']:
            if ev_model in ['Seed Mass', 'Clonality', 'APM'] and case == 'with_extinct':
                continue

            for missing_type in ['mcar']:
                if ev_model == 'APM' and missing_type == 'phyloNa':
                    continue
                run_dict = evaluate_model_outputs_for_specific_case(case, ev_model, tag, missing_type, scorer=scorer, signal=signal)
                run_df = pd.DataFrame(run_dict, index=[tag])
                full_df = pd.concat([full_df, run_df])
    if out_dir is None:
        out_dir = os.path.join('compiled_score_outputs', ev_model)
    output_results_from_df(full_df, out_dir, ev_model, scorer_label=scorer.__name__ if scorer is not None else 'Loss')

    return full_df

def main():
    for ev_model in large_ev_models:
        collate_simulation_outputs(ev_model, signal=False, range_to_eval=10, out_dir = os.path.join('large_compiled_score_outputs', ev_model))



if __name__ == '__main__':
    main()
