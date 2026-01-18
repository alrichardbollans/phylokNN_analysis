import sys

from sklearn.linear_model import LogisticRegression, LinearRegression
from tqdm import tqdm
from xgboost import XGBClassifier, XGBRegressor

from imputation.run_encodings_predictions import get_umap_data, add_y_to_data, get_eigenvectors, get_PEM_eigenvectors, get_autoencoded_data, \
    logit_init_kwargs, logit_grid_search_params, fit_and_output, xgb_clf_init_kwargs, xgb_clf_grid_search_params, get_semi_supervised_umap_data, \
    get_semi_supervised_autoencoded_data, linear_init_kwargs, linear_grid_search_params, xgb_rgr_init_kwargs, xgb_rgr_grid_search_params

sys.path.append('../..')
from data.helper_functions import number_of_simulation_iterations
from imputation.helper_functions import phylnn_predict, get_prediction_data_paths, get_bin_or_cont_from_ev_model

large_ev_models = ['large_BiSSE', 'large_BMT']

def main():
    for iteration in tqdm(range(1, 10 + 1)):
        for m in ['mcar']:
            for case in ['ultrametric']:
                for ev_model in large_ev_models:
                    phylnn_predict(case, ev_model, iteration, m)
    for iteration in tqdm(range(1, 10 + 1)):
        for m in ['mcar']:
            for case in ['ultrametric']:
                for ev_model in large_ev_models:
                    umap_X = get_umap_data(case, ev_model, iteration)
                    umap_df, umap_encoding_vars, umap_target_name = add_y_to_data(umap_X, case, ev_model, iteration, m)

                    eigen_X = get_eigenvectors(case, ev_model, iteration)
                    eigen_df, eigen_encoding_vars, eigen_target_name = add_y_to_data(eigen_X, case, ev_model, iteration, m)

                    # PEM_X = get_PEM_eigenvectors(case, ev_model, iteration, m)
                    # PEM_df, PEM_encoding_vars, PEM_target_name = add_y_to_data(PEM_X, case, ev_model, iteration, m)

                    autoenc_X = get_autoencoded_data(case, ev_model, iteration)
                    autoenc_df, autoenc_encoding_vars, autoenc_target_name = add_y_to_data(autoenc_X, case, ev_model, iteration, m)
                    out_dir = get_prediction_data_paths(case, ev_model, iteration, m)

                    bin_or_cont = get_bin_or_cont_from_ev_model(ev_model)
                    if bin_or_cont == 'binary':
                        # Compare logistic regression and XGBoost models i.e. for modelling simpler relationships and complex relationships
                        clf_instance = LogisticRegression(**logit_init_kwargs)
                        fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_umap', umap_df, umap_encoding_vars, umap_target_name,
                                       bin_or_cont)

                        clf_instance = LogisticRegression(**logit_init_kwargs)
                        fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_eigenvecs', eigen_df, eigen_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        # clf_instance = LogisticRegression(**logit_init_kwargs)
                        # fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_PEM', PEM_df, PEM_encoding_vars, PEM_target_name,
                        #                bin_or_cont)
                        #
                        clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_umap', umap_df, umap_encoding_vars, umap_target_name,
                                       bin_or_cont)

                        clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_eigenvecs', eigen_df, eigen_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        # clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        # fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_PEM', PEM_df, PEM_encoding_vars, PEM_target_name,
                        #                bin_or_cont)
                        #
                        # # # ### Semisupervised umap
                        semi_supervised_umap_df, semi_sup_umap_encoding_vars, semi_sup_umap_target_name = get_semi_supervised_umap_data(case,
                                                                                                                                        ev_model,
                                                                                                                                        iteration, m)
                        clf_instance = LogisticRegression(**logit_init_kwargs)
                        fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_umap_supervised', semi_supervised_umap_df,
                                       semi_sup_umap_encoding_vars,
                                       semi_sup_umap_target_name, bin_or_cont)
                        clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_umap_supervised', semi_supervised_umap_df,
                                       semi_sup_umap_encoding_vars,
                                       semi_sup_umap_target_name, bin_or_cont)

                        # ### autoencoder
                        clf_instance = LogisticRegression(**logit_init_kwargs)
                        fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_autoencoded', autoenc_df, autoenc_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_autoencoded', autoenc_df, autoenc_encoding_vars,
                                       umap_target_name,
                                       bin_or_cont)

                        ## Semisupervised autoenc
                        semi_supervised_autoenc_df, semi_sup_autoenc_encoding_vars, semi_sup_autoenc_target_name = get_semi_supervised_autoencoded_data(
                            case, ev_model,
                            iteration, m)

                        clf_instance = LogisticRegression(**logit_init_kwargs)
                        fit_and_output(clf_instance, logit_grid_search_params, out_dir, 'logit_autoenc_supervised', semi_supervised_autoenc_df,
                                       semi_sup_autoenc_encoding_vars,
                                       semi_sup_autoenc_target_name, bin_or_cont)
                        clf_instance = XGBClassifier(**xgb_clf_init_kwargs)
                        fit_and_output(clf_instance, xgb_clf_grid_search_params, out_dir, 'xgb_autoenc_supervised', semi_supervised_autoenc_df,
                                       semi_sup_autoenc_encoding_vars,
                                       semi_sup_autoenc_target_name, bin_or_cont)

                    elif bin_or_cont == 'continuous':
                        clf_instance = LinearRegression(**linear_init_kwargs)
                        fit_and_output(clf_instance, linear_grid_search_params, out_dir, 'linear_umap', umap_df, umap_encoding_vars, umap_target_name,
                                       bin_or_cont)

                        clf_instance = LinearRegression(**linear_init_kwargs)
                        fit_and_output(clf_instance, linear_grid_search_params, out_dir, 'linear_eigenvecs', eigen_df, eigen_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        # clf_instance = LinearRegression(**linear_init_kwargs)
                        # fit_and_output(clf_instance, linear_grid_search_params, out_dir, 'linear_PEM', PEM_df, PEM_encoding_vars, PEM_target_name,
                        #                bin_or_cont)

                        clf_instance = XGBRegressor(**xgb_rgr_init_kwargs)
                        fit_and_output(clf_instance, xgb_rgr_grid_search_params, out_dir, 'xgb_umap', umap_df, umap_encoding_vars, umap_target_name,
                                       bin_or_cont)

                        clf_instance = XGBRegressor(**xgb_rgr_init_kwargs)
                        fit_and_output(clf_instance, xgb_rgr_grid_search_params, out_dir, 'xgb_eigenvecs', eigen_df, eigen_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        # clf_instance = XGBRegressor(**xgb_rgr_init_kwargs)
                        # fit_and_output(clf_instance, xgb_rgr_grid_search_params, out_dir, 'xgb_PEM', PEM_df, PEM_encoding_vars, PEM_target_name,
                        #                bin_or_cont)

                        ### autoencoder
                        clf_instance = LinearRegression(**linear_init_kwargs)
                        fit_and_output(clf_instance, linear_grid_search_params, out_dir, 'linear_autoencoded', autoenc_df, autoenc_encoding_vars,
                                       eigen_target_name, bin_or_cont)

                        clf_instance = XGBRegressor(**xgb_rgr_init_kwargs)
                        fit_and_output(clf_instance, xgb_rgr_grid_search_params, out_dir, 'xgb_autoencoded', autoenc_df, autoenc_encoding_vars,
                                       umap_target_name,
                                       bin_or_cont)


if __name__ == '__main__':
    main()
