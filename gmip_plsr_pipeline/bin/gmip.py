########################################################################################################################

import argparse
import os
import sys
import scipy
import random
import pickle
import pandas as pd
import numpy as np
from sklearn import preprocessing
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
from sklearn.decomposition import IncrementalPCA
from sklearn.decomposition import PCA
from sklearn.metrics import silhouette_score
import logging

########################################################################################################################

def get_args(args=None):
  parser = argparse.ArgumentParser(description='Run GMIP: GWAS & Multiomics Integration Pipeline in python')
  parser.add_argument("--trainFile",
    help="File with predictor matrix X and Y after processing from pops. (.traindata after covariates projected out.)")
  parser.add_argument("--matFile",
    help="File with predictor matrix X after processing from pops. (.traindata before covariates projected out.)")
  parser.add_argument("--marginalFile",
    help="File with marginal OLS results after processing from pops. (.marginals after covariates projected out.)")
  parser.add_argument("--xIDColNum", default=0,
    help="default:0 (based on pops traindata output)")
  parser.add_argument("--yColName", default="Y_train",
    help="default:'Y_train' (based on magma output)")
  parser.add_argument("--noFeatureClustering", dest='doFeatureClustering', action='store_false',
    help="Do not perform feature clustering before regression.")
  parser.set_defaults(doFeatureClustering=True)
  parser.add_argument("--outPrefix",
    help="Prefix to be added to the file name. Cannot be NA.")
  parser.add_argument("--methodName", default="RidgeCV",
    help="Name of the method to use. Using the exact name is important. Available methods: RidgeCV, PC_RidgeCV")
  parser.add_argument("--randomSeed", type=int, default=42,
    help="RandomSeed for reproducibility. default:42")
  return parser.parse_args(args)

########################################################################################################################

def cluster_features2(x, explained_variance=0.8, batch_size=None):
  x_t=x.transpose()
  ipca, n_components, x_reduced_df=reduction_using_PCA(x_t, explained_variance=explained_variance, batch_size=batch_size)
  corr_matrix=x_reduced_df.transpose().corr().values
  corr_matrix=np.clip(corr_matrix, -1, 1)
  squared_corr_matrix=corr_matrix**2
  squared_corr_matrix[squared_corr_matrix < 1e-5] = 0
  distance_matrix=1-squared_corr_matrix
  condensed_distance_matrix=squareform(distance_matrix, checks=False)
  linked=linkage(condensed_distance_matrix, method='complete')
  threshold_range = np.linspace(0.01, 1, 10)
  max_silhouette_score = -1
  best_threshold = 0.2
  square_distance_matrix = squareform(condensed_distance_matrix)
  max_clusters=x_reduced_df.shape[0]-1
  for threshold in threshold_range:
    clusters = fcluster(linked, t=threshold, criterion='distance')
    num_clusters = len(np.unique(clusters))
    if num_clusters >=max_clusters:
      continue
    if len(np.unique(clusters)) > 1:  # At least 2 clusters are needed to compute silhouette score
      score = silhouette_score(square_distance_matrix, clusters, metric='precomputed')
      if score > max_silhouette_score:
        max_silhouette_score = score
        best_threshold = threshold
  clusters=fcluster(linked, t=best_threshold, criterion='distance')
  print(clusters.max())
  feature_names=x.columns
  feature_to_cluster=dict(zip(feature_names, clusters))
  cluster_row=pd.Series(feature_to_cluster, name='Cluster')
  cluster_row_reindexed=cluster_row.reindex(x_t.index)
  cluster_row_reindexed=cluster_row_reindexed.reset_index()
  return cluster_row_reindexed, best_threshold, max_silhouette_score

########################################################################################################################

def reduction_using_PCA(x, explained_variance=0.8, batch_size=None):
  ipca=IncrementalPCA(batch_size=batch_size)
  x_reduced=ipca.fit_transform(x)
  # Selecting the number of components that explain the predefined variance
  cumulative_variance=np.cumsum(ipca.explained_variance_ratio_)
  n_components=np.where(cumulative_variance >= explained_variance)[0][0] + 1
  x_reduced=x_reduced[:, :n_components]
  component_names=[f'PC{i+1}' for i in range(n_components)]
  x_reduced_df=pd.DataFrame(x_reduced, columns=component_names, index=x.index)
  x_reduced_df=x_reduced_df.astype(np.float32)
  return ipca, n_components, x_reduced_df

########################################################################################################################

def gmip_predict(mat, rows, cols, coefs_df):
    pred = mat.dot(coefs_df.loc[cols].beta.values)
    preds_df = pd.DataFrame([rows, pred]).T
    preds_df.columns = ["ENSGID", "GMIP_Score"]
    return preds_df

########################################################################################################################
def main(config):
  logging.basicConfig(format="%(levelname)s: %(message)s", level=logging.DEBUG)
  logging.info("Seeting random seed.")
  np.random.seed(config["randomSeed"])
  random.seed(config["randomSeed"])
  logging.info("Reading trainFile.")
  xy_df=pd.read_csv(config["trainFile"], sep="\t", index_col=config['xIDColNum'])
  mat_df=pd.read_csv(config["matFile"], sep="\t", index_col=config['xIDColNum'])
  marginal_df=pd.read_csv(config["marginalFile"], sep="\t", index_col=config['xIDColNum'])
  X_train=xy_df.drop(config['yColName'], axis=1)
  Y_train=xy_df[config['yColName']]
  mat_df=mat_df[X_train.columns]
  if config['doFeatureClustering'] == True:
    logging.info("Feature Clustering will be performed and used.")
    cluster_row_reindexed, best_threshold, max_silhouette_score=cluster_features2(X_train)
    feature_df=pd.merge(left=marginal_df, right=cluster_row_reindexed, how="outer", left_index=True,
        right_on='index')
    feature_df_clean=feature_df.dropna(subset=['Cluster'])
    feature_df_clean=feature_df_clean.copy()
    feature_df_clean.loc[:, 'abs_beta']=feature_df_clean['beta'].abs()
    final_features_df=feature_df_clean.sort_values('abs_beta', ascending=False).groupby('Cluster').head(1)
    features_selected_list=final_features_df['index'].to_list()
    X_train=X_train[features_selected_list]
    mat_df=mat_df[features_selected_list]
    mat_df=mat_df[X_train.columns]
  if config['methodName'] == 'RidgeCV':
    from sklearn.linear_model import RidgeCV
    alphas=np.logspace(-2, 10, num=25)
    reg=RidgeCV(fit_intercept=False, alphas=alphas)
    cols=X_train.columns.to_list()
    rows=mat_df.index.to_list()
    mat_df=mat_df[X_train.columns]
    logging.info("Model = RidgeCV with 25 alphas, generalized leave-one-out cross-validation, NMSE as scoring metric.")
    reg.fit(X_train, Y_train)
    coefs_df=pd.DataFrame([["METHOD", "RidgeCV"], ["SELECTED_CV_ALPHA", reg.alpha_], ["BEST_CV_SCORE", reg.best_score_]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, reg.coef_]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)
  elif config['methodName'] == 'PC_RidgeCV':
    from sklearn.linear_model import RidgeCV
    alphas=np.logspace(-2, 10, num=25)
    reg=RidgeCV(fit_intercept=False, alphas=alphas)
    ipca, n_components, X_train_reduced_df=reduction_using_PCA(X_train)
    mat_reduced=ipca.transform(mat_df)
    mat_reduced=mat_reduced[:, :n_components]
    component_names=[f'PC{i+1}' for i in range(n_components)]
    mat_reduced_df=pd.DataFrame(mat_reduced, columns=component_names, index=mat_df.index)
    mat_reduced_df=mat_reduced_df.astype(np.float32)
    cols=mat_reduced_df.columns.to_list()
    rows=mat_reduced_df.index.to_list()
    logging.info("Model = PC_RidgeCV with 25 alphas, generalized leave-one-out cross-validation, NMSE as scoring metric.")
    reg.fit(X_train_reduced_df, Y_train)
    coefs_df=pd.DataFrame([["METHOD", "PC_RidgeCV"], ["SELECTED_CV_ALPHA", reg.alpha_], ["BEST_CV_SCORE", reg.best_score_]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, reg.coef_]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_reduced_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)
  elif config['methodName'] == 'OrthogonalMatchingPursuitCV':
    from sklearn.linear_model import OrthogonalMatchingPursuitCV
    reg=OrthogonalMatchingPursuitCV(fit_intercept=False)
    cols=X_train.columns.to_list()
    rows=mat_df.index.to_list()
    mat_df=mat_df[X_train.columns]
    logging.info("Model = OrthogonalMatchingPursuitCV.")
    reg.fit(X_train, Y_train)
    coefs_df=pd.DataFrame([["METHOD", "OrthogonalMatchingPursuitCV"]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, reg.coef_]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)
  elif config['methodName'] == 'PLSRegression_nc1':
    from sklearn.cross_decomposition import PLSRegression
    reg=PLSRegression(n_components=1)
    cols=X_train.columns.to_list()
    rows=mat_df.index.to_list()
    mat_df=mat_df[X_train.columns]
    logging.info("Model = PLSRegression_nc1.")
    reg.fit(X_train, Y_train)
    coefs_df=pd.DataFrame([["METHOD", "PLSRegression"], ["SELECTED_n_components", 1]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, reg.coef_[0]]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)
  elif config['methodName'] == 'PLSRegression_nc2':
    from sklearn.cross_decomposition import PLSRegression
    reg=PLSRegression(n_components=2)
    cols=X_train.columns.to_list()
    rows=mat_df.index.to_list()
    mat_df=mat_df[X_train.columns]
    logging.info("Model = PLSRegression_nc2.")
    reg.fit(X_train, Y_train)
    coefs_df=pd.DataFrame([["METHOD", "PLSRegression"], ["SELECTED_n_components", 2]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, reg.coef_[0]]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)
  elif config['methodName'] == 'PLSRegression_nc3':
    from sklearn.cross_decomposition import PLSRegression
    reg=PLSRegression(n_components=3)
    cols=X_train.columns.to_list()
    rows=mat_df.index.to_list()
    mat_df=mat_df[X_train.columns]
    logging.info("Model = PLSRegression_nc3.")
    reg.fit(X_train, Y_train)
    coefs_df=pd.DataFrame([["METHOD", "PLSRegression"], ["SELECTED_n_components", 3]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, reg.coef_[0]]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)
  elif config['methodName'] == 'PLSRegression_nc5':
    from sklearn.cross_decomposition import PLSRegression
    reg=PLSRegression(n_components=5)
    cols=X_train.columns.to_list()
    rows=mat_df.index.to_list()
    mat_df=mat_df[X_train.columns]
    logging.info("Model = PLSRegression_nc5.")
    reg.fit(X_train, Y_train)
    coefs_df=pd.DataFrame([["METHOD", "PLSRegression"], ["SELECTED_n_components", 5]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, reg.coef_[0]]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)
  elif config['methodName'] == 'PLSRegression_nc10':
    from sklearn.cross_decomposition import PLSRegression
    reg=PLSRegression(n_components=10)
    cols=X_train.columns.to_list()
    rows=mat_df.index.to_list()
    mat_df=mat_df[X_train.columns]
    logging.info("Model = PLSRegression_nc10.")
    reg.fit(X_train, Y_train)
    coefs_df=pd.DataFrame([["METHOD", "PLSRegression"], ["SELECTED_n_components", 10]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, reg.coef_[0]]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)
  elif config['methodName'] == 'PLSRegressionCV':
    from sklearn.model_selection import GridSearchCV
    from sklearn.cross_decomposition import PLSRegression
    from sklearn.pipeline import Pipeline
    cols=X_train.columns.to_list()
    rows=mat_df.index.to_list()
    mat_df=mat_df[X_train.columns]
    logging.info("Model = PLSRegressionCV.")
    max_components=min(15, X_train.shape[1])
    parameters_pls = {'pls__n_components': np.arange(1, max_components + 1)}
    plsr_pipe=Pipeline([('pls', PLSRegression())])
    grid_search=GridSearchCV(plsr_pipe, parameters_pls, cv=3, scoring='neg_mean_squared_error', verbose=3, n_jobs=-1)
    grid_search.fit(X_train, Y_train)
    # Best parameters and score
    print("Best parameters:", grid_search.best_params_)
    print("Best score:", grid_search.best_score_)
    best_model = grid_search.best_estimator_
    coefs_df=pd.DataFrame([["METHOD", "PLSRegression"], ["SELECTED_n_components", grid_search.best_params_['pls__n_components']]])
    coefs_df=pd.concat([coefs_df, pd.DataFrame([cols, best_model.named_steps['pls'].coef_[0]]).T])
    coefs_df.columns=["parameter", "beta"]
    coefs_df=coefs_df.set_index("parameter")
    coefs_df.to_csv(config['outPrefix'] + ".coefs", sep="\t")
    preds_df=gmip_predict(mat_df, rows, cols, coefs_df)
    preds_df.to_csv(config['outPrefix'] + ".preds", sep="\t", index=False)

########################################################################################################################

if __name__ == '__main__':
  config=vars(get_args())
  main(config)

# config={
#   'trainFile': 'test_dir/2_pops/pops.test_chr1.traindata', 'matFile': 'test_dir/2_pops/pops.test_chr1.matdata',
#   'marginalFile': 'test_dir/2_pops/pops.test_chr1.marginals',
#   'outPrefix': 'gmip.withoutFC.PLSRegressionCV',
#   'methodName': 'PLSRegressionCV',
#   'doFeatureClustering': False,
#   'xIDColNum': 0, 'yColName': 'Y_train','randomSeed': 42
# }
