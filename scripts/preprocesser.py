import umap.umap_ as umap
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn.preprocessing import StandardScaler

#select top n attributes with high variance  
def select_highly_variable(data, cutoff):
   variance = data.var()
   attributes = list(variance.sort_values(ascending=False).index)
   lowVarianceAttributes = attributes[cutoff:]
   data.drop(lowVarianceAttributes, axis = 1, inplace=True) 
   return data

def standardize(data):
   columnNames = data.columns.values.tolist()
   data = StandardScaler().fit_transform(data)
   data = pd.DataFrame(data, columns = columnNames)
   return data
   
def normalize(data, vectorNorm):
   columnNames = data.columns.values.tolist()
   data = preprocessing.normalize(data,norm = vectorNorm)
   data = pd.DataFrame(data, columns = columnNames)
   return data

def run_pca(data, numComponents):
   pca = PCA(n_components = numComponents)
   columnNames = []
   for i in range(0,numComponents):
      columnNames.append('pc' + str(i+1))
   data = pd.DataFrame(data = pca.fit_transform(data), columns = columnNames)
   return data
   
def run_umap(data):
   columnNames = ['umap1','umap2']
   data = umap.UMAP().fit_transform(data)
   data = pd.DataFrame(data, columns = columnNames)
   return data

def scale(data, lower, upper):
   data = np.array(data).reshape(-1,1)
   min_max_scaler = preprocessing.MinMaxScaler(feature_range=(lower, upper))
   data = min_max_scaler.fit_transform(data)
   data = data.reshape(len(data)).tolist()
   return data
