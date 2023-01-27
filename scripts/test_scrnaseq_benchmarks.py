import os, sys

from reader import *
from clusterer import *
from validator import *
from visualizer import *
from preprocesser import *

def main(dirName):
   report = pd.DataFrame(columns = ['Benchmark', 'Instance #', 'Attribute #',\
         'Cluster #','ga-sw', 'k-sw', 'phg-sw', 'scp-sw'])

   numCenters_list = []
   phgMembership_list = []
   for dataset in os.listdir(dirName):
      dataFile = dirName + dataset
      data, cellId = read_scrnaseq_data(dataFile)
      phgMembership = run_phenograph(data)
      phgMembership_list.append(phgMembership)
      numCenters = len(np.unique(phgMembership))
      numCenters_list.append(numCenters)
   
   numCenters_list = scale(numCenters_list, 3, 10)
   numCenters_list = list(map(int, numCenters_list))
   
   index = 0
   for dataset in os.listdir(dirName):   
      dataFile = dirName + dataset
      print(dataFile)
      data, cellId = read_scrnaseq_data(dataFile)
      numCenters = numCenters_list[index]
      
      ga = run_mo_gakmeans(datasets = data, numCenters = numCenters, numGen = 350, sizePop = 600)
      ga.assign_memberships(data)
      
      kmeans = run_kmeans(data, numCenters,maxiter = 350)
      kmeans.assign_memberships(data)

      scanpy_memberships = run_scanpy_clustering(dataFile)
      
      silhouette_ga = silhouette(data, ga.memberships)
      silhouette_k  = silhouette(data, kmeans.memberships) 
      silhouette_phg = silhouette(data, phgMembership_list[index])
      silhouette_scp = silhouette(data, scanpy_memberships)
      
      report.loc[len(report)] = [dataset, data.shape[0],data.shape[1], numCenters,"{:.2f}".format(silhouette_ga), \
         "{:.2f}".format(silhouette_k), "{:.2f}".format(silhouette_phg), "{:.2f}".format(silhouette_scp)]

      index += 1
      break #only run for one dataset for demo
   fileName = '../output/test_scrna_benchmarks.csv'
   report.to_csv(fileName)
   
if __name__ == '__main__':
   dirName = '../data/scrna_benchmarks_umap/'
   main(dirName)
