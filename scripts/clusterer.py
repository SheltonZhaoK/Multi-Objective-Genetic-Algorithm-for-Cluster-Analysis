import random, time, phenograph, yaml, multiprocessing, sys, os

import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from deap import base, creator, tools, algorithms
from sklearn.cluster import KMeans
from math import dist

sys.path.append(os.path.dirname(__file__))
from validator import *

class ClusteringObject:
   def __init__(self, centers,  intraclusterDistance, time, fitness):
      self.centers = centers
      self.memberships = []
      self.intraclusterDistance = intraclusterDistance
      self.time = time
      self.fitness = fitness

   def print_report(self):
      print('fitness: %s \ntime(second): %f ' % (str(self.intraclusterDistance), self.time))
      print("centers:", self.centers, '\n' )

   def assign_memberships(self, data):
      data = data.values.tolist()
      for i in range(0, len(data)):
         distanceList = []
         for j in range(0,len(self.centers)):
            distanceList.append(dist(data[i],self.centers[j]))
         minDistance = min(distanceList)
         self.memberships.append(distanceList.index(minDistance))

with open('/deac/csc/khuriGrp/zhaok220/drug/scripts/MOGA/scripts/ga_parameters.yaml') as file:
   params = yaml.load(file, Loader=yaml.FullLoader)['MOEA']

def reshape_chromosome(chromosome, numCenter):
   listOfCenters = np.array(chromosome).reshape(numCenter,int(len(chromosome)/numCenter)).tolist()
   return listOfCenters

def compute_total_intracluster_distance(data, listOfCenters):
   intraclusterDistance = 0
   for i in range(0, len(data)):
      distanceList = []
      for j in range(0,len(listOfCenters)):
         distanceList.append(dist(data[i],listOfCenters[j]))
      minDistance = min(distanceList)
      intraclusterDistance += minDistance
   return intraclusterDistance

def compactness(pop):
   listOfCenters = reshape_chromosome(pop, numCenter)
   intraclusterDistance = compute_total_intracluster_distance(data, listOfCenters)
   return intraclusterDistance,

def separation(pop):
   intraclusterDistance = 0
   interclusterDistance = 0
   listOfCenters = reshape_chromosome(pop, numCenter)
   k = len(listOfCenters)
   counts = [0] * k
   for i in range(0, len(data)):
         distanceList = []
         for j in range(0,len(listOfCenters)):
            distanceList.append(dist(data[i],listOfCenters[j]))
         minDistance = min(distanceList)
         intraclusterDistance += minDistance
         center = distanceList.index(minDistance)
         counts[center] += 1

   for i in range(0, k):
      interclusterDistance += counts[i] * dist(listOfCenters[i], dataCenter)
   return interclusterDistance,
   
def compactness_separation(pop, numCenter, data, dataCenter):
   intraclusterDistance = 0
   interclusterDistance = 0
   listOfCenters = reshape_chromosome(pop, numCenter)
   k = len(listOfCenters)
   counts = [0] * k
   for i in range(0, len(data)):
         distanceList = []
         for j in range(0,len(listOfCenters)):
            distanceList.append(dist(data[i],listOfCenters[j]))
         minDistance = min(distanceList)
         intraclusterDistance += minDistance
         center = distanceList.index(minDistance)
         counts[center] += 1
   
   for i in range(0, k):
      interclusterDistance += counts[i] * dist(listOfCenters[i], dataCenter)
   return intraclusterDistance, interclusterDistance, 

def assign_memberships(pop):
   memberships = []
   centers = reshape_chromosome(pop, numCenter)
   for i in range(0, len(data)):
         distanceList = []
         for j in range(0,len(centers)):
            distanceList.append(dist(data[i],centers[j]))
         memberships.append(distanceList.index(min(distanceList)))
   return memberships

def select_best(pop):
   dbi_scores=[]
   for ind in pop:
      membership = assign_memberships(ind)
      dbi_scores.append(dbi(data, membership))
   index = dbi_scores.index(min(dbi_scores))
   best_ind = pop[index]
   return best_ind
   
def setup_so_ga(numAttribute, minNum, maxNum, indPb, mutPb, cxPb, numGen, sizePop, crowding, sizeTour):
   '''
   Creator-build new classes at run-time
   '''
   creator.create("FitnessMin", base.Fitness, weights = (-1.0,) )
   creator.create("Individual", list, fitness=creator.FitnessMin)

   '''
   Toolbox-create method
   '''
   toolbox = base.Toolbox()
   pool = multiprocessing.Pool()
   toolbox.register("map", pool.map)  #multi-processing
   #toolbox.register("map", map)      #single-processing
   toolbox.register("attr_bool", random.uniform, minNum, maxNum)
   toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, numAttribute * numCenter)
   toolbox.register("population", tools.initRepeat, list, toolbox.individual)
   toolbox.register("evaluate", compactness)
   toolbox.register('mate', tools.cxOnePoint)
   toolbox.register('mutate', tools.mutPolynomialBounded, eta = crowding, low = minNum, up = maxNum, indpb = mutPb)
   toolbox.register("select", tools.selTournament, tournsize=sizeTour)
   return toolbox
    
def setup_mo_ga(numAttribute, minNum, maxNum, indPb, mutPb, cxPb, numGen, sizePop, crowding, sizeTour, numCenter, data, dataCenter):
   creator.create("Fitness", base.Fitness, weights = (-1.0,1.0,) )
   creator.create("Individual", list, fitness=creator.Fitness)

   toolbox = base.Toolbox()
   ref_points = tools.uniform_reference_points(2, sizePop)
   pool = multiprocessing.Pool()
   toolbox.register("map", pool.map) #multi-processing
   #toolbox.register("map", map)     #single-processing
   toolbox.register("attr_bool", random.uniform, minNum, maxNum)
   toolbox.register("individual", tools.initRepeat, creator.Individual, toolbox.attr_bool, numAttribute * numCenter)
   toolbox.register("population", tools.initRepeat, list, toolbox.individual)
   toolbox.register("evaluate", compactness_separation, numCenter = numCenter, data =data, dataCenter = dataCenter)
   toolbox.register('mate', tools.cxOnePoint)
   toolbox.register('mutate', tools.mutPolynomialBounded, eta = crowding, low = minNum, up = maxNum, indpb = mutPb)
   toolbox.register("NSGAselect", tools.selNSGA3, ref_points=ref_points)
   toolbox.register("tourSelect", tools.selTournament, tournsize=sizeTour)
   return toolbox
 
def run_so_gakmeans(datasets, numCenters, indPb = params['indPb'], mutPb = params['mutPb'], cxPb = params['cxPb'], numGen = params['numGen'], sizePop = params['sizePop'], crowding = params['crowding']):
   global numCenter
   global data
   numCenter = numCenters
   data = datasets.values.tolist()

   fitness = []
   sizeTour = int(0.2 * sizePop)  
   numAttribute = len(data[0])
   minNum = min([entry for sub in data for entry in sub])
   maxNum = max([entry for sub in data for entry in sub])

   toolbox = setup_so_ga(numAttribute, minNum, maxNum, indPb, mutPb, cxPb, numGen, sizePop, crowding, sizeTour) 
   '''
   Evolution starts
   '''
   start = time.time()
   print("-- Begin evolution --")
   pop = toolbox.population(n=sizePop)
   fitnesses = list(toolbox.map(toolbox.evaluate, pop))
   for ind, fit in zip(pop, fitnesses):
      ind.fitness.values = fit
   fits = [ind.fitness.values[0] for ind in pop]
   gen = 0
   for gen in range(1,numGen+1):
      print("-- Generation %i --" % gen)
      offspring = toolbox.select(pop, len(pop))
      
      offspring = list(toolbox.map(toolbox.clone, offspring))
      
      for child1, child2 in zip(offspring[::2], offspring[1::2]):
         if random.random() < cxPb:
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values

      for mutant in offspring:
         if random.random() < indPb:
            toolbox.mutate(mutant)
            del mutant.fitness.values
            
      invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
      fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
      for ind, fit in zip(invalid_ind, fitnesses):
         ind.fitness.values = fit
   
      pop[:] = offspring
      fits = [ind.fitness.values[0] for ind in pop]
      fitness.append(min(fits))
   end = time.time()
   print("-- End of (successful) evolution --")
   best_ind = tools.selBest(pop, 1)[0] 
   centers = reshape_chromosome(best_ind, numCenter)
   intraclusterDistance = best_ind.fitness.values[0]
   return co.ClusteringObject(centers, intraclusterDistance, end - start, fitness)

def run_mo_gakmeans(datasets, numCenters, indPb = params['indPb'], mutPb = params['mutPb'], cxPb = params['cxPb'], numGen = params['numGen'], sizePop = params['sizePop'], crowding = params['crowding']):
   global numCenter
   global data
   global dataCenter
   dataCenter = datasets.mean().to_list()
   numCenter = numCenters
   data = datasets.values.tolist()
   
   sizeTour = int(0.2 * sizePop)
   numAttribute = len(data[0])
   minNum = min([entry for sub in data for entry in sub])
   maxNum = max([entry for sub in data for entry in sub])

   toolbox = setup_mo_ga(numAttribute, minNum, maxNum, indPb, mutPb, cxPb, numGen, sizePop, crowding, sizeTour, numCenter, data, dataCenter)
   '''
   Evolution starts
   '''
   start = time.time()
   print("-- Begin evolution --")
   pop = toolbox.population(n=sizePop)
   fitnesses = list(toolbox.map(toolbox.evaluate, pop))
   for ind, fit in zip(pop, fitnesses):
      ind.fitness.values = fit
   fits = [ind.fitness.values for ind in pop] 
   
   for gen in range(1, numGen+1):
      print("-- Generation %i --" % gen)
         
      offspring = toolbox.tourSelect(pop, len(pop))
      offspring = list(toolbox.map(toolbox.clone, offspring))   
      for child1, child2 in zip(offspring[::2], offspring[1::2]):
         if random.random() < cxPb:
            toolbox.mate(child1, child2)
            del child1.fitness.values
            del child2.fitness.values

      for mutant in offspring:
         if random.random() < indPb:
            toolbox.mutate(mutant)
            del mutant.fitness.values
      
      invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
      fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
      for ind, fit in zip(invalid_ind, fitnesses):
         ind.fitness.values = fit
      
      pop = toolbox.NSGAselect(pop + offspring, sizePop)

   print("-- End of (successful) evolution --")
   best_ind = select_best(pop)
   end = time.time()
   centers = reshape_chromosome(best_ind, numCenter)
   
   intraclusterDistance = best_ind.fitness.values[0]
   interclusterDistance = best_ind.fitness.values[1] 
   return ClusteringObject(centers, [intraclusterDistance,interclusterDistance], end-start , [])

def run_kmeans(datasets, numCenters, maxiter = 100, random_state = None):
   start = time.time()
   kmeans = KMeans(n_clusters = numCenters, max_iter = maxiter).fit(datasets)
   end = time.time()
   centers = kmeans.cluster_centers_.tolist()
   intraclusterDistance = compute_total_intracluster_distance(datasets.values.tolist(), centers)
   return ClusteringObject(centers, intraclusterDistance, end - start, [])

def run_phenograph(datasets):
   memberships, graph, Q = phenograph.cluster(datasets)
   return memberships

def run_scanpy_clustering(data):
    adata = sc.AnnData(data)
    sc.pp.neighbors(adata)
    sc.tl.leiden(adata, resolution = 0.2)
    return adata.obs['leiden'].to_list()
