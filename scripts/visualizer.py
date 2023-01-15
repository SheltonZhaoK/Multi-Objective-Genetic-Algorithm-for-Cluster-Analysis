import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

matplotlib.rcParams['svg.fonttype'] = 'none'

def plot_data(data, fileName):
   plt.figure(fileName)
   x,y = data.to_numpy().T
   plt.scatter(x,y)
   plt.xlabel("umap_1")
   plt.ylabel("umap_2")
   plt.savefig(fileName)

def plot_data_clustering(data, memberships, dataType ,fileName):
   data['memberships'] = memberships
   if dataType == 'umap':
      sns.FacetGrid(data, hue="memberships").map(plt.scatter, "UMAP_1", "UMAP_2").add_legend()
   else:
      sns.FacetGrid(data, hue="memberships").map(plt.scatter, "tSNE_1", "tSNE_2").add_legend()
   plt.savefig(fileName)
