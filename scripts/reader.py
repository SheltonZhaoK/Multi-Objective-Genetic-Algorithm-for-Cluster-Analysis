import subprocess
import io
import pandas as pd

def read_keel_data(filename):
   f = open(filename,'r')
   lines = f.readlines()
   data_index = -1
   col_names = []
   for line, i in zip(lines, range(len(lines))):
      if '@data' in line:
         data_index = i + 1
         break
      if '@inputs' in line:
         line = line.replace('@inputs', '').strip()
         col_names = line.split(',')
         
   data = pd.read_csv(filename, sep=",", header=None,on_bad_lines='skip',skiprows = data_index, engine='python')
   data.columns = col_names 
   if 'ID' in col_names:
      del data['ID']
   data.insert(0,'did',range(0,0+len(data)))
   f.close
   return data


def read_scrnaseq_data(filename):
   data = pd.read_csv(filename)
   cellId = data['Unnamed: 0'].to_list()
   del data['Unnamed: 0']
   return data, cellId

def read_scrnaseq_from_r():
   datasets = ["10X_LLU_A_cellranger2.0","10X_LLU_A_cellranger3.1","10X_LLU_A_umitools","10X_LLU_A_zumi",\
            "10X_LLU_B_cellranger2.0","10X_LLU_B_cellranger3.1","10X_LLU_B_umitools", "10X_LLU_B_zumi",\
            "10X_NCI_A_cellranger2.0", "10X_NCI_A_cellranger3.1", "10X_NCI_A_umitools", "10X_NCI_A_zumi",\
            "10X_NCI_B_cellranger2.0", "10X_NCI_B_cellranger3.1", "10X_NCI_B_umitools", "10X_NCI_B_zumi",\
            "10X_NCI_M_A_cellranger2.0", "10X_NCI_M_A_cellranger3.1", "10X_NCI_M_A_umitools", "10X_NCI_M_A_zumi",\
            "10X_NCI_M_A_umitools", "10X_NCI_M_A_zumi", "10X_NCI_M_B_cellranger2.0", "10X_NCI_M_B_cellranger3.1",\
            "10X_NCI_M_B_umitools", "10X_NCI_M_B_zumi", "C1_FDA_HT_A_featureCounts", "C1_FDA_HT_A_kallisto",\
            "C1_FDA_HT_A_rsem", "C1_FDA_HT_B_featureCounts", "C1_FDA_HT_B_kallisto", "C1_FDA_HT_B_rsem",\
            "C1_LLU_A_featureCounts", "C1_LLU_A_kallisto", "C1_LLU_A_rsem", "C1_LLU_B_featureCounts",\
            "C1_LLU_B_kallisto", "C1_LLU_B_rsem", "ICELL8_PE_A_featureCounts", "ICELL8_PE_A_kallisto",\
            "ICELL8_PE_A_rsem", "ICELL8_PE_B_featureCounts", "ICELL8_PE_B_kallisto", "ICELL8_PE_B_rsem",\
            "ICELL8_SE_A_featureCounts", "ICELL8_SE_A_kallisto", "ICELL8_SE_A_rsem", "ICELL8_SE_B_featureCounts",\
            "ICELL8_SE_B_kallisto", "ICELL8_SE_B_rsem"]
   
   command ='Rscript'
   path2script ='./embeddings.R'
   
   datasets_list = []
   datasets_cellId = []
   for dataset in datasets:
      cmd = [command, path2script, dataset]
      data = subprocess.check_output(cmd, universal_newlines=True)
      embeddings = pd.read_csv(io.StringIO(data), sep='\s+', lineterminator='\n', skiprows=4, skipfooter = 2,  header=None, engine='python')
      embeddings.columns = ['cell_id', 'umap_1', 'umap_2']
      cellID = embeddings['cell_id'].to_list()
      del embeddings['cell_id']
      
      datasets_cellId.append(cellID)
      datasets_list.append(embeddings)
      break 
   return datasets, datasets_cellId, datasets_list
