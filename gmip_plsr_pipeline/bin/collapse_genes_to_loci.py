########################################################################################################################

# This script is a preprocessing script for output of MAGMA
# It collapses genes which are close to each other into a single locus
# For every locus, gene with most significant P-value or largest Z-score is chosen
# The output is similar to MAGMA input with additional columns at the end.
# Above output file can be directly used by GIMP pipeline

########################################################################################################################

import argparse
import pandas as pd
from tqdm.auto import tqdm

########################################################################################################################

def get_args(args=None):
  parser = argparse.ArgumentParser(description='Run GMIP: Integrated Network Analysis Pipeline in python')
  parser.add_argument("--yfile", required=True,   type=str, help="Path to file which has response values. Expected to have ID column to be used as row identifiers, Score column with response values & Chromosome column with chromosome names in it. Identifiers in ID column should match identifiers in above xFile ID column. All values to be used as response should be numeric and NA values are not allowed.")
  parser.add_argument("--id",    default="GENE",  type=str, help="default:'GENE'  (based on magma output)")
  parser.add_argument("--chr",   default="CHR",   type=str, help="default:'CHR'   (based on magma output)")
  parser.add_argument("--start", default="START", type=str, help="default:'START' (based on magma output)")
  parser.add_argument("--end",   default="STOP",  type=str, help="default:'STOP'  (based on magma output)")
  parser.add_argument("--score", default="ZSTAT", type=str, help="default:'ZSTAT' (based on magma output)")
  parser.add_argument("--wd_sz", default=100000,  type=int, help="default:100000")
  return parser.parse_args(args)

########################################################################################################################

def main(config):
  print('\n##########\nINFO: config args used.')
  print(config)
  print('\n##########\nINFO: Reading yFile.')
  df=pd.read_csv(config["yfile"], delim_whitespace=True)
  print('\n##########\nINFO: Performing regions collapse.')
  df['LOCUS_START']=df[config['start']] - config['wd_sz']
  df['LOCUS_END']=df[config['end']] + config['wd_sz']
  df.loc[df['LOCUS_START'] < 0, 'LOCUS_START']=0
  df.loc[df['LOCUS_END'] < 0, 'LOCUS_END']=0
  # Sort genomic regions by chromosome & start position
  sorted_regions=df.sort_values(by=[config['chr'], 'LOCUS_START'])
  # Initialize empty dictionary to store selected genes
  selected_genes={}
  locus_gene_mapping=[]
  locus_id=1
  # Iterate through each region
  for _, current_gene in tqdm(sorted_regions.iterrows()):
    overlap_found=False
    # Iterate through selected genes dictionary
    for gene, selected_gene in selected_genes.items():
      # Check for overlap with the selected gene
      if current_gene[config['chr']]==selected_gene[config['chr']] and \
        current_gene['LOCUS_START']<=selected_gene['LOCUS_END'] and \
        selected_gene['LOCUS_START']<=current_gene['LOCUS_END']:
        overlap_found=True
        current_gene['LOCUS_ID'] = selected_gene['LOCUS_ID']
        # Update the selected gene if the current gene has a higher score
        if current_gene[config['score']] > selected_gene[config['score']]:
          selected_genes[gene]=current_gene
        break
    # If no overlap is found, add the current gene to the selected genes list
    if not overlap_found:
      selected_genes[current_gene[config['id']]]=current_gene
      current_gene['LOCUS_ID'] = f"locus{locus_id}"
      locus_id += 1
    locus_gene_mapping.append(current_gene)
  # Convert the selected genes dictionary back to a DataFrame
  selected_genes_df=pd.DataFrame(list(selected_genes.values()))
  locus_gene_mapping_df = pd.DataFrame(locus_gene_mapping)
  locus_gene_mapping_df.drop(['LOCUS_START', 'LOCUS_END'], axis=1, inplace=True)
  locus_gene_mapping_df.to_csv(config['yfile'] + '.loci', index=False, header=True, sep="\t")
  selected_genes_df.to_csv(config['yfile'] + '.collapsed.genes.out', index=False, header=True, sep="\t")

########################################################################################################################

if __name__ == '__main__':
  config=vars(get_args())
  main(config)

# config={
#  'yfile': '/scratch/04179/mshabb/final_GNAP_testings/mnt/efs/fs1/reference_files/GNAP_v4.0_2024_02_12/data/benchmarker_paper/results_2024_03_23_01/LDL_with_Effect.tbl/1_magma/LDL_with_Effect.tbl.magma.genes.out',
#  'id': 'GENE', 'chr': 'CHR', 'start': 'START', 'end': 'STOP', 'score': 'ZSTAT', 'wd_sz': 100000
#  }
