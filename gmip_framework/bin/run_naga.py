import nbgwas
from nbgwas import Nbgwas
import pandas as pd
import igraph
import pickle
import argparse

def get_args(args=None):
  parser = argparse.ArgumentParser(description='Run NAGA: Network Assisted Genomic Analysisn')
  parser.add_argument("--YtrainFileName", help="")
  parser.add_argument("--XnetworkFileName", help="")
  parser.add_argument("--ScoreColname", default="neg_log_P", help="")
  parser.add_argument("--IDColname", default="ENSGID", help="")
  parser.add_argument("--Prefix", default="out", help="")
  return parser.parse_args(args)

########################################################################################################################
def main(config):
  with open(config['XnetworkFileName'], 'rb') as file:
    Net=pickle.load(file)
  #
  g=Nbgwas(network=Net)
  g.genes.from_file(config['YtrainFileName'], sep='\t', pval_col=config['ScoreColname'], name_col=config['IDColname'])
  g.genes.table.head()
  g.map_to_node_table(columns=config['ScoreColname'])
  g.diffuse(method='random_walk', node_attribute=config['ScoreColname'], result_name='PredictedScore')
  g.network.refresh_node_attributes()
  g.map_to_gene_table(columns=['PredictedScore'])
  g.genes.table.to_csv('%s.genes.nagaRes.tsv' % config['Prefix'], index=False, sep="\t")
  g.genes.table[[config['IDColname'], 'PredictedScore']].to_csv('%s.genes.nagaRes.lim.tsv' % config['Prefix'], index=False, sep="\t")
  

########################################################################################################################

if __name__ == '__main__':
  config=vars(get_args())
  main(config)

# config={
#   'YtrainFileName': '/scratch/04179/mshabb/GMIP2/2024_09_26/NAFLD/naga/outdir/2_preprocess/out_Ytrain_chrom10.tsv',
#   'XnetworkFileName': '/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files/naga/PCnet_matrix_full_for_NAGA.igraph.pkl',
#   'ScoreColname': 'neg_log_P', 'IDColname': 'ENSGID', 'Prefix': 'naga'
# }
