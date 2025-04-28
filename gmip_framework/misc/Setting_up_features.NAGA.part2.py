import nbgwas
from nbgwas import Nbgwas
import pandas as pd
import igraph
import pickle

AdjMat="/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files/naga/PCnet_matrix_full_for_NAGA.txt.gz"
AdjMat_df=pd.read_csv(AdjMat, sep='\t', header=0, index_col=0)
AdjMat_df.shape
AdjMat=AdjMat_df.values
Net=igraph.Graph.Adjacency((AdjMat > 0).tolist())
Net.es['weight']=AdjMat[AdjMat.nonzero()]
Net.vs['name']=AdjMat_df.columns
with open('/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files/naga/PCnet_matrix_full_for_NAGA.igraph.pkl', 'wb') as file:
  pickle.dump(Net, file)

AdjMat="/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files/netwas/netwas_global_top_matrix_full_for_NAGA.txt.gz"
AdjMat_df=pd.read_csv(AdjMat, sep='\t', header=0, index_col=0)
AdjMat_df.shape
AdjMat=AdjMat_df.values
Net=igraph.Graph.Adjacency((AdjMat > 0).tolist())
Net.es['weight']=AdjMat[AdjMat.nonzero()]
Net.vs['name']=AdjMat_df.columns
with open('/work2/04179/mshabb/common/04488464786465dwkljhefkhbefjuh/refdir/2_pops_files/netwas/netwas_global_top_matrix_full_for_NAGA.igraph.pkl', 'wb') as file:
  pickle.dump(Net, file)
