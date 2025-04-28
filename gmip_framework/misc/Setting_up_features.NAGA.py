#mamba create -n GMIP_naga_v2 python==3.8.10 conda-forge::pip==24.0 -y 
#conda activate GMIP_naga_v2
#pip install naga-gwas

import nbgwas
from nbgwas import Nbgwas
import pandas as pd
import igraph

g=Nbgwas()
g.network.from_ndex('f93f402c-86d4-11e7-a10d-0ac135e8bacf') #PCNet
csr=g.network.adjacency_matrix
num_rows, num_cols = csr.shape
row_names=g.network.node_table['name'].to_list()
col_names=g.network.node_table['name'].to_list()

# Write directly to a file using non-zero elements only
with open("csr_matrix_with_default_labels.txt", "w") as f:
    # Write the column headers
    f.write("ID\t" + "\t".join(col_names) + "\n")
    for i in range(num_rows):
        f.write(row_names[i])  # Write the row name
        # Get the non-zero indices and data for the current row
        row_start = csr.indptr[i]
        row_end = csr.indptr[i + 1]
        row_data = csr.data[row_start:row_end]
        row_indices = csr.indices[row_start:row_end]
        # Write the matrix row, filling in zeros where there are no non-zero values
        current_col = 0
        for col in range(num_cols):
            if current_col < len(row_indices) and row_indices[current_col] == col:
                f.write(f"\t{row_data[current_col]}")
                current_col += 1
            else:
                f.write("\t0")
        f.write("\n")

