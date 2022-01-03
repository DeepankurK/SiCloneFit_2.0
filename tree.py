from Bio import Phylo
from io import StringIO

#sumtrees.py -s mcct --rooted --no-annotations -F newick --support-as-labels -b0 -o result.tre best_trees.tre 

f=open("/home/deepank/Downloads/Prof_Hamim/Data/P342/Try_01_final2/samples/best/result.tre","r")
tree=f.read()
tree=tree[5:]

handle = StringIO(tree)
graph = Phylo.read(handle, "newick")
file_out=open("/home/deepank/Downloads/Prof_Hamim/Data/P342/Try_01_final2/samples/best/tree_plot.txt","w")
Phylo.draw_ascii(graph,file=file_out)
