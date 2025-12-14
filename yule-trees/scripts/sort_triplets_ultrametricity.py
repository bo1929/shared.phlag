import sys
import fileinput
import treeswift as ts
import numpy as np

if __name__ == "__main__":
    try:
        tree = ts.read_tree_newick(sys.argv[1])
        nd_to_lbl = tree.label_to_node(selection="all")
    except AttributeError:
        pass

    um_ins = []
    for ix, line in enumerate(fileinput.input()):
        if ix == 0:
            continue
        nd1, nd2 = line.strip().split(" ")
        ndp = tree.mrca({nd1, nd2})
        nd1 = nd_to_lbl[nd1]
        nd2 = nd_to_lbl[nd2]
        if nd1.get_parent().get_label() != ndp.get_label():
            nd_tmp = nd2
            nd2 = nd1
            nd1 = nd_tmp
        tree_sub = tree.extract_subtree(ndp)
        dist_l = [dist for nd, dist in tree_sub.distances_from_root(leaves=True)]
        # um_ins.append((line.strip().split(" "), sum(dist_l) / len(dist_l)))
        um_ins.append(((nd1.get_label(), nd2.get_label()), np.var(dist_l), nd2.get_parent().get_edge_length()))
    if um_ins:
        um_ins = list(sorted(um_ins, key=lambda x: x[1]))
        print("\n".join(list(map(lambda x: f"{x[0][0]}\t{x[0][1]}\t{x[1]}\t{x[2]}", um_ins))))
