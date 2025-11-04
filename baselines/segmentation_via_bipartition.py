import sys
import argparse
import copy
import pathlib

import jax
import jax.numpy as jnp
import jax.random as jr
import treeswift as ts

from dynamax.hidden_markov_model import BernoulliHMM
from jaxtyping import Array, Float, Int
from typing import Union
from scipy.stats import entropy
from fastroot.MinVar import MV00_Tree


PRNGKeyT = Array
Scalar = Union[float, Float[Array, ""]]
IntScalar = Union[int, Int[Array, ""]]


def is_float(val):
    try:
        return float(val) == float(val)
    except (ValueError, TypeError):
        return False


def partition_tree(tree):
    for nd in tree.traverse_postorder(leaves=True, internal=True):
        nd.set_edge_length(1.0)
    mv_tree = MV00_Tree(ddpTree=tree)
    mv_tree.Reroot()
    root = mv_tree.get_root()
    subtrees = [mv_tree.ddpTree.extract_subtree(nd) for nd in root.child_nodes()]
    return subtrees


def label_tree(tree):
    is_labeled = True
    i = 0
    labels = set()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            continue
        if not node.label or (node.label in labels) or is_float(node.label):
            is_labeled = False
            node.set_label("I" + str(i))
            i += 1
        labels.add(node.label)
    return tree


class SegmentationViaBipartition:
    def __init__(self, args):
        self.args = args
        self.output_file = self.args.output_file

        # Read the species tree and the gene trees
        self.species_tree = label_tree(ts.read_tree_newick(self.args.species_tree))
        self.label_to_node = self.species_tree.label_to_node(selection="all")
        with open(self.args.gene_trees) as f:
            self.gene_trees = [ts.read_tree_newick(line.strip()) for line in f]
        self.gc = len(self.gene_trees)
        self.mask = jnp.ones(self.gc, dtype=bool)

        self.output_str = f"# {' '.join(sys.argv)}"
        self.clade_label = self.args.clade_label
        self.output_str += "\n# " + self.species_tree.newick()

        if self.clade_label:
            obs = self.get_bp(self.clade_label)
            obs = jnp.array(obs)[:, None]
        else:
            selected_clades = self.select_clades(self.args.num_clades)
            obs = []
            for clade_label in selected_clades:
                obs.append(self.get_bp(clade_label))
            obs = jnp.stack(obs, axis=1)
        hmm = BernoulliHMM(
            2,
            obs.shape[1],
            initial_probs_concentration=1.1,
            transition_matrix_concentration=1.1,
            transition_matrix_stickiness=0.0,
            emission_prior_concentration0=1.1,
            emission_prior_concentration1=1.1,
        )
        params, props = hmm.initialize()
        em_params, log_probs = hmm.fit_em(params, props, obs, num_iters=500)
        most_likely_states = hmm.most_likely_states(em_params, obs)
        self.output_str += "\n" + ",".join(map(lambda x: str(x), most_likely_states.astype(int).tolist()))
        self.output_str += "\n" + ",".join(map(lambda x: str(x), obs[:, :].sum(axis=1).tolist()))
        with open(self.output_file, "w") as f:
            f.write(self.output_str)

    def get_bp(self, clade_label):
        obs = []
        clbl_s = set([nd.label for nd in self.label_to_node[clade_label].traverse_leaves()])
        p_nd = self.label_to_node[clade_label].get_parent()

        for nd in p_nd.child_nodes():
            if nd.label != clade_label:
                tlbl = nd.label

        for ix in range(self.gc):
            glbl_s = [nd.label for nd in self.gene_trees[ix].mrca(clbl_s).traverse_leaves()]
            if tlbl in glbl_s:
                obs.append(1)
            else:
                obs.append(0)
        return jnp.array(obs)

    def select_best_clade(self, tree):
        nd_l, entropy_l = [], []
        for nd, _ in tree.distances_from_parent(internal=True, leaves=False):
            if (not (self.label_to_node.get(nd.label, ""))) or nd.is_root() or (nd.get_parent() is None):
                continue
            obs = self.get_bp(nd.label)
            p = sum(obs) / len(obs)
            ent = entropy([p, 1 - p], base=2)
            if is_float(ent):
                nd_l.append(nd)
                entropy_l.append(ent)
        if not nd_l:
            best_clade = None
        else:
            best_idx = jnp.argmax(jnp.array(entropy_l))
            best_clade = nd_l[best_idx].label
        return best_clade

    def select_clades(self, num_clades):
        partitions = [copy.deepcopy(self.species_tree)]
        selected_clades = []
        while len(selected_clades) < num_clades:
            partitions = [subtree for partition in partitions for subtree in partition_tree(partition)]
            selected_clades = [self.select_best_clade(tree) for tree in partitions if tree.num_nodes(leaves=True, internal=False) > 1]
            selected_clades = [label for label in selected_clades if label is not None]
        return selected_clades


def parse_arguments():
    parser = argparse.ArgumentParser(description="SegmentationViaBipartition")
    parser.add_argument("-s", "--species-tree", type=pathlib.Path, required=True, help="Path to species tree in Newick format")
    parser.add_argument("-g", "--gene-trees", type=pathlib.Path, required=True, help="Path to file for ordered gene trees (one Newick tree per line)")
    parser.add_argument("-c", "--clade-label", type=str, required=False, default="", help="Specific clade labels to analyze (auto-selected if not provided)")
    parser.add_argument("--num-clades", type=int, default=3, help="Minimum number of clades to automatically select (default: 3)")
    parser.add_argument("-o", "--output-file", type=pathlib.Path, required=True, help="Path to output file")
    return parser.parse_args()


def main():
    args = parse_arguments()
    svb = SegmentationViaBipartition(args)
    # svb.run()


if __name__ == "__main__":
    main()
