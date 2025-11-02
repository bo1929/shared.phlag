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


PRNGKeyT = Array
Scalar = Union[float, Float[Array, ""]]
IntScalar = Union[int, Int[Array, ""]]


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


class IntrogressionViaBipartition:
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

        obs = []
        clbl_s = set([nd.label for nd in self.label_to_node[self.clade_label].traverse_leaves()])
        p_nd = self.label_to_node[self.clade_label].get_parent()

        for nd in p_nd.child_nodes():
            if nd.label != self.clade_label:
                tlbl = nd.label

        for ix in range(self.gc):
            glbl_s = [nd.label for nd in self.gene_trees[ix].mrca(clbl_s).traverse_leaves()]
            if tlbl in glbl_s:
                obs.append(1)
            else:
                obs.append(0)

        obs = jnp.array(obs)[:, None]
        hmm = BernoulliHMM(
            2, 1, initial_probs_concentration=1.1, transition_matrix_concentration=1.1, transition_matrix_stickiness=0.0, emission_prior_concentration0=1.1, emission_prior_concentration1=1.1
        )
        params, props = hmm.initialize(key=jr.PRNGKey(7))
        em_params, log_probs = hmm.fit_em(params, props, obs, num_iters=500)
        most_likely_states = hmm.most_likely_states(em_params, obs)
        self.output_str += "\n" + ",".join(map(lambda x: str(x), most_likely_states.astype(int).tolist()))
        self.output_str += "\n" + ",".join(map(lambda x: str(x), obs[:, 0].tolist()))
        with open(self.output_file, "w") as f:
            f.write(self.output_str)


def parse_arguments():
    parser = argparse.ArgumentParser(description="IntrogressionViaBipartition")
    parser.add_argument("-s", "--species-tree", type=pathlib.Path, required=True, help="Path to species tree in Newick format")
    parser.add_argument("-g", "--gene-trees", type=pathlib.Path, required=True, help="Path to file for ordered gene trees (one Newick tree per line)")
    parser.add_argument("-c", "--clade-label", type=str, required=True, help="Specific clade labels to analyze (auto-selected if not provided)")
    parser.add_argument("-o", "--output-file", type=pathlib.Path, required=True, help="Path to output file")
    return parser.parse_args()


def main():
    args = parse_arguments()
    ivb = IntrogressionViaBipartition(args)
    # ivb.run()


if __name__ == "__main__":
    main()
