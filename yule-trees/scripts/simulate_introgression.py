import os
import random
import argparse
import treeswift as ts
import numpy as np
from pathlib import Path


def count_trees(filepath):
    try:
        with open(filepath, "r") as f:
            line_count = sum(1 for line in f)
        return line_count
    except Exception as e:
        print(f"An error occurred: {e}")
        raise e


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


def switch_taxa(tree, recipient, donor):
    lbl2node = tree.label_to_node(selection="all")
    try:
        lbl2node[recipient].label = donor
    except KeyError:
        raise KeyError(f"Reciepient {recipient} not found in tree")
    try:
        lbl2node[donor].label = recipient
    except KeyError:
        raise KeyError(f"Donor {donor} not found in tree")
    return tree.newick()


def simulate_introgression_event(gene_trees_l, dstart, dend, donor, recipient):
    for i, j in zip(dstart, dend):
        for ix in range(i, j):
            gene_trees_l[ix] = switch_taxa(ts.read_tree_newick(gene_trees_l[ix].strip()), recipient, donor) + "\n"
    return gene_trees_l


def positive_integers_with_sum(n, total):
    ls = [0]
    rv = []
    while len(ls) < n:
        c = random.randint(0, total)
        ls.append(c)
    ls = sorted(ls)
    ls.append(total)
    for i in range(1, len(ls)):
        rv.append(ls[i] - ls[i - 1])
    return rv


# See: https://stackoverflow.com/questions/60171086/generate-multiple-sets-of-random-non-overlapping-intervals-within-a-range
# by https://peteroupc.github.io/
def f(z, w, a, b):
    rv = []
    indices = [x + w for x in positive_integers_with_sum(z, (b - a) - z * w)]
    start = a
    for i in indices:
        i_start = random.randint(start, i + start - w)
        rv.append([i_start, i_start + w - 1])
        start += i
    return rv


def simulate_separate_blocks(gene_trees, p, b):
    bb = np.random.poisson(b) + 1
    gc = count_trees(gene_trees)
    assert p < 0.5
    size = int(gc * p)
    size_b = size // (bb) + 1
    assert size_b > 2
    dstart, dend = list(zip(*f(bb, size_b, 0, gc - 1)))
    return dstart, dend


def save_event(output_dir, gene_trees_l, metadata):
    with open(output_dir / "emission.gtrees", "w") as f:
        for i, gt in enumerate(gene_trees_l):
            f.write(gt)
    with open(output_dir / "info.txt", "w") as f:
        f.write(f"type: separate_blocks\n")
        f.write(f"start: {metadata['start']}\n")
        f.write(f"end: {metadata['end']}\n")
        f.write(f"gc: {metadata['gc']}\n")
        f.write(f"gene_trees: {metadata['gene_trees']}\n")
        f.write(f"p: {metadata['p']}\n")
        f.write(f"b: {metadata['b']}\n")
        f.write(f"donor: {metadata['donor']}\n")
        f.write(f"recipient: {metadata['recipient']}\n")
        f.write(f"clade: {metadata['clade']}")


def get_target_clade(species_tree, donor, recipient):
    tree = label_tree(ts.read_tree_newick(species_tree))
    nd_to_lbl = tree.label_to_node(selection="all")
    nd_parent = tree.mrca({donor, recipient})
    donor = nd_to_lbl[donor]
    recipient = nd_to_lbl[recipient]
    if donor.get_parent().get_label() != nd_parent.get_label():
        nd_tmp = recipient
        recipient = donor
        donor = nd_tmp
    return recipient.get_parent().get_label()


def main(args):
    output_dir = args.output_dir
    gene_trees = args.gene_trees
    species_tree = args.species_tree
    recipient = args.recipient
    donor = args.donor
    p = args.discordant_portion
    b = args.num_blocks

    assert p < 0.5
    assert b >= 1

    os.makedirs(output_dir, exist_ok=True)

    dstart, dend = simulate_separate_blocks(gene_trees, p, b)

    with open(gene_trees, "r") as f:
        gene_trees_l = simulate_introgression_event(f.readlines(), dstart, dend, donor, recipient)
    metadata = {"type": "separate_blocks", "start": dstart, "end": dend, "gc": count_trees(gene_trees), "gene_trees": gene_trees, "p": p, "b": b, "recipient": recipient, "donor": donor}
    metadata["clade"] = get_target_clade(species_tree, donor, recipient)
    save_event(output_dir, gene_trees_l, metadata)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-s", "--species-tree", required=True, help="File for species trees to select the target branch.")
    parser.add_argument("-g", "--gene-trees", required=True, help="File for gene trees to modify.")
    parser.add_argument("-p", "--discordant-portion", required=True, type=float, help="Portion of the region of the genome with recombination suppression.")

    parser.add_argument("-b", "--num-blocks", required=False, type=int, default=1, help="Number of blocks to distribute the introgression event.")
    parser.add_argument("-r", "--recipient", required=True, help="Label of the recipient taxon.")
    parser.add_argument("-d", "--donor", required=True, help="Label of the donor taxon.")
    parser.add_argument("-o", "--output-dir", required=True, type=Path, help="Output directory.")
    args = parser.parse_args()

    main(args)
