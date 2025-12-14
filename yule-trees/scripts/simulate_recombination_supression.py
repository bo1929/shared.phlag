import os
import argparse
import random
import treeswift as ts

from pathlib import Path


def count_trees(filepath):
    try:
        with open(filepath, "r") as f:
            line_count = sum(1 for line in f)
        return line_count
    except Exception as e:
        print(f"An error occurred: {e}")
        raise e


def get_random_internal_node(tree):
    nd_list = list(tree.traverse_postorder(leaves=False, internal=True))
    while True:
        random.shuffle(nd_list)
        nd = nd_list[0]
        if not nd.is_root() and len(nd.get_parent().child_nodes()) > 1:
            break
    return tree, nd


def drop_nni(nd):
    assert not nd.is_root()
    assert len(nd.get_parent().child_nodes()) > 1
    nd_parent = nd.get_parent()
    nd_child = random.sample(nd.child_nodes(), 1)[0]
    for nd_curr in nd_parent.child_nodes():
        if nd_curr != nd:
            nd_sibling = nd_curr
    nd_parent.remove_child(nd_sibling)
    nd_parent.add_child(nd_child)
    nd.remove_child(nd_child)
    nd.add_child(nd_sibling)
    nd_child.set_parent(nd_parent)
    nd_sibling.set_parent(nd)
    return nd


def simulate_independent_region(gene_trees, p):
    gc = count_trees(gene_trees)
    assert p < 0.5
    size = int(gc * p)
    dstart = random.randint(1, gc - size)
    dend = dstart + size
    return dstart, dend


def get_noncontigous_idx(dstart, dend, density):
    noncontigous_idx = random.sample(
        range(dstart, dend), int((1.0 - density) * (dend - dstart))
    )
    noncontigous_idx.sort()
    return noncontigous_idx


def fixed_suppressed_trees(gt1, gt2, i, dstart, dend, r):
    if i > dstart and i < dend:
        if random.random() < r:
            return gt2, 0
        else:
            return gt1, i
    else:
        return gt1, 0


def random_suppressed_trees(gt1, gt2, i, dstart, dend, r):
    if i > dstart and i < dend:
        if random.random() < r:
            ntree, nd = get_random_internal_node(ts.read_tree_newick(gt2))
            nd = drop_nni(nd)
            return ntree.newick() + "\n", 0
        else:
            return gt1, i
    else:
        return gt1, 0


def support_suppressed_trees(gt1, gt2, i, dstart, dend, r):
    c = 0
    if i > dstart and i < dend:
        ntree = ts.read_tree_newick(gt2)
        for nd in ntree.traverse_postorder(leaves=False, internal=True):
            if (
                (not nd.is_root())
                and (len(nd.get_parent().child_nodes()) > 1)
                and (nd.get_label() is not None)
            ):
                s = float(nd.get_label())
                if random.random() < (1.0 - (r * s)):
                    nd = drop_nni(nd)
                    c += 1
        return ntree.newick() + "\n", c
    else:
        return gt1, 0


def simulate_suppression_event(gene_trees_l, dstart, dend, option, r):
    vl = []
    for i, gt in enumerate(gene_trees_l):
        if option == "fixed":
            gt, v = fixed_suppressed_trees(
                gene_trees_l[i], gene_trees_l[dstart], i, dstart, dend, r
            )
        elif option == "random":
            gt, v = random_suppressed_trees(
                gene_trees_l[i], gene_trees_l[dstart], i, dstart, dend, r
            )
        elif option == "support":
            gt, v = support_suppressed_trees(
                gene_trees_l[i], gene_trees_l[dstart], i, dstart, dend, r
            )
        if v > 0:
            vl.append(v)
        gene_trees_l[i] = gt
    return gene_trees_l, vl


def save_event(output_dir, gene_trees_l, metadata):
    with open(output_dir / "emission.gtrees", "w") as f:
        for i, gt in enumerate(gene_trees_l):
            f.write(gt)
    with open(output_dir / "info.txt", "w") as f:
        f.write(f"type: single_independent\n")
        f.write(f"start: {metadata['start']}\n")
        f.write(f"end: {metadata['end']}\n")
        f.write(f"gc: {metadata['gc']}\n")
        f.write(f"gene_trees: {metadata['gene_trees']}\n")
        f.write(f"p: {metadata['p']}\n")
        f.write(f"r: {metadata['r']}\n")
        f.write(f"v: {metadata['v']}\n")
        f.write(f"option: {metadata['option']}")


def main(args):
    output_dir = args.output_dir
    gene_trees = args.gene_trees
    p = args.discordant_portion
    r = args.rate
    option = args.option

    assert p < 0.5
    assert r <= 1.0
    s = p / r
    assert s < 1.0

    os.makedirs(output_dir, exist_ok=True)

    dstart, dend = simulate_independent_region(gene_trees, p)

    with open(gene_trees, "r") as f:
        gene_trees_l, vl = simulate_suppression_event(
            f.readlines(), dstart, dend, option, r
        )
    metadata = {
        "type": "recombination_suppression",
        "start": dstart,
        "end": dend,
        "gc": count_trees(gene_trees),
        "gene_trees": gene_trees,
        "p": p,
        "r": r,
        "v": vl,
        "option": option,
    }
    save_event(output_dir, gene_trees_l, metadata)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-g",
        "--gene-trees",
        required=True,
        type=Path,
        help="File for gene trees to use.",
    )
    parser.add_argument(
        "-p",
        "--discordant-portion",
        required=True,
        type=float,
        help="Portion of the region of the genome with recombination suppression.",
    )
    parser.add_argument(
        "-r",
        "--rate",
        required=False,
        type=float,
        default=1.0,
        help="Rate of the gene tree mimicking recombination suppression.",
    )
    parser.add_argument(
        "-o", "--output-dir", required=True, type=Path, help="Output directory."
    )
    parser.add_argument(
        "--option",
        required=False,
        type=str,
        default="fixed",
        choices=["fixed", "random", "support"],
        help="Option for the trees in the suppressed region: fixed, random, support.",
    )
    args = parser.parse_args()

    main(args)
