import os
import argparse
from treeswift import *
import dendropy
import subprocess

from tqdm import tqdm


def is_float(val):
    try:
        return float(val) == float(val)
    except (ValueError, TypeError):
        return False


def __label_tree__(tree_obj):
    is_labeled = True
    i = 0
    labels = set()
    for node in tree_obj.traverse_postorder():
        if node.is_leaf():
            continue
        if not node.label or node.label in labels or is_float(node.label):
            is_labeled = False
            node.label = "I" + str(i)
            i += 1
        labels.add(node.label)
    return is_labeled


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument("-i", "--input-tree", required=True, help="Input Species Tree.")
    parser.add_argument("-g", "--gene-trees", required=True, help="Gene tree file.")
    parser.add_argument("-a", "--annot", required=False, help="Annotation file.")
    parser.add_argument(
        "-n", "--num-genes", default=1000, help="Number gene trees to simulate. [1000]"
    )
    parser.add_argument("-t", "--num-threads", default=8, help="Number of threads. [8]")
    parser.add_argument(
        "-o", "--outdir", required=False, default="./", help="Output directory. [CWD]"
    )
    args = parser.parse_args()

    # getting species tree in CU unit
    cu_tree_path = os.path.join(args.outdir, "cu_tree.tree")
    cmd = [
        "astral4_coalescent_unit",
        "-C",
        "-c",
        args.input_tree,
        "-i",
        args.gene_trees,
        "-o",
        cu_tree_path,
        "-t",
        args.num_threads,
    ]
    try:
        subprocess.run(
            cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL
        )
        print("CU tree completed successfully.")
    except subprocess.CalledProcessError as e:
        print(f"Error: Command failed with return code {e.returncode}")
    except FileNotFoundError:
        print("Error: 'astral4_coalescent_unit' not found in PATH.")

    with open(cu_tree_path, "r") as f:
        tree = f.read().strip().split("\n")[0]

    tree_obj = read_tree_newick(tree)
    is_labeled = __label_tree__(tree_obj)
    labels = [n.label for n in tree_obj.traverse_internal() if not n.is_root()]
    if not is_labeled:
        tree_obj.write_tree_newick(os.path.join(args.outdir, "labelled_cu_tree.tree"))

    # simulating gene trees
    tns = dendropy.TaxonNamespace()
    species_tree = dendropy.Tree.get(
        data=tree_obj.newick(), schema="newick", taxon_namespace=tns
    )
    gene_to_species_map = dendropy.TaxonNamespaceMapping.create_contained_taxon_mapping(
        containing_taxon_namespace=tns, num_contained=1
    )

    for i in tqdm(range(int(args.num_genes))):
        gene_tree = dendropy.simulate.treesim.contained_coalescent_tree(
            species_tree, gene_to_species_map
        )
        gt = read_tree_newick(gene_tree.as_string(schema="newick"))
        for n in gt.traverse_leaves():
            n.label = "_".join(n.label.split("_")[:-1])
        with open(os.path.join(args.outdir, "simulated.gtrees"), "a") as file:
            file.write(f"{gt.newick()}\n")


if __name__ == "__main__":
    main()
