import sys
import argparse
import pathlib
import numpy as np
from sklearn.metrics import confusion_matrix
import treeswift as ts
from scipy.spatial import distance


DEFAULT_GC = 2000


def is_float(val):
    try:
        return float(val) == float(val)
    except (ValueError, TypeError):
        return False


def label_tree(tree):
    is_labeled = True
    i = 0
    labels = set()
    for node in tree.traverse_postorder():
        if node.is_leaf():
            node.set_edge_length(0)
            continue
        if not node.label or node.label in labels or is_float(node.label):
            is_labeled = False
            node.label = "I" + str(i)
            i += 1
        labels.add(node.label)
    return tree


def count_trees(filepath):
    try:
        with open(filepath, "r") as f:
            line_count = sum(1 for line in f)
        return line_count
    except Exception as e:
        print(f"An error occurred: {e}")
        raise e


def read_info(info_file):
    info_dict = {}
    with open(info_file, "r") as f:
        for line in f:
            k, v = line.strip().split(":")
            k, v = k.strip(), v.strip()
            try:
                info_dict[k] = eval(v)
            except (NameError, SyntaxError):
                info_dict[k] = v
    return info_dict


def get_labels(info_dict):
    # assert info_dict["p"] < 0.5
    labels = np.zeros(info_dict.get("gc", DEFAULT_GC), dtype=int)
    labels[info_dict["start"] : info_dict["end"]] = 1
    labels[info_dict["v"]] = 0
    return labels


def get_phylter_pred(input_file, info_dict):
    pred = np.zeros(info_dict.get("gc", DEFAULT_GC), dtype=int)
    pos = []
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("# Outlier gene(s) detected: 1"):
                pos = list(map(lambda x: int(x) - 1, line[28:].split(";")))
            if line.startswith("#"):
                continue
            else:
                break
    if len(pos) > 0:
        pred[pos] = 1
    return pred


def get_inference_stats(input_file):
    with open(input_file, "r") as f:
        for idx, line in enumerate(f):
            if line.startswith("#"):
                if idx == 0:
                    invocation = line[2:]
                elif idx == 1:
                    input_tree = ts.read_tree_newick(line[2:])
                elif idx == 2:
                    reestimated_tree = ts.read_tree_newick(line[2:])
                elif idx == 3:
                    clades_to_blen = eval(line[2:])
                    clades = list(clades_to_blen.keys())
                elif idx == 4:
                    default_emissions = eval(line[15:])
                elif idx == 5:
                    outlier_emissions = eval(line[15:])
                else:
                    continue
            else:
                break
    negblen = sum(
        [
            1 if blen <= 0 else 0
            for blen in reestimated_tree.branch_lengths(terminal=False, internal=True)
        ]
    )
    d = distance.jensenshannon(default_emissions, outlier_emissions, base=2)
    return d, negblen


def main(args):
    input_file = args.input_file
    info_file = args.info_file
    describe = args.describe
    method = args.method

    info_dict = read_info(info_file)
    true = get_labels(info_dict)

    if method == "gosh":
        pred = np.loadtxt(input_file, delimiter="\t", comments="#")
        tn, fp, fn, tp = confusion_matrix(true, pred).ravel().tolist()
        d, negblen = get_inference_stats(input_file)
        print("TN\tFP\tFN\tTP\tMp\tMr\tMv\tDjs\tCnb", file=sys.stderr)
        print(
            f"{tn}\t{fp}\t{fn}\t{tp}\t{info_dict['p']}\t{info_dict['r']}\t{np.mean(info_dict['v'])}\t{d:.4f}\t{negblen}",
            end=None,
            file=sys.stdout,
        )
    elif method == "phylter":
        pred = get_phylter_pred(input_file, info_dict)
        tn, fp, fn, tp = confusion_matrix(true, pred).ravel().tolist()
        print("TN\tFP\tFN\tTP\tMp\tMr\tMv", file=sys.stderr)
        print(
            f"{tn}\t{fp}\t{fn}\t{tp}\t{info_dict['p']}\t{info_dict['r']}\t{np.mean(info_dict['v'])}",
            end=None,
            file=sys.stdout,
        )
    else:
        raise ValueError(f"Invalid method: {args.method}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--input-file", type=pathlib.Path, required=True)
    parser.add_argument("-y", "--info-file", type=pathlib.Path, required=True)
    parser.add_argument(
        "--method", type=str, required=False, choices=["gosh", "phylter"]
    )
    parser.add_argument("--describe", action="store_true", required=False)
    args = parser.parse_args()
    main(args)
