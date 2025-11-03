import sys
import argparse
import pathlib
import numpy as np
from sklearn.metrics import confusion_matrix
import treeswift as ts
from scipy.spatial import distance
from collections import defaultdict


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


def get_labels(info_dict, stype):
    # assert info_dict["p"] < 0.5
    labels = np.zeros(info_dict.get("gc", DEFAULT_GC), dtype=int)
    if stype:
        for i in range(len(info_dict["start"])):
            labels[info_dict["start"][i] : info_dict["end"][i]] = 1
            # labels[info_dict["v"]] = 0
    else:
        labels[info_dict["start"] : info_dict["end"]] = 1
        # labels[info_dict["v"]] = 0
    return labels


def get_phlag_pred(input_file, info_dict, stype, alt):
    pred = np.zeros(info_dict.get("gc", DEFAULT_GC), dtype=int)
    pos = []
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("#"):
                continue
            else:
                pred = np.array(list(map(lambda x: int(x), line.strip().split(","))), dtype=int)
                if alt:
                    pass
                else:
                    break
    if len(pos) > 0:
        pred[pos] = 1
    return pred


def get_phylter_pred(input_file, info_dict, stype):
    a_dict = defaultdict(list)
    pred = np.zeros(info_dict.get("gc", DEFAULT_GC), dtype=int)
    pos = []
    with open(input_file, "r") as f:
        for line in f:
            if line.startswith("# Outlier gene(s) detected: 1"):
                pos = list(map(lambda x: int(x)-1, line[28:].split(";")))
            if line.startswith("#"):
                continue
            else:
                if stype:
                    o = line.strip().split("\t")
                    g = o[0]
                    s = o[1]
                    a_dict[s].append(int(g)-1)
                else:
                    break
    if stype:
        donor_pos = a_dict.get(info_dict['donor'], [])
        recipient_pos = a_dict.get(info_dict['recipient'], [])
        if donor_pos:
            pred[donor_pos] = 1
        if recipient_pos:
            pred[recipient_pos] = 1
    else:
        if len(pos) > 0:
            pred[pos] = 1
    return pred


def main(args):
    input_file = args.input_file
    info_file = args.info_file
    stype = args.stype
    revert = args.revert
    method = args.method
    alt = args.alt

    info_dict = read_info(info_file)
    true = get_labels(info_dict, stype)

    if method == "phlag":
        pred = get_phlag_pred(input_file, info_dict, stype, alt)
    elif method == "phylter":
        pred = get_phylter_pred(input_file, info_dict, stype)
    else:
        raise ValueError(f"Invalid method: {args.method}")
    if revert:
        if (np.sum(pred) > (pred.shape[0]/2)):
            pred = 1 - pred
            pred = pred.astype(int)
    tn, fp, fn, tp = confusion_matrix(true, pred).ravel().tolist()
    if stype:
        # r = info_dict.get('b', 0)
        b = len(info_dict["end"])
        p = info_dict['p']
        D = info_dict['donor']
        R = info_dict['recipient']
        ST = info_dict['gene_trees'].split("/")[-1].split("_")[0]
        print("TN\tFP\tFN\tTP\tMp\tMb\tD\tR\tST", file=sys.stderr)
        print(
            f"{tn}\t{fp}\t{fn}\t{tp}\t{p}\t{b}\t{D}\t{R}\t{ST}",
            end=None,
            file=sys.stdout,
        )
    else:
        r = info_dict.get('r', 0)
        p = info_dict['p']
        print("TN\tFP\tFN\tTP\tMp\tMr", file=sys.stderr)
        print(
            f"{tn}\t{fp}\t{fn}\t{tp}\t{p}\t{r}",
            end=None,
            file=sys.stdout,
        )


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--input-file", type=pathlib.Path, required=True)
    parser.add_argument("-y", "--info-file", type=pathlib.Path, required=True)
    parser.add_argument("--method", type=str, required=False, choices=["phlag", "phylter"])
    parser.add_argument("--stype", action="store_true", required=False)
    parser.add_argument("--alt", action="store_true", required=False)
    parser.add_argument("--revert", action="store_true", required=False)
    args = parser.parse_args()
    main(args)
