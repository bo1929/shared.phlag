import os
import random
import argparse
from pathlib import Path


def count_trees(filepath):
    try:
        with open(filepath, "r") as f:
            line_count = sum(1 for line in f)
        return line_count
    except Exception as e:
        print(f"An error occurred: {e}")
        raise e


def simulate_independent_region(default_gtrees, discordant_gtrees, p, rate):
    default_gc = count_trees(default_gtrees)
    discordant_gc = count_trees(discordant_gtrees)
    assert default_gc == discordant_gc
    gc = default_gc = discordant_gc
    assert p < 0.5
    assert rate <= 1.0
    r = p / rate
    assert r < 1.0
    size = int(gc * r)
    dstart = random.randint(1, gc - size)
    dend = dstart + size
    vl = random.sample(
        range(dstart, dend), int((1.0 - rate) * (dend - dstart))
    )
    vl.sort()
    return dstart, dend, vl


def save_event(
    output_dir, default_gtrees, discordant_gtrees, dstart, dend, vl
):
    default_f = open(default_gtrees, "r")
    discordant_f = open(discordant_gtrees, "r")
    default_gtrees_l = default_f.readlines()
    discordant_gtrees_l = discordant_f.readlines()
    default_len = len(default_gtrees_l)
    discordant_len = len(discordant_gtrees_l)
    assert default_len == discordant_len

    with open(output_dir / "emission.gtrees", "w") as f:
        for i, gt in enumerate(default_gtrees_l):
            if i in vl:
                gt = default_gtrees_l[i]
            elif i >= dstart and i < dend:
                gt = discordant_gtrees_l[i]
            f.write(gt)

    default_f.close()
    discordant_f.close()
    with open(output_dir / "info.txt", "w") as f:
        f.write(f"type: single_independent\n")
        f.write(f"start: {dstart}\n")
        f.write(f"end: {dend}\n")
        f.write(f"default_gtrees: {default_gtrees}\n")
        f.write(f"discordant_gtrees: {discordant_gtrees}\n")
        f.write(f"p: {(dend - dstart)/float(default_len)}\n")
        f.write(f"r: {1.0-len(vl)/float(dend-dstart)}\n")
        f.write(f"v: {vl}")


def main(args):
    output_dir = args.output_dir
    default_gtrees = args.default_gene_trees
    discordant_gtrees = args.discordant_gene_trees
    p = args.discordant_portion
    r = args.rate

    os.makedirs(output_dir, exist_ok=True)
    output_dir = Path(output_dir)

    dstart, dend, vl = simulate_independent_region(
        default_gtrees, discordant_gtrees, p, r
    )
    save_event(
        output_dir,
        default_gtrees,
        discordant_gtrees,
        dstart,
        dend,
        vl,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-x",
        "--default-gene-trees",
        required=True,
        help="Newick file for default gene (ordered) trees.",
    )
    parser.add_argument(
        "-y",
        "--discordant-gene-trees",
        required=True,
        help="Newick file for discordant (ordered) trees.",
    )
    parser.add_argument(
        "-p",
        "--discordant-portion",
        type=float,
        required=True,
        help="Portion of the discordant segment.",
    )
    parser.add_argument(
        "-r",
        "--rate",
        type=float,
        required=False,
        default=1.0,
        help="Desired rate of the discordant segment.",
    )
    parser.add_argument("-o", "--output-dir", required=True, help="Output directory.")
    args = parser.parse_args()

    main(args)
