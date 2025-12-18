import sys
import collections
import dendropy
import treeswift as ts
import msprime
import tskit
from pathlib import Path
from dendropy.simulate import treesim
import daiquiri
import argparse

daiquiri.setup(level="DEBUG")


def scale_time(name, time):
    pop = name_to_pop[name]
    return (
        pop_to_stime[pop]
        + (time - pop_to_time[pop]) / pop_to_dur[pop] * pop_to_sdur[pop]
    )
    # return pop_to_stime[pop] + (time - pop_to_time[pop]) * pop_to_rmult[pop]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input-tree", required=True, type=Path, help="Input species tree."
    )
    parser.add_argument(
        "-p",
        "--pop-sizes",
        required=True,
        type=Path,
        help="Population to size mapping.",
    )
    parser.add_argument(
        "-s",
        "--sub-units",
        required=True,
        type=Path,
        help="Population to substitution units mapping.",
    )
    parser.add_argument(
        "-g",
        "--num-generations",
        required=True,
        type=Path,
        help="Population to number of generations mapping.",
    )
    parser.add_argument(
        "-o", "--output-path", required=True, type=Path, help="Output path."
    )
    parser.add_argument(
        "-l", "--seq-len", required=True, type=int, help="Sequence length."
    )
    parser.add_argument(
        "-r",
        "--recombination-rate",
        required=False,
        default=5e-9,
        type=float,
        help="Recombination rate.",
    )
    parser.add_argument("--seed", required=True, type=int, help="Random seed.")
    parser.add_argument(
        "--root-population",
        required=False,
        default="root",
        type=str,
        help="Name of the root population.",
    )
    args = parser.parse_args()

    print("--input-tree", args.input_tree, file=sys.stderr)
    print("--pop-sizes", args.pop_sizes, file=sys.stderr)
    print("--sub-units", args.sub_units, file=sys.stderr)
    print("--num-generations", args.num_generations, file=sys.stderr)
    print("--output-path", args.output_path, file=sys.stderr)
    print("--seq-len", args.seq_len, file=sys.stderr)
    print("--recombination-rate", args.recombination_rate, file=sys.stderr)
    print("--seed", args.seed, file=sys.stderr)
    print("--root-population", args.root_population, file=sys.stderr)

    t = dendropy.Tree.get(
        path=args.input_tree, schema="newick", preserve_underscores=True
    )

    samples = []
    height = t.max_distance_from_root()
    name_to_pop = {}
    for ix, node in enumerate(t.leaf_nodes()):
        # st = height - node.distance_from_root()
        # node.edge_length += st
        st = 0
        samples.append(
            msprime.SampleSet(1, ploidy=1, time=st, population=node.taxon.label)
        )
        name_to_pop[ix] = node.taxon.label

    initial_size = {}
    with open(args.pop_sizes, "r") as f:
        initial_size = dict(map(lambda x: x.strip().split("\t"), f.readlines()[1:]))
    initial_size = {k: float(v) / 2 for k, v in initial_size.items()}

    tree = t.as_string(
        schema="newick", suppress_rooting=True, unquoted_underscores=True
    )
    demography = msprime.Demography.from_species_tree(tree, initial_size)
    demography.sort_events()

    for ix, pop in enumerate(demography.populations):
        name_to_pop[ix] = pop.name

    tseq = msprime.sim_ancestry(
        samples=samples,
        ploidy=2,
        demography=demography,
        recombination_rate=args.recombination_rate,
        sequence_length=args.seq_len,
        random_seed=args.seed,
    )
    tables = tseq.dump_tables()

    nd_to_pop = {}
    for name, nd in enumerate(tables.nodes):
        nd_to_pop[name] = tables.populations[nd.population].metadata["name"]
    retained_pop_s = set(nd_to_pop.values())

    pop_to_sunit = {}
    with open(args.sub_units, "r") as f:
        for ix, l in enumerate(f):
            if ix > 0:
                pop, sunit = l.strip().split("\t")
                pop_to_sunit[pop] = float(sunit)

    pop_to_ngen = {}
    with open(args.num_generations, "r") as f:
        for ix, l in enumerate(f):
            if ix > 0:
                pop, ngen = l.strip().split("\t")
                ngen = eval(ngen)
                if ngen is not None:
                    pop_to_ngen[pop] = float(ngen)

    pop_to_rmult = {}
    pop_to_srate = {}
    mean_srate = 0
    ix = 0
    for pop, ngen in pop_to_ngen.items():
        pop_to_srate[pop] = pop_to_sunit[pop] / ngen
        if pop in retained_pop_s:
            mean_srate += pop_to_srate[pop]
            ix = ix + 1

    mean_srate /= ix
    for pop, srate in pop_to_srate.items():
        pop_to_rmult[pop] = srate / mean_srate
    pop_to_rmult["root"] = 1
    pop_to_srate["root"] = mean_srate

    tmain = ts.read_tree_newick(args.input_tree)
    pop_to_stime = {}
    pop_to_sdur = {}
    pop_to_time = {}
    pop_to_dur = {}

    for nd in tmain.traverse_levelorder(internal=True, leaves=True):
        l = nd.get_edge_length()
        if l is None:
            l = 0
        pop_to_time[nd.get_label()] = tmain.extract_subtree(nd).height() - l
        pop_to_dur[nd.get_label()] = l

    for nd in tmain.traverse_postorder(internal=True, leaves=False):
        for nd_child in nd.child_nodes():
            nd_child.set_edge_length(
                nd_child.get_edge_length() * pop_to_rmult[nd_child.get_label()]
            )

    h = tmain.height()
    for nd, dist in tmain.distances_from_root(
        leaves=True, internal=True, unlabeled=True, weighted=True
    ):
        hp = h - dist
        pop_to_stime[nd.get_label()] = hp
        pop_to_sdur[nd.get_label()] = nd.get_edge_length()
    pop_to_sdur[args.root_population] = 1
    pop_to_dur[args.root_population] = 1

    tables.nodes.time = [scale_time(nd.population, nd.time) for nd in tables.nodes]
    tables.sort()
    ntseq = tables.tree_sequence()
    ntseq.dump(args.output_path)
