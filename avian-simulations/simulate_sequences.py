import sys
import collections
import dendropy
import msprime
import tskit
import argparse
import daiquiri
import treeswift as ts
import numpy as np

from pathlib import Path


daiquiri.setup(level="DEBUG")
RS=110

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(
        "-i", "--input-ts", type=Path, required=True, help="Input tree sequence object."
    )
    parser.add_argument(
        "-o", "--output-fasta", type=Path, required=True, help="Output FASTA alignment."
    )
    parser.add_argument(
        "-m",
        "--mutation-rate",
        type=float,
        required=False,
        default=1.9086e-08,
        help="Mutation rate.",
    )
    args = parser.parse_args()

    np.random.seed(RS)

    mrate = args.mutation_rate  # 1.908590331970381e-08  # mean mutation rate

    ts = tskit.load(args.input_ts)
    seq_len = int(ts.sequence_length)

    rate_gamma_shape = np.random.lognormal(1.5, 1)
    lrate_gamma_shape = rate_gamma_shape * 20

    rate_l = mrate * np.random.gamma(
        lrate_gamma_shape, 1 / lrate_gamma_shape, size=seq_len
    )
    lastp = 0
    while lastp < seq_len:
        rlen = np.random.geometric(1 / 1e4)  # mean local length is 1e4
        curp = min(lastp + rlen + 1, seq_len)
        rate_l[lastp:curp] *= np.random.gamma(rate_gamma_shape, 1 / rate_gamma_shape)
        lastp = curp
    rate_map = msprime.RateMap(position=range(seq_len + 1), rate=rate_l)

    eq_freq = np.random.dirichlet([36, 26, 28, 32])
    model = msprime.GTR(
        relative_rates=np.random.dirichlet([5, 15, 3, 6, 16, 5]),
        equilibrium_frequencies=eq_freq,
    )
    mts = msprime.sim_mutations(ts, rate=rate_map, model=model, random_seed=RS)

    dseq = np.random.choice(["A", "C", "G", "T"], seq_len, p=eq_freq)
    with open(args.output_fasta, "w") as f:
        for seq, pi in zip(mts.alignments(), mts.samples()):
            pop = mts.populations()[mts.nodes()[pi].population].metadata["name"]
            f.write(f">{pop}\n")
            nseq_l = [dseq[ix] if c == "N" else c for ix, c in enumerate(seq)]
            f.write("".join(nseq_l) + "\n")
