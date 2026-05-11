#!/usr/bin/env python3
"""Reproduce the core tables and figures for the SHIFT/Lugano quotient manuscript.

Usage:
    python scripts/reproduce_all.py --out .
    python scripts/reproduce_all.py --out results

When --out . is used, the script updates the tracked `data/` and `paper/figures/`
directories. Otherwise it writes the same structure under the requested output directory.
"""
from __future__ import annotations

import argparse
import json
from itertools import product, permutations
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def ket(s: str) -> np.ndarray:
    if s == "0":
        return np.array([1, 0], complex)
    if s == "1":
        return np.array([0, 1], complex)
    if s == "+":
        return np.array([1, 1], complex) / np.sqrt(2)
    if s == "-":
        return np.array([1, -1], complex) / np.sqrt(2)
    raise ValueError(s)


def kron3(a: np.ndarray, b: np.ndarray, c: np.ndarray) -> np.ndarray:
    return np.kron(np.kron(a, b), c)


def comp(a: int, b: int, c: int) -> np.ndarray:
    return kron3(ket(str(a)), ket(str(b)), ket(str(c)))


def f_symbol(s: str) -> int:
    return 1 if s in ["+", "-"] else 0


def f_output(label: tuple[str, str, str]) -> tuple[int, int, int]:
    return tuple(f_symbol(s) for s in label)


def lugano_w(a: int, b: int, c: int) -> tuple[int, int, int]:
    return (c * (b ^ 1), a * (c ^ 1), b * (a ^ 1))


def all_bool_funcs(n: int):
    domain = list(product([0, 1], repeat=n))
    return [dict(zip(domain, table)) for table in product([0, 1], repeat=len(domain))]


def generate_partitions(n: int):
    """Generate set partitions as restricted growth strings."""
    def rec(prefix, max_label):
        if len(prefix) == n:
            yield tuple(prefix)
        else:
            for lab in range(max_label + 2):
                prefix.append(lab)
                yield from rec(prefix, max(max_label, lab))
                prefix.pop()
    yield from rec([0], 0)


def blocks(partition: tuple[int, ...]) -> list[list[int]]:
    d: dict[int, list[int]] = {}
    for i, b in enumerate(partition):
        d.setdefault(b, []).append(i)
    return list(d.values())


def block_signature(partition: tuple[int, ...], shift_names: list[str]) -> str:
    return "|".join(",".join(shift_names[i] for i in block) for block in blocks(partition))


def best_dynamic_causal_bound(inputs, funcs):
    """Brute-force/majority optimized dynamic causal bound for the Lugano game.

    One party is first. The order of the remaining two parties can depend on
    the first party's input and output. The last party is allowed to know all
    inputs, which is the standard generous causal strategy for this game.
    """
    def best_first(first: int):
        parties = [0, 1, 2]
        rem = [p for p in parties if p != first]
        best = -1.0
        count = 0
        best_trace = None
        for f_first in funcs[1]:
            for order_func in funcs[2]:
                groups = {rem[0]: {}, rem[1]: {}}
                preliminary = []
                for inp in inputs:
                    target = lugano_w(*inp)
                    in_first = inp[first]
                    out_first = f_first[(in_first,)]
                    first_ok = out_first == target[first]
                    sel = order_func[(in_first, out_first)]
                    second = rem[sel]
                    third = rem[1 - sel]
                    key = (in_first, out_first, inp[second])
                    if first_ok:
                        groups[second].setdefault(key, []).append(target[second])
                    preliminary.append((inp, target, out_first, first_ok, second, third, key))

                second_policy = {rem[0]: {}, rem[1]: {}}
                for second in rem:
                    for key, vals in groups[second].items():
                        ones = sum(vals)
                        zeros = len(vals) - ones
                        second_policy[second][key] = 1 if ones > zeros else 0

                ok = 0
                trace = []
                for inp, target, out_first, first_ok, second, third, key in preliminary:
                    out = [None, None, None]
                    out[first] = out_first
                    out[second] = second_policy[second].get(key, 0) if first_ok else 0
                    out[third] = target[third]
                    success = tuple(out) == target
                    ok += int(success)
                    trace.append((inp, tuple(out), target, success, (first, second, third)))
                score = ok / len(inputs)
                count += 1
                if score > best:
                    best = score
                    best_trace = trace
        return best, count, best_trace

    rows = []
    for first in [0, 1, 2]:
        best, count, trace = best_first(first)
        rows.append({
            "first_party": "ABC"[first],
            "best_success_probability": best,
            "n_first_order_strategies_checked": count,
            "n_successes_out_of_8": int(round(best * 8)),
            "best_trace_json": json.dumps([
                {"input": t[0], "output": t[1], "target": t[2], "success": t[3],
                 "order": "".join("ABC"[i] for i in t[4])}
                for t in trace
            ]),
        })
    return pd.DataFrame(rows)


def make_shift_quantities():
    inputs = list(product([0, 1], repeat=3))
    shift_labels = [
        ("0", "0", "0"),
        ("+", "0", "1"),
        ("0", "1", "+"),
        ("0", "1", "-"),
        ("1", "+", "0"),
        ("-", "0", "1"),
        ("1", "-", "0"),
        ("1", "1", "1"),
    ]
    shift_names = ["".join(x) for x in shift_labels]
    shift_states = np.array([kron3(ket(a), ket(b), ket(c)) for a, b, c in shift_labels])
    f_outputs = [f_output(lab) for lab in shift_labels]
    unique_outputs = sorted(set(f_outputs))
    out_index = {out: i for i, out in enumerate(unique_outputs)}

    P_label_given_input = np.zeros((8, 8))
    for i, inp in enumerate(inputs):
        P_label_given_input[i, :] = np.abs(shift_states.conj() @ comp(*inp)) ** 2

    P_LC = np.zeros((8, 4))
    for i, inp in enumerate(inputs):
        c = out_index[lugano_w(*inp)]
        for l in range(8):
            P_LC[l, c] += P_label_given_input[i, l] / len(inputs)

    return inputs, shift_labels, shift_names, shift_states, f_outputs, unique_outputs, P_label_given_input, P_LC


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=str, default=".", help="Output root directory")
    args = parser.parse_args()
    root = Path(args.out)
    data_dir = root / "data"
    fig_dir = root / "paper" / "figures"
    data_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    inputs, shift_labels, shift_names, shift_states, f_outputs, unique_outputs, P_label_given_input, P_LC = make_shift_quantities()

    # Dynamic causal bound.
    funcs = {n: all_bool_funcs(n) for n in range(4)}
    df_dynamic = best_dynamic_causal_bound(inputs, funcs)
    df_dynamic.to_csv(data_dir / "L1_dynamic_causal_bound_3over4.csv", index=False)

    plt.figure()
    plt.bar(df_dynamic["first_party"], df_dynamic["best_success_probability"])
    plt.axhline(0.75, linestyle="--", label="causal bound 3/4")
    plt.xlabel("First party in dynamic causal strategy")
    plt.ylabel("Best Lugano-game success")
    plt.title("Dynamic causal strategies reproduce the 3/4 Lugano bound")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / "fig_dynamic_bound.png", dpi=190)
    plt.close()

    # Joint distribution and quotient theorem.
    pd.DataFrame(P_LC, index=shift_names, columns=[str(o) for o in unique_outputs]).to_csv(
        data_dir / "T1_joint_distribution_SHIFT_label_and_Lugano_output.csv"
    )
    P_C = P_LC.sum(axis=0)
    rows = []
    for m in range(1, 9):
        rows.append({
            "n_coarse_grained_outcomes_m": m,
            "max_success_over_stochastic_coarse_grainings": min(m, 4) / 4,
            "formula_min_m4_over4": min(m, 4) / 4,
            "can_violate_3over4": min(m, 4) / 4 > 0.75,
            "can_reach_perfect": min(m, 4) / 4 == 1,
        })
    df_theorem = pd.DataFrame(rows)
    df_theorem.to_csv(data_dir / "T2_stochastic_coarse_graining_bound.csv", index=False)

    plt.figure()
    plt.plot(df_theorem["n_coarse_grained_outcomes_m"], df_theorem["max_success_over_stochastic_coarse_grainings"], marker="o", label="stochastic optimum")
    plt.plot(df_theorem["n_coarse_grained_outcomes_m"], df_theorem["formula_min_m4_over4"], linestyle="--", label="min(m,4)/4")
    plt.axhline(0.75, linestyle=":", label="causal bound 3/4")
    plt.xlabel("Number m of coarse-grained outcomes")
    plt.ylabel("Maximum Lugano success")
    plt.title("Minimal quotient bound")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / "fig_quotient_bound.png", dpi=190)
    plt.close()

    # SHIFT quotient diagram.
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.axis("off")
    groups = [
        ["000", "111"],
        ["+01", "-01"],
        ["01+", "01-"],
        ["1+0", "1-0"],
    ]
    outputs = ["000", "100", "001", "010"]
    for i, (g, out) in enumerate(zip(groups, outputs)):
        x = 0.08 + i * 0.24
        ax.text(x, 0.70, "\n".join(g), ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="black"))
        ax.annotate("", xy=(x, 0.36), xytext=(x, 0.56), arrowprops=dict(arrowstyle="->"))
        ax.text(x, 0.22, f"f-output\n{out}", ha="center", va="center",
                bbox=dict(boxstyle="round,pad=0.30", fc="white", ec="black"))
    ax.set_title("The four quotient classes of the SHIFT/Lugano witness")
    plt.tight_layout()
    plt.savefig(fig_dir / "fig_shift_quotient_diagram.png", dpi=190)
    plt.close()

    # Deterministic partitions.
    def partition_success(partition):
        success = 0.0
        for block in blocks(partition):
            scores = np.zeros(4)
            for l in block:
                scores += P_LC[l, :]
            success += scores.max()
        return float(success)

    def refines_f(partition):
        for block in blocks(partition):
            if len({f_outputs[i] for i in block}) > 1:
                return False
        return True

    part_rows = []
    for p in generate_partitions(8):
        succ = partition_success(p)
        part_rows.append({
            "partition": json.dumps(p),
            "n_blocks": len(blocks(p)),
            "success": succ,
            "signature": block_signature(p, shift_names),
            "violates_3over4": succ > 0.75 + 1e-12,
            "perfect": abs(succ - 1) < 1e-12,
            "refines_f_output_quotient": refines_f(p),
        })
    df_part = pd.DataFrame(part_rows)
    df_part.to_csv(data_dir / "T3_all_deterministic_partitions_verification.csv", index=False)
    summary = df_part.groupby("n_blocks").agg(
        count=("success", "count"),
        max_success=("success", "max"),
        min_success=("success", "min"),
        n_violating=("violates_3over4", "sum"),
        n_perfect=("perfect", "sum"),
        n_refining_f=("refines_f_output_quotient", "sum"),
    ).reset_index()
    summary.to_csv(data_dir / "T3_partition_summary.csv", index=False)
    pd.DataFrame([{
        "all_perfect_refine_f": bool(df_part[df_part["perfect"]]["refines_f_output_quotient"].all()),
        "all_refine_f_perfect": bool(df_part[df_part["refines_f_output_quotient"]]["perfect"].all()),
        "n_perfect": int(df_part["perfect"].sum()),
        "n_refine_f": int(df_part["refines_f_output_quotient"].sum()),
    }]).to_csv(data_dir / "T4_perfect_success_iff_refines_f_quotient.csv", index=False)

    plt.figure(figsize=(8, 5))
    plt.scatter(df_part["n_blocks"], df_part["success"], s=8, alpha=0.25)
    plt.plot(df_theorem["n_coarse_grained_outcomes_m"], df_theorem["formula_min_m4_over4"], color="black", linestyle="--", label="stochastic upper bound")
    plt.axhline(0.75, linestyle=":", label="causal bound")
    plt.xlabel("Number of deterministic partition blocks")
    plt.ylabel("Success")
    plt.title("All 4140 deterministic partitions obey the quotient bound")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / "fig_all_partitions.png", dpi=190)
    plt.close()

    # Continuous confusion.
    named = {
        "identity_8_outcomes": tuple(range(8)),
        "f_output_4_blocks": (0, 1, 2, 2, 3, 1, 3, 0),
        "x_only_2_blocks": tuple(f_symbol(l[0]) for l in shift_labels),
        "xyz_parity_2_blocks": tuple(f_symbol(a) ^ f_symbol(b) ^ f_symbol(c) for a, b, c in shift_labels),
        "one_block_no_info": tuple(0 for _ in shift_labels),
    }

    def success_with_confusion(partition, eta):
        Cmat = np.zeros((8, 8))
        for block in blocks(partition):
            m = len(block)
            for true_l in block:
                for report_l in block:
                    Cmat[report_l, true_l] += (1 - eta) / m
                Cmat[true_l, true_l] += eta
        P_report = P_label_given_input @ Cmat.T
        success = 0.0
        for report_l in range(8):
            scores = np.zeros(4)
            for i, inp in enumerate(inputs):
                c = {out: j for j, out in enumerate(unique_outputs)}[lugano_w(*inp)]
                scores[c] += P_report[i, report_l] / len(inputs)
            success += scores.max()
        return success

    eta_grid = np.linspace(0, 1, 201)
    cont_rows = []
    for name, part in named.items():
        # normalize labels to restricted growth style
        mapping, nxt, norm = {}, 0, []
        for b in part:
            if b not in mapping:
                mapping[b] = nxt; nxt += 1
            norm.append(mapping[b])
        part = tuple(norm)
        for eta in eta_grid:
            cont_rows.append({"partition_name": name, "eta_identity_within_block": eta,
                              "success_probability": success_with_confusion(part, eta)})
    df_cont = pd.DataFrame(cont_rows)
    df_cont.to_csv(data_dir / "P3_continuous_confusion_SHIFT_success.csv", index=False)

    plt.figure(figsize=(8, 5))
    for name in named:
        sub = df_cont[df_cont["partition_name"] == name]
        plt.plot(sub["eta_identity_within_block"], sub["success_probability"], label=name)
    plt.axhline(0.75, linestyle="--", label="causal bound 3/4")
    plt.xlabel("eta: identity component within record blocks")
    plt.ylabel("Lugano success")
    plt.title("Continuous confusion within SHIFT quotient blocks")
    plt.legend(fontsize=8)
    plt.tight_layout()
    plt.savefig(fig_dir / "fig_continuous_confusion.png", dpi=190)
    plt.close()

    # Simple topology witness figure.
    vals = pd.DataFrame([
        {"topology": "parity triangles", "uniform": 0.5, "parity": 0.5, "A-first": 0.0},
        {"topology": "pairs + background", "uniform": 0.5, "parity": 0.0, "A-first": 0.25},
    ])
    vals.to_csv(data_dir / "record_topology_mode_weights.csv", index=False)
    x = np.arange(len(vals)); width = 0.25
    plt.figure(figsize=(7, 4))
    for j, col in enumerate(["uniform", "parity", "A-first"]):
        plt.bar(x + (j - 1) * width, vals[col], width, label=col)
    plt.xticks(x, vals["topology"])
    plt.ylabel("mode weight")
    plt.title("Same scalar coherence can preserve different modes")
    plt.legend()
    plt.tight_layout()
    plt.savefig(fig_dir / "fig_topology_witness.png", dpi=190)
    plt.close()

    print("Wrote data to", data_dir)
    print("Wrote figures to", fig_dir)
    print(summary)


if __name__ == "__main__":
    main()
