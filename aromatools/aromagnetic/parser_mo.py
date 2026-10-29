#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import numpy as np
from .amg_read_input import read_input

def find_outs(prefix, exts=("out", "log")):
    exts = tuple(e.lower() for e in exts)
    rx = re.compile(
        rf"^{re.escape(prefix)}(?:_(\d{{4}}))?\.(\w+)$",
        re.IGNORECASE
    )
    items = []
    for fn in os.listdir(os.getcwd()):
        m = rx.match(fn)
        if not m:
            continue
        ext = m.group(2).lower()
        if ext not in exts:
            continue
        chunk = m.group(1)
        idx = int(chunk) if chunk is not None else -1
        items.append((idx, fn))
    items.sort(key=lambda x: x[0])

    return [fn for _, fn in items]

def _expand_tokens_to_mos(tokens):
    mos = []
    for tok in tokens:
        tok = str(tok).strip()
        if not tok:
            continue
        if tok.upper() in ("TRUE", "FALSE", "YES", "NO", "ON", "OFF"):
            continue
        if "-" in tok:
            a, b = tok.split("-", 1)
            a_i, b_i = int(a), int(b)
            if b_i < a_i:
                a_i, b_i = b_i, a_i
            mos.extend(range(a_i, b_i + 1))
        else:
            mos.append(int(tok))
    return sorted(set(mos))


def parse_orbitals_from_cfg(cfg, block_key="ORBITALS"):
    orb = cfg.get(block_key, {})
    spec = {}
    others_labels = []
    for key, val in orb.items():
        key = str(key).strip()
        if not key.startswith("MO_"):
            continue
        label = key[3:]
        if val is None or isinstance(val, bool):
            tokens = []
        elif isinstance(val, (list, tuple)):
            tokens = [str(x) for x in val]
        else:
            tokens = str(val).split()
        clean = []
        for t in tokens:
            if isinstance(t, bool):
                continue
            s = str(t).strip()
            if not s:
                continue
            if s.upper() in ("TRUE", "FALSE", "YES", "NO", "ON", "OFF"):
                continue
            clean.append(s)
        mos = _expand_tokens_to_mos(clean)
        if mos:
            spec[label] = mos
        else:
            others_labels.append(label)

    return spec, others_labels

def build_mo_to_labels(spec, others_labels=None):
    mo2l = {}
    for label, mos in spec.items():
        for m in mos:
            mo2l.setdefault(m, []).append(label)
    return mo2l, list(others_labels or [])

def load_cfg(input_path="INPUT.txt"):
    cfg = read_input(input_path)
    if not cfg:
        raise FileNotFoundError(f"Missing {input_path}")
    return cfg

_re_bind_head = re.compile(
    r"Full Cartesian NMR shielding tensor \(ppm\) for atom\s+(?:gh|bq)\(\s*(\d+)\)",
    re.I,
)
_re_mo_row = re.compile(r"^\s*(\d+)\.\s+(.*)$")
_re_float = re.compile(
    r"[-+]?\d*\.\d+(?:[EeDd][-+]?\d+)?|[-+]?\d+(?:\.\d*)?(?:[EeDd][-+]?\d+)?"
)
_re_iso_head = re.compile(r"Summary of isotropic NMR chemical shielding", re.I)
_re_iso_canon = re.compile(r"Canonical MO contributions:\s*\(ppm\)", re.I)
_re_gh_col = re.compile(r"\b(?:gh|bq)\s+(\d+)\b", re.I)

def parse_gaussian_bind(prefix, expected, spec, others_labels=None, Bext=(0.0, 0.0, 1.0), da_max=30, component="z"):
    outs = sorted(find_outs(prefix))
    if not outs:
        raise FileNotFoundError(f"No Gaussian outputs found for prefix: {prefix}")
    expected = int(expected)
    Bext = np.array(Bext, float).reshape(3,)
    mo2l, others_labels = build_mo_to_labels(spec, others_labels)
    all_labels = list(spec.keys()) + list(others_labels)
    if not all_labels:
        raise ValueError("No orbital labels were defined in spec / others_labels.")
    comp = component.lower()
    if comp == "xyz":
        out = {
            lbl: (np.zeros(expected), np.zeros(expected), np.zeros(expected))
            for lbl in all_labels
        }
    else:
        if comp not in {"x", "y", "z"}:
            raise ValueError("component must be one of: 'x', 'y', 'z', 'xyz'")
        comp_idx = {"x": 0, "y": 1, "z": 2}[comp]
        out = {lbl: np.zeros(expected, float) for lbl in all_labels}
    touched = set()
    for file_idx, path in enumerate(outs):
        file_offset = file_idx * int(da_max)
        first_gh_local = None

        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            in_atom = False
            in_table = False
            acc = None
            gh_num = None
            for raw in fh:
                m = _re_bind_head.search(raw)
                if m:
                    gh_num = int(m.group(1))
                    if first_gh_local is None:
                        first_gh_local = gh_num
                    local_idx = gh_num - first_gh_local
                    idx = file_offset + local_idx
                    if not (0 <= idx < expected):
                        in_atom = False
                        in_table = False
                        acc = None
                        continue

                    acc = {lbl: np.zeros((3, 3), float) for lbl in all_labels}
                    in_atom = True
                    in_table = False
                    continue

                if not in_atom:
                    continue

                if "Canonical MO contributions" in raw:
                    in_table = True
                    continue

                if not in_table:
                    continue

                s = raw.strip()
                if not s:
                    continue
                if s.startswith("MO") and ("XX" in s and "YY" in s and "ZZ" in s):
                    continue
                if set(s) <= {"=", "-", " "}:
                    continue

                if s.startswith("Total"):
                    idx = file_offset + (gh_num - first_gh_local)
                    if 0 <= idx < expected:
                        touched.add(idx)
                        for lbl in all_labels:
                            Bind = -acc[lbl] @ Bext
                            if comp == "xyz":
                                bx, by, bz = out[lbl]
                                bx[idx], by[idx], bz[idx] = Bind
                            else:
                                out[lbl][idx] = Bind[comp_idx]
                    in_atom = False
                    in_table = False
                    acc = None
                    gh_num = None
                    continue

                mr = _re_mo_row.match(raw)
                if not mr:
                    continue

                mo = int(mr.group(1))

                labels = mo2l.get(mo)
                if not labels:
                    if others_labels:
                        labels = others_labels
                    else:
                        continue

                vals = [float(x.replace("D", "E").replace("d", "e")) for x in _re_float.findall(mr.group(2))]
                if len(vals) < 9:
                    continue

                sigma_mo = np.array(
                    [
                        [vals[0], vals[1], vals[2]],
                        [vals[3], vals[4], vals[5]],
                        [vals[6], vals[7], vals[8]],
                    ],
                    float,
                )

                for lbl in labels:
                    acc[lbl] += sigma_mo

    if not touched:
        raise ValueError(
            f"No MO tensor contributions were parsed from Gaussian outputs with prefix '{prefix}'."
        )

    if len(touched) != expected:
        print(
            f"WARNING: parsed tensor MO contributions for {len(touched)} / {expected} points "
            f"(prefix='{prefix}')."
        )

    return out

def parse_gaussian_nics(prefix, expected, spec, others_labels=None, da_max=30):
    outs = sorted(find_outs(prefix))
    if not outs:
        raise FileNotFoundError(f"No Gaussian outputs found for prefix: {prefix}")

    expected = int(expected)
    mo2l, others_labels = build_mo_to_labels(spec, others_labels)
    all_labels = list(spec.keys()) + list(others_labels)

    if not all_labels:
        raise ValueError("No orbital labels were defined in spec / others_labels.")
    out = {lbl: np.zeros(expected, float) for lbl in all_labels}
    touched = set()
    for file_idx, path in enumerate(outs):
        file_offset = file_idx * int(da_max)
        first_gh_local = None
        in_iso = False
        in_table = False
        gh_atom_nums = []
        gh_col_pos = []

        with open(path, "r", encoding="utf-8", errors="ignore") as fh:
            for raw in fh:
                if _re_iso_head.search(raw):
                    in_iso = True
                    in_table = False
                    gh_atom_nums = []
                    gh_col_pos = []
                    continue

                if not in_iso:
                    continue

                if _re_iso_canon.search(raw):
                    in_table = True
                    gh_atom_nums = []
                    gh_col_pos = []
                    continue

                if not in_table:
                    continue

                s = raw.strip()
                if not s:
                    continue

                if s.startswith("****") or s.lower().startswith("population") or s.lower().startswith("nbo"):
                    break

                if s.startswith("MO"):
                    toks = s.split()
                    gh_atom_nums = []
                    gh_col_pos = []

                    labels_line = []
                    i = 1
                    while i < len(toks):
                        if toks[i].lower() in ("gh", "bq") and i + 1 < len(toks):
                            labels_line.append(f"{toks[i]} {toks[i+1]}")
                            i += 2
                        else:
                            if i + 1 < len(toks) and toks[i+1].isdigit():
                                labels_line.append(f"{toks[i]} {toks[i+1]}")
                                i += 2
                            else:
                                labels_line.append(toks[i])
                                i += 1

                    for pos, lab in enumerate(labels_line):
                        m = _re_gh_col.search(lab)
                        if m:
                            gh_num = int(m.group(1))
                            gh_atom_nums.append(gh_num)
                            gh_col_pos.append(pos)

                    if gh_atom_nums and first_gh_local is None:
                        first_gh_local = min(gh_atom_nums)

                    continue

                if set(s) <= {"-", "=", " "}:
                    continue

                if s.startswith("Total"):
                    continue

                mr = _re_mo_row.match(raw)
                if not mr:
                    continue

                mo = int(mr.group(1))
                labels = mo2l.get(mo)
                if not labels:
                    if others_labels:
                        labels = others_labels
                    else:
                        continue

                vals = [float(x.replace("D", "E").replace("d", "e")) for x in _re_float.findall(mr.group(2))]
                if not gh_col_pos:
                    continue

                for gh_num, pos in zip(gh_atom_nums, gh_col_pos):
                    if pos >= len(vals):
                        continue
                    if first_gh_local is None:
                        first_gh_local = gh_num
                    idx = file_offset + (gh_num - first_gh_local)
                    if not (0 <= idx < expected):
                        continue

                    touched.add(idx)
                    v = vals[pos]
                    for lbl in labels:
                        out[lbl][idx] += v

    if not touched:
        raise ValueError(
            f"No isotropic MO contributions were parsed from Gaussian outputs with prefix '{prefix}'."
        )

    if len(touched) != expected:
        print(
            f"WARNING: parsed isotropic MO contributions for {len(touched)} / {expected} points "
            f"(prefix='{prefix}')."
        )

    for lbl in out:
        out[lbl] = -out[lbl]

    return out

