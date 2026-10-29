#!/usr/bin/python3
# -*- coding: utf-8 -*-

import math
import pandas as pd
import networkx as nx
import plotly.graph_objects as go
from ase.data import covalent_radii, atomic_numbers, chemical_symbols


BANNER = r"""
======================================================================
                           AROGEOMETRIC
            Geometric Criteria for Aromaticity Evaluation
======================================================================

AroGeometric implements widely used geometric descriptors of aromaticity:

 • HOMA93 (Harmonic Oscillator Model of Aromaticity with experimental data)
 • HOMAc (Harmonic Oscillator Model of Aromaticity with computational data)
 • HOMER (excited-state aromaticity)

Capabilities:

 • Automatic ring detection from molecular geometry
 • Quantification of bond-length equalization (BLA)
 • Analysis of mono- and polycyclic systems
 • Graph-based ring identification
 • Interactive 3D visualization of rings

Geometric indices complement magnetic criteria by providing structural
evidence of electron delocalization in cyclic molecules.

----------------------------------------------------------------------
Developed by:

   Fernando Martínez-Villarino
   Gabriel Merino

   Cinvestav – Mérida
   Centro de Investigación y de Estudios Avanzados del IPN

----------------------------------------------------------------------
Citation:

   Martínez-Villarino, F.; Orozco-Ic, M.; Merino, G.
   AromaTools 1.0, Cinvestav Mérida, Mérida, Mexico, 2025.

======================================================================
"""

R_opt = {'HOMA93': {'CC': 1.388, 'CN': 1.334, 'NN': 1.309, 'CO': 1.265},
         'HOMER':  {'CC': 1.437, 'CN': 1.39,  'NN': 1.375, 'CO': 1.379},
         'HOMAC':  {'CC': 1.392, 'CN': 1.333, 'NN': 1.318, 'CO': 1.315,
                    'SiSi': 2.163, 'CSi': 1.752}}

alpha = {'HOMA93': {'CC': 257.7,  'CN': 93.52,  'NN': 130.33, 'CO': 157.38},
         'HOMER':  {'CC': 950.74, 'CN': 506.43, 'NN': 187.36, 'CO': 164.96},
         'HOMAC':  {'CC': 153.37, 'CN': 111.83, 'NN': 98.99,  'CO': 335.16,
                    'SiSi': 325.6, 'CSi': 115.41}}


def build_graph(atoms, tol=0.45):
    G = nx.Graph()
    for idx, (sym, _, _, _) in enumerate(atoms):
        G.add_node(idx, element=sym)
    for i in range(len(atoms)):
        sym1, x1, y1, z1 = atoms[i]
        r1 = covalent_radii[atomic_numbers[sym1]]
        for j in range(i + 1, len(atoms)):
            sym2, x2, y2, z2 = atoms[j]
            r2 = covalent_radii[atomic_numbers[sym2]]
            dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
            if dist <= r1 + r2 + tol:
                bond_type = ''.join(sorted([sym1, sym2]))
                G.add_edge(i, j, length=dist, bond_type=bond_type)
    return G


def find_rings(G, max_size=8):
    rings = set()
    sizes = tuple(range(3, max_size + 1))

    def canon_ring(path):
        seq1 = list(path)
        seq2 = list(reversed(path))
        def normalize(seq):
            m = min(seq)
            i = seq.index(m)
            return tuple(seq[i:] + seq[:i])
        return min(normalize(seq1), normalize(seq2))

    def dfs(path, depth):
        if depth > max_size:
            return
        current = path[-1]
        for neighbor in G.neighbors(current):
            if neighbor == path[0] and len(path) in sizes:
                rings.add(canon_ring(path))
            elif neighbor not in path:
                dfs(path + [neighbor], depth + 1)

    for start in G.nodes():
        dfs([start], 1)

    return [list(r) for r in rings]


def order_ring_nodes(G, ring):
    ordered = [ring[0]]
    used = set(ordered)
    while len(ordered) < len(ring):
        current = ordered[-1]
        for neighbor in G.neighbors(current):
            if neighbor in ring and neighbor not in used:
                ordered.append(neighbor)
                used.add(neighbor)
                break
        else:
            raise ValueError(f"Could not consistently order ring nodes: {ring}")
    return ordered


def compute_index(G, atoms, index_type):
    rings = find_rings(G)
    ring_data = []
    bond_rows = []
    coords = {i: (x, y, z) for i, (_, x, y, z) in enumerate(atoms)}
    symbols = {i: sym for i, (sym, _, _, _) in enumerate(atoms)}

    if index_type == "HOMA":
        sub_indices = ["HOMA93", "HOMAC"]
    elif index_type == "HOMER":
        sub_indices = ["HOMER"]
    else:
        raise ValueError("Invalid index_type. Choose 'HOMA' or 'HOMER'.")

    for idx, ring in enumerate(rings, start=1):
        ordered_ring = order_ring_nodes(G, ring)
        ring_atoms = ordered_ring + [ordered_ring[0]]
        distances = []
        bond_info = []
        completo = True

        for i in range(len(ring)):
            a1, a2 = ring_atoms[i], ring_atoms[i + 1]
            if a1 not in coords or a2 not in coords:
                completo = False
                break
            x1, y1, z1 = coords[a1]
            x2, y2, z2 = coords[a2]
            sym1, sym2 = symbols[a1], symbols[a2]
            dist = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
            bond_type = ''.join(sorted([sym1, sym2]))
            bond_info.append({"type": bond_type, "length": dist})
            distances.append(dist)

        if not completo or len(distances) != len(ring):
            continue

        for b in bond_info:
            bond_rows.append({
                "Ring": idx,
                "Type_Bond": b["type"],
                "Length (Å)": round(b["length"], 4)
            })

        row = {"Ring": f"Ring {idx}", "Bonds": [round(d, 4) for d in distances]}

        for sub_index in sub_indices:
            contribs = []
            for b in bond_info:
                bt = b["type"]
                if bt in R_opt[sub_index]:
                    contribs.append(alpha[sub_index][bt] * (b["length"] - R_opt[sub_index][bt])**2)
                else:
                    print(f"[!] Bond type {bt} not supported for {sub_index}")
            if contribs:
                row[sub_index] = round(1.0 - sum(contribs) / len(contribs), 4)

        if len(distances) >= 3:
            n = len(distances)
            bla = sum(abs(distances[i] - distances[(i + 1) % n]) for i in range(n)) / n
            row["BLA"] = round(bla, 4)

        ring_data.append(row)

    df_main = pd.DataFrame(ring_data)
    df_main.index += 1
    df_bonds = pd.DataFrame(bond_rows)
    df_bonds.index += 1
    ring_indices = [
        (i, row.get("HOMAC", row.get("HOMA93", row.get("HOMER", 0.0))))
        for i, row in enumerate(ring_data)
    ]

    return rings, ring_indices, df_main, df_bonds


def get_3Dgraph(G, atoms, rings, ring_indices):
    if not isinstance(ring_indices, dict):
        ring_indices = dict(ring_indices)

    node_xyz = [(x, y, z) for _, x, y, z in atoms]
    symbols = [atoms[i][0] for i in range(len(atoms))]

    x_edges, y_edges, z_edges = [], [], []
    for u, v in G.edges():
        x_edges += [node_xyz[u][0], node_xyz[v][0], None]
        y_edges += [node_xyz[u][1], node_xyz[v][1], None]
        z_edges += [node_xyz[u][2], node_xyz[v][2], None]

    node_trace = go.Scatter3d(
        x=[c[0] for c in node_xyz],
        y=[c[1] for c in node_xyz],
        z=[c[2] for c in node_xyz],
        mode='markers+text',
        text=symbols,
        textposition='top center',
        marker=dict(size=5, color='lightblue'),
        name="Atoms"
    )
    edge_trace = go.Scatter3d(
        x=x_edges, y=y_edges, z=z_edges,
        mode='lines',
        line=dict(width=2, color='gray'),
        name="Bonds"
    )

    cx_list, cy_list, cz_list = [], [], []
    value_labels, ring_labels = [], []

    for i, ring in enumerate(rings):
        ring_coords = [node_xyz[a] for a in ring]
        cx = sum(c[0] for c in ring_coords) / len(ring_coords)
        cy = sum(c[1] for c in ring_coords) / len(ring_coords)
        cz = sum(c[2] for c in ring_coords) / len(ring_coords)
        cx_list.append(cx)
        cy_list.append(cy)
        cz_list.append(cz)
        value = ring_indices.get(i)
        value_labels.append(f"{value:.3f}" if value is not None else "")
        ring_labels.append(f"Ring {i + 1}")

    value_trace = go.Scatter3d(
        x=cx_list, y=cy_list, z=cz_list,
        mode='text',
        text=value_labels,
        textfont=dict(size=10, color='red'),
        textposition='middle center',
        name="Indices"
    )
    label_trace = go.Scatter3d(
        x=cx_list, y=cy_list, z=cz_list,
        mode='text',
        text=ring_labels,
        textfont=dict(size=12, color='black'),
        textposition='top center',
        hoverinfo='none',
        showlegend=False
    )

    fig = go.Figure(data=[edge_trace, node_trace, value_trace, label_trace])
    fig.update_layout(
        margin=dict(l=0, r=0, b=0, t=30),
        showlegend=False,
        scene=dict(
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
            zaxis=dict(visible=False)
        )
    )
    return fig


def extract_atoms_from_output(fileout):
    atoms = []
    with open(fileout, 'r') as f:
        lines = f.readlines()

    if fileout.endswith(".xyz"):
        for line in lines[2:]:
            parts = line.strip().split()
            if len(parts) == 4:
                symbol = parts[0]
                x, y, z = map(float, parts[1:])
                atoms.append((symbol, x, y, z))
    elif any("Entering Gaussian System" in line for line in lines):
        for i in range(len(lines) - 5):
            if "Standard orientation" in lines[i]:
                atoms = []
                j = i + 5
                while "---" not in lines[j]:
                    parts = lines[j].split()
                    symbol = chemical_symbols[int(parts[1])]
                    x, y, z = map(float, parts[3:6])
                    atoms.append((symbol, x, y, z))
                    j += 1
    elif any("O   R   C   A" in line for line in lines):
        start = False
        for line in lines:
            if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
                atoms = []
                start = True
                continue
            if start:
                if "---" in line or not line.strip():
                    start = False
                    continue
                parts = line.strip().split()
                if len(parts) == 4:
                    symbol = parts[0]
                    x, y, z = map(float, parts[1:])
                    atoms.append((symbol, x, y, z))

    return pd.DataFrame(atoms, columns=["Element", "X", "Y", "Z"])


def main():
    print(BANNER)

    fileout = input("Output file (Gaussian or ORCA): ").strip()
    index_type = input("Index type [HOMA/HOMER]: ").strip().upper()

    df_atoms = extract_atoms_from_output(fileout)
    atoms_list = [
        (row.Element, row.X, row.Y, row.Z)
        for row in df_atoms.itertuples(index=False)
        if row.Element != 'H'
    ]

    G = build_graph(atoms_list)
    rings, ring_indices, df_main, df_bonds = compute_index(G, atoms_list, index_type)

    if index_type == "HOMA":
        cols = [c for c in ["Ring", "HOMA93", "HOMAC", "BLA"] if c in df_main.columns]
    else:
        cols = [c for c in ["Ring", "HOMER", "BLA"] if c in df_main.columns]

    print("\n" + "=" * 50)
    print(df_main[cols].to_string())
    print("\n" + "=" * 50)
    print(df_bonds.to_string())

    base = fileout.rsplit(".", 1)[0]
    df_main[cols].to_csv(f"{base}_indices.csv", index=False)
    df_bonds.to_csv(f"{base}_bonds.csv", index=False)
    print(f"\nSaved: {base}_indices.csv")
    print(f"Saved: {base}_bonds.csv")

    show = input("\nShow 3D visualization? [y/n]: ").strip().lower()
    if show == 'y':
        fig = get_3Dgraph(G, atoms_list, rings, dict(ring_indices))
        fig.show()


if __name__ == "__main__":
    main()

