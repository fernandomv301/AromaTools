import sys
import numpy as np
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning, message='.*matmul.*')

from esipy.make_aoms import make_aoms
from esipy.indicators import compute_iring, sequential_mci, compute_pdi, compute_flu, compute_av1245
from esipy.tools import find_di, format_partition

# Ring search defaults taken from ESIpy
RING_MIN_LENGTH = 6
RING_MAX_LENGTH = 12
RING_THRESHOLD = 0.3     # minimum DI (electrons) for two atoms to be connected
MCI_MAX_RING_SIZE = 12   # ESIpy only computes the MCI for rings up to 12 atoms
PDI_RING_SIZE = 6        # ESIpy only computes the PDI for 6-membered rings
AV1245_MIN_RING_SIZE = 6 # ESIpy only computes AV1245/AVmin for rings of 6 atoms or more

PARTITION_NAMES = {'NAO': 'nao', 'IAO': 'iao', 'mate-LOWDIN': 'meta-lowdin'}
NOT_IMPLEMENTED_PARTITIONS = ('QTAIM', 'MULLIKEN', 'LOWDIN')

# ============================================================
# INPUT.txt and geometry
# ============================================================

def read_input(input_path):
    data = {}

    file = open(input_path, 'r')
    lines = file.readlines()
    file.close()

    for line in lines:
        hash_position = line.find('#')
        if hash_position != -1:
            line = line[:hash_position]

        line = line.strip()

        if line == '':
            continue

        words = line.split()
        key = words[0]
        value = ' '.join(words[1:])

        data[key] = value

    return data

def is_requested(input_data, key):
    value = input_data.get(key, '').strip().upper()
    return value == 'YES'

def read_lines(path):
    try:
        file = open(path, 'r')
    except OSError:
        sys.exit("Could not open '" + path + "'.")
    lines = file.readlines()
    file.close()
    return lines

def read_xyz(xyz_path):
    lines = read_lines(xyz_path)

    number_of_atoms = int(lines[0].strip())

    symbols = []
    coordinates = []

    for i in range(2, 2 + number_of_atoms):
        parts = lines[i].split()
        symbols.append(parts[0])
        coordinates.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return symbols, coordinates

# ============================================================
# PySCF calculation
# ============================================================

def run_pyscf(input_data, symbols, coordinates, charge, multiplicity):
    from pyscf import gto, dft

    functional_basis = input_data['FUNCTIONAL/BASIS-SET'].split('/')
    functional = functional_basis[0].replace('-', '').lower()
    basis = functional_basis[1].lower()

    atom_text = ''
    for symbol, xyz in zip(symbols, coordinates):
        atom_text = atom_text + symbol + ' ' + str(xyz[0]) + ' ' + str(xyz[1]) + ' ' + str(xyz[2]) + '\n'

    mol = gto.Mole()
    mol.atom = atom_text
    mol.basis = basis
    mol.charge = charge
    mol.spin = multiplicity - 1
    mol.build()

    is_open_shell = multiplicity != 1
    mf = dft.UKS(mol) if is_open_shell else dft.RKS(mol)
    mf.xc = functional
    mf.kernel()

    return mol, mf, is_open_shell

# ============================================================
# Populations, charges and DI, straight from the AOM (works for any
# partition, not just NAO/IAO/meta-Lowdin -- it is the same formula ESIpy
# itself uses in rest.py/unrest.py's population analysis)
# ============================================================

def combine_channel_value(values, is_open_shell):
    if is_open_shell:
        return values[0] + values[1]
    return 2 * values[0]

def combine_channel_list(lists, is_open_shell):
    if is_open_shell:
        return list(np.add(lists[0], lists[1]))
    return [2 * v for v in lists[0]]

def charges_and_DI(mol, aoms_by_spin, is_open_shell):
    natoms = mol.natm

    N = np.array([combine_channel_value([np.trace(ch[A]) for ch in aoms_by_spin], is_open_shell)
                  for A in range(natoms)])
    Z = np.array([mol.atom_charge(A) for A in range(natoms)])
    charges = Z - N
    # print(charges)

    DI = np.zeros((natoms, natoms))
    for A in range(natoms):
        for B in range(A + 1, natoms):
            value = combine_channel_value([find_di(ch, A + 1, B + 1) for ch in aoms_by_spin], is_open_shell)
            DI[A][B] = value
            DI[B][A] = value

    return charges, DI

# ============================================================
# Ring detection (same algorithm as ESIpy, atoms are 1-indexed)
# ============================================================

def filter_connec(connec):
    if len(connec) < 2:
        return connec
    filtered_connec = {}
    for key, values in connec.items():
        if len(values) > 1:
            filtered_vals = [v for v in values if len(connec[v]) > 1]
            if filtered_vals:
                filtered_connec[key] = filtered_vals
    return filtered_connec

def build_connectivity(DI, threshold):
    natoms = len(DI)
    connec = {i: [] for i in range(1, natoms + 1)}

    for i in range(1, natoms):
        for j in range(i + 1, natoms + 1):
            if DI[i - 1][j - 1] >= threshold:
                connec[i].append(j)
                connec[j].append(i)

    return filter_connec({k: v for k, v in connec.items() if v})

def find_middle_nodes(connec):
    return [key for key, values in connec.items() if len(values) > 2]

def find_rings(connec, minlen, maxlen):
    def dfs(start):
        stack = [(start, [start])]
        while stack:
            v, path = stack.pop()
            if len(path) > maxlen:
                continue
            if (len(path) >= minlen and start in connec.get(path[-1], [])
                    and path[1] < path[-1] and path[0] == start):
                yield path
            for next_atom in connec.get(v, []):
                if next_atom not in path:
                    stack.append((next_atom, path + [next_atom]))

    def unique(path, all_paths):
        for shift in range(len(path)):
            rotated = np.roll(path, shift).tolist()
            if rotated in all_paths or rotated[::-1] in all_paths:
                return False
        return True

    all_paths = []
    starts = find_middle_nodes(connec)
    if not starts:
        starts = sorted(connec)
    for start in starts:
        for path in dfs(start):
            if unique(path, all_paths):
                all_paths.append(path)

    return all_paths

def get_rings(input_data, DI):
    ring_value = input_data.get('RING', '').strip()

    if ring_value != '':                                # Manual
        return [[int(number) for number in ring_value.split()]]

    connec = build_connectivity(DI, RING_THRESHOLD)     # Automatic, from the DIs
    if not connec:
        return []
    return find_rings(connec, RING_MIN_LENGTH, RING_MAX_LENGTH)

# ============================================================
# Aromaticity indices (Iring, MCI, PDI, FLU, AV1245/AVmin)
# Identical logic to mulliken_lowdin.py -- see the comment there for why
# the closed/open-shell combination rules are what they are.
# ============================================================

def nth_root(value, n):
    if value < 0:
        return -(abs(value) ** (1.0 / n))
    return value ** (1.0 / n)

def calc_iring_mci(aoms_by_spin, ring, partition, is_open_shell, want_iring, want_mci):
    result = {}

    if want_iring:
        irings = [compute_iring(ring, aom) for aom in aoms_by_spin]
        result['iring'] = combine_channel_value(irings, is_open_shell)
        if is_open_shell:
            result['iring_spin'] = irings

    if want_mci:
        if len(ring) > MCI_MAX_RING_SIZE:
            result['mci'] = None
        else:
            mcis = [sequential_mci(ring, aom, partition) for aom in aoms_by_spin]
            result['mci'] = combine_channel_value(mcis, is_open_shell)
            if is_open_shell:
                result['mci_spin'] = mcis

    return result

def calc_pdi(aoms_by_spin, ring, is_open_shell):
    if len(ring) != PDI_RING_SIZE:
        return None

    per_channel = [compute_pdi(ring, aom) for aom in aoms_by_spin]   # (value, [DI_14, DI_25, DI_36])
    values = [channel[0] for channel in per_channel]
    lists = [channel[1] for channel in per_channel]

    pdi = combine_channel_value(values, is_open_shell)
    pdi_list = combine_channel_list(lists, is_open_shell)
    return pdi, pdi_list

def calc_av1245(aoms_by_spin, ring, partition, is_open_shell):
    if len(ring) < AV1245_MIN_RING_SIZE:
        return None

    per_channel = [compute_av1245(ring, aom, partition) for aom in aoms_by_spin]   # (mean, avmin, products)
    products_per_channel = [channel[2] for channel in per_channel]

    products = combine_channel_list(products_per_channel, is_open_shell)
    av1245 = float(np.mean(products))
    avmin = min(products, key=abs)
    return av1245, avmin, products

def calc_flu(aoms_by_spin, ring, molinfo, partition, is_open_shell, flurefs=None):
    values = [compute_flu(ring, molinfo, aom, flurefs, partition) for aom in aoms_by_spin]
    if any(value is None for value in values):
        return None
    if is_open_shell:
        return values[0] + values[1]
    return values[0]

# ============================================================
# Output files: {molecule_name}_{partition}_{criterion}.csv
# ============================================================

def atom_label(symbols, atom):
    return symbols[atom - 1] + str(atom)

def ring_bonds(ring):
    bonds = []
    for i in range(len(ring)):
        A = ring[i]
        B = ring[(i + 1) % len(ring)]
        bonds.append((A, B))
    return bonds

def write_csv(molecule_name, partition, criterion, header, rows):
    file_name = molecule_name + '_' + partition + '_' + criterion + '.csv'
    file = open(file_name, 'w')

    file.write(header + '\n')
    for row in rows:
        file.write(','.join(str(item) for item in row) + '\n')

    file.close()
    return file_name

def write_info_txt(molecule_name, mol, is_open_shell, charge, multiplicity, functional_basis,
                    output_files, skipped, rings_by_partition):
    system_type = 'open-shell' if is_open_shell else 'closed-shell'
    equations = 'UHF/UKS' if is_open_shell else 'RHF/RKS'

    text = 'The ' + molecule_name + ' is ' + system_type + ', therefore, we use ' + equations + '.\n\n'
    text = text + 'Program: PySCF\n'
    text = text + 'Functional/basis-set: ' + functional_basis + '\n'
    text = text + 'Charge: ' + str(charge) + '\n'
    text = text + 'Multiplicity: ' + str(multiplicity) + '\n\n'

    text = text + 'Rings found:\n'
    for partition in rings_by_partition:
        rings = rings_by_partition[partition]
        if not rings:
            text = text + '  ' + partition.capitalize() + ': none found (ring size ' + \
                str(RING_MIN_LENGTH) + '-' + str(RING_MAX_LENGTH) + ', DI >= ' + str(RING_THRESHOLD) + \
                '; set RING in INPUT.txt to give it manually)\n'
        else:
            for number, ring in enumerate(rings, start=1):
                ring_text = '-'.join(str(atom) for atom in ring)
                text = text + '  ' + partition.capitalize() + ' ring ' + str(number) + ': ' + ring_text + '\n'
    text = text + '\n'

    text = text + 'Output files:\n'
    for file_name in output_files:
        text = text + '  - ' + file_name + '\n'

    if skipped:
        text = text + '\nRequested but not implemented (skipped): ' + ', '.join(skipped) + '\n'

    file = open('info.txt', 'w')
    file.write(text)
    file.close()
    print(text)

# ============================================================
# MAIN PROGRAM
# ============================================================

def analyse_partition(partition, input_data, mol, mf, is_open_shell, molecule_name,
                      want_di, want_iring, want_mci, want_pdi, want_flu, want_av1245,
                      output_files, rings_by_partition):
    symbols = [mol.atom_symbol(A) for A in range(mol.natm)]
    molinfo = {'symbols': symbols}
    print('=== ' + partition.capitalize() + ' partition ===')

    aom = make_aoms(mol, mf, partition, myhf=None, save=None)
    aoms_by_spin = aom if is_open_shell else [aom]

    charges, DI = charges_and_DI(mol, aoms_by_spin, is_open_shell)

    rows = []
    for A in range(len(symbols)):
        rows.append((symbols[A] + str(A + 1), charges[A]))
    output_files.append(write_csv(molecule_name, partition, 'charges',
                                  'Atom,' + partition.capitalize() + ' charges', rows))

    want_aromaticity = want_di or want_iring or want_mci or want_pdi or want_flu or want_av1245
    if not want_aromaticity:
        rings_by_partition[partition] = []
        return

    rings = get_rings(input_data, DI)
    rings_by_partition[partition] = rings

    if not rings:
        print('No ring found (ring size ' + str(RING_MIN_LENGTH) + '-' + str(RING_MAX_LENGTH) +
              ', DI >= ' + str(RING_THRESHOLD) + '). Set RING in INPUT.txt to give it manually.')
        return

    for number, ring in enumerate(rings, start=1):
        print('Ring ' + str(number) + ': ' + '-'.join(str(atom) for atom in ring))

    if want_di:
        rows = []
        for number, ring in enumerate(rings, start=1):
            for A, B in ring_bonds(ring):
                label = atom_label(symbols, A) + '-' + atom_label(symbols, B)
                rows.append((number, label, DI[A - 1][B - 1]))
        output_files.append(write_csv(molecule_name, partition, 'di',
                                      'Ring,Bond,delta(A,B)', rows))

    spin_header = ',alpha,beta' if is_open_shell else ''
    iring_rows, mci_rows, pdi_rows, flu_rows, av1245_rows = [], [], [], [], []

    for number, ring in enumerate(rings, start=1):
        ring_text = '-'.join(str(atom) for atom in ring)
        n = len(ring)

        result = calc_iring_mci(aoms_by_spin, ring, partition, is_open_shell, want_iring, want_mci)

        if want_iring:
            row = [number, ring_text, n, result['iring'], nth_root(result['iring'], n)]
            if is_open_shell:
                row = row + result['iring_spin']
            iring_rows.append(row)
            print('Iring ' + str(number) + ' = ' + str(round(result['iring'], 6)))

        if want_mci:
            if result['mci'] is None:
                print('MCI ' + str(number) + ' skipped: ring size ' + str(n) +
                      ' > ' + str(MCI_MAX_RING_SIZE))
            else:
                row = [number, ring_text, n, result['mci'], nth_root(result['mci'], n)]
                if is_open_shell:
                    row = row + result['mci_spin']
                mci_rows.append(row)
                print('MCI   ' + str(number) + ' = ' + str(round(result['mci'], 6)))

        if want_pdi:
            pdi_result = calc_pdi(aoms_by_spin, ring, is_open_shell)
            if pdi_result is None:
                print('PDI   ' + str(number) + ' skipped: ring size ' + str(n) +
                      ' != ' + str(PDI_RING_SIZE))
            else:
                pdi, pdi_list = pdi_result
                for k, (A, B) in enumerate(((ring[0], ring[3]), (ring[1], ring[4]), (ring[2], ring[5]))):
                    label = atom_label(symbols, A) + '-' + atom_label(symbols, B)
                    pdi_rows.append([number, pdi, label, pdi_list[k]])
                print('PDI   ' + str(number) + ' = ' + str(round(pdi, 6)))

        if want_flu:
            flu = calc_flu(aoms_by_spin, ring, molinfo, partition, is_open_shell)
            if flu is None:
                print('FLU   ' + str(number) + ' skipped: missing reference for a bond type in the ring')
            else:
                flu_rows.append([number, ring_text, n, flu])
                print('FLU   ' + str(number) + ' = ' + str(round(flu, 6)))

        if want_av1245:
            av_result = calc_av1245(aoms_by_spin, ring, partition, is_open_shell)
            if av_result is None:
                print('AV1245 ' + str(number) + ' skipped: ring size ' + str(n) +
                      ' < ' + str(AV1245_MIN_RING_SIZE))
            else:
                av1245, avmin, products = av_result
                for k in range(n):
                    quadruplet = (ring[k % n], ring[(k + 1) % n], ring[(k + 3) % n], ring[(k + 4) % n])
                    label = '-'.join(atom_label(symbols, atom) for atom in quadruplet)
                    av1245_rows.append([number, av1245, avmin, label, products[k]])
                print('AV1245 ' + str(number) + ' = ' + str(round(av1245, 6)) +
                      ' | AVmin ' + str(number) + ' = ' + str(round(avmin, 6)))

    if want_iring:
        output_files.append(write_csv(molecule_name, partition, 'iring',
                                      'Ring,Atoms,Size,Iring,Iring**(1/n)' + spin_header, iring_rows))
    if want_mci:
        output_files.append(write_csv(molecule_name, partition, 'mci',
                                      'Ring,Atoms,Size,MCI,MCI**(1/n)' + spin_header, mci_rows))
    if want_pdi:
        output_files.append(write_csv(molecule_name, partition, 'pdi',
                                      'Ring,PDI,Bond,DI', pdi_rows))
    if want_flu:
        output_files.append(write_csv(molecule_name, partition, 'flu',
                                      'Ring,Atoms,Size,FLU', flu_rows))
    if want_av1245:
        output_files.append(write_csv(molecule_name, partition, 'av1245',
                                      'Ring,AV1245,AVmin,Quadruplet,Value', av1245_rows))

def main():
    input_path = sys.argv[1] if len(sys.argv) > 1 else 'INPUT.txt'
    input_data = read_input(input_path)

    program = input_data.get('PROGRAM', '').strip().upper()
    if program not in ('', 'PYSCF'):
        sys.exit("This script only supports PROGRAM = PySCF (found '" + program +
                  "'). Use mulliken_lowdin.py for Gaussian/ORCA (Mulliken/Lowdin only).")

    partitions = [format_partition(PARTITION_NAMES[key]) for key in PARTITION_NAMES if is_requested(input_data, key)]
    skipped = [key for key in NOT_IMPLEMENTED_PARTITIONS if is_requested(input_data, key)]

    if not partitions:
        print('None of NAO, IAO, mate-LOWDIN is set to YES in INPUT.txt; nothing to do.')
        return
    if skipped:
        print('Not implemented in nao_iao_metalowdin.py, skipping: ' + ', '.join(skipped) +
              ' (use mulliken_lowdin.py for Mulliken/Lowdin)')

    for key in ('NAME_MOLECULE', 'FUNCTIONAL/BASIS-SET', 'CHARGE', 'MULTIPLICITY', 'MOLECULAR_COORDINARES'):
        if key not in input_data:
            sys.exit(key + ' is missing in INPUT.txt (required for PySCF).')

    molecule_name = input_data['NAME_MOLECULE']
    charge = int(input_data['CHARGE'])
    multiplicity = int(input_data['MULTIPLICITY'])
    symbols, coordinates = read_xyz(input_data['MOLECULAR_COORDINARES'])

    want_di = is_requested(input_data, 'DI')
    want_iring = is_requested(input_data, 'Iring')
    want_mci = is_requested(input_data, 'MCI')
    want_pdi = is_requested(input_data, 'PDI')
    want_flu = is_requested(input_data, 'FLU')
    want_av1245 = is_requested(input_data, 'AV1245')

    mol, mf, is_open_shell = run_pyscf(input_data, symbols, coordinates, charge, multiplicity)

    output_files = []
    rings_by_partition = {}
    for partition in partitions:
        analyse_partition(partition, input_data, mol, mf, is_open_shell, molecule_name,
                          want_di, want_iring, want_mci, want_pdi, want_flu, want_av1245,
                          output_files, rings_by_partition)

    write_info_txt(molecule_name, mol, is_open_shell, charge, multiplicity,
                    input_data['FUNCTIONAL/BASIS-SET'], output_files, skipped, rings_by_partition)


if __name__ == '__main__':
    main()


    
