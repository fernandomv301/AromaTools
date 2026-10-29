import re
import sys
import numpy as np
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning, message='.*matmul.*')

from esipy.indicators import compute_iring, sequential_mci, compute_pdi, compute_flu, compute_av1245

# Ring search defaults taken from ESIpy
RING_MIN_LENGTH = 6
RING_MAX_LENGTH = 12
RING_THRESHOLD = 0.3     # minimum DI (electrons) for two atoms to be connected
MCI_MAX_RING_SIZE = 12   # ESIpy only computes the MCI for rings up to 12 atoms
PDI_RING_SIZE = 6        # ESIpy only computes the PDI for 6-membered rings
AV1245_MIN_RING_SIZE = 6 # ESIpy only computes AV1245/AVmin for rings of 6 atoms or more

ELEMENT_SYMBOLS = ('H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni '
                   'Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe').split()

PARTITION_NAMES = {'MULLIKEN': 'mulliken', 'LOWDIN': 'lowdin'}
NOT_IMPLEMENTED_PARTITIONS = ('QTAIM', 'META-LOWDIN', 'mate-LOWDIN', 'NAO', 'IAO')

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

def get_program(input_data):
    program = input_data.get('PROGRAM', '').strip().upper()
    if program not in ('PYSCF', 'GAUSSIAN', 'ORCA'):
        sys.exit("PROGRAM in INPUT.txt must be PYSCF, GAUSSIAN or ORCA (found '" + program + "').")
    return program

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

def element_charge(symbol):
    return ELEMENT_SYMBOLS.index(symbol.capitalize()) + 1

# ============================================================
# S and P from each program
# data = {S, P_list, atom_of_ao, Z, symbols, is_open_shell}
# P_list = [P_total] (closed-shell) or [P_alpha, P_beta]
# ============================================================

def get_pyscf_data(input_data, symbols, coordinates, charge, multiplicity):
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
    if is_open_shell:
        om = dft.UKS(mol)
    else:
        om = dft.RKS(mol)
    om.xc = functional
    om.kernel()

    S = mol.intor('int1e_ovlp')
    dm = om.make_rdm1()
    if is_open_shell:
        P_list = [dm[0], dm[1]]
    else:
        P_list = [dm]

    ao_labels = mol.ao_labels(fmt=None, base=0)
    atom_of_ao = np.array([label[0] for label in ao_labels])

    Z = [mol.atom_charge(A) for A in range(mol.natm)]

    return {'S': S, 'P_list': P_list, 'atom_of_ao': atom_of_ao, 'Z': Z,
            'symbols': [mol.atom_symbol(A) for A in range(mol.natm)],
            'is_open_shell': is_open_shell, 'charge': charge, 'multiplicity': multiplicity}

# ---------------- Gaussian ----------------

def find_last_line(lines, title, message):
    position = None
    for i in range(len(lines)):
        if lines[i].strip() == title:
            position = i
    if position is None:
        sys.exit(message)
    return position

def g16_number_of_basis_functions(lines):
    for line in lines:
        match = re.match(r'\s*(\d+)\s+basis functions,', line)
        if match:
            return int(match.group(1))
    sys.exit("Could not find the number of basis functions in the Gaussian output.")

def read_g16_overlap(lines, n):
    start = find_last_line(
        lines, '*** Overlap ***',
        "'*** Overlap ***' not found in the Gaussian output. Rerun g16 printing it (e.g. IOp(3/33=4)).")

    matrix = np.zeros((n, n))
    i = start + 1
    columns_done = 0

    while columns_done < n:
        columns = [int(x) for x in lines[i].split()]
        i = i + 1

        for row in range(columns[0], n + 1):
            values = lines[i].split()[1:]
            for k in range(len(values)):
                value = float(values[k].replace('D', 'E'))
                matrix[row - 1][columns[k] - 1] = value
                matrix[columns[k] - 1][row - 1] = value
            i = i + 1

        columns_done = columns_done + len(columns)

    return matrix

def read_g16_density(lines, n, title):

    start = find_last_line(
        lines, title,
        "'" + title + "' not found in the Gaussian output. Use POP=FULL (open-shell: Alpha/Beta blocks).")

    matrix = np.zeros((n, n))
    atom_of_ao = [0] * n
    atom_symbols = []
    current_atom = -1

    i = start + 1
    columns_done = 0

    while columns_done < n:
        columns = [int(x) for x in lines[i].split()]
        i = i + 1

        for row in range(columns[0], n + 1):
            parts = lines[i].split()
            number_of_values = min(row, columns[-1]) - columns[0] + 1
            head = parts[:len(parts) - number_of_values]
            values = parts[len(parts) - number_of_values:]

            if columns_done == 0:
                if len(head) >= 4 and head[1].isdigit():
                    current_atom = int(head[1]) - 1
                    atom_symbols.append(head[2])
                atom_of_ao[row - 1] = current_atom

            for k in range(len(values)):
                matrix[row - 1][columns[k] - 1] = float(values[k])
                matrix[columns[k] - 1][row - 1] = float(values[k])
            i = i + 1

        columns_done = columns_done + len(columns)

    return matrix, np.array(atom_of_ao), atom_symbols

def get_gaussian_charge_multiplicity(lines):
    pattern = re.compile(r'Charge\s*=\s*(-?\d+)\s+Multiplicity\s*=\s*(\d+)')
    for line in lines:
        match = pattern.search(line)
        if match:
            return int(match.group(1)), int(match.group(2))
    sys.exit("Could not find 'Charge = ... Multiplicity = ...' in the Gaussian output.")

def get_gaussian_data(output_path):
  
    lines = read_lines(output_path)
    n = g16_number_of_basis_functions(lines)
    S = read_g16_overlap(lines, n)
    charge, multiplicity = get_gaussian_charge_multiplicity(lines)

    is_open_shell = multiplicity != 1
    if is_open_shell:
        P_alpha, atom_of_ao, atom_symbols = read_g16_density(lines, n, 'Alpha Density Matrix:')
        P_beta, _, _ = read_g16_density(lines, n, 'Beta Density Matrix:')
        P_list = [P_alpha, P_beta]
    else:
        P_total, atom_of_ao, atom_symbols = read_g16_density(lines, n, 'Density Matrix:')
        P_list = [P_total]

    return {'S': S, 'P_list': P_list, 'atom_of_ao': atom_of_ao,
            'Z': [element_charge(s) for s in atom_symbols], 'symbols': atom_symbols,
            'is_open_shell': is_open_shell, 'charge': charge, 'multiplicity': multiplicity}

# ---------------- ORCA ----------------

def read_orca_matrix(lines, title, message):
    # Full square matrix, printed in column blocks (usually 6 columns)
    start = None
    for i in range(len(lines) - 1):
        if lines[i].strip() == title and lines[i + 1].strip() != '' and set(lines[i + 1].strip()) == {'-'}:
            start = i
    if start is None:
        sys.exit(message)

    values = {}
    columns = []
    i = start + 2
    while i < len(lines):
        parts = lines[i].split()

        if len(parts) == 0:
            i = i + 1
            continue

        if all(part.isdigit() for part in parts):
            columns = [int(part) for part in parts]
        elif len(columns) > 0 and parts[0].isdigit() and len(parts) == len(columns) + 1:
            for k in range(len(columns)):
                values[(int(parts[0]), columns[k])] = float(parts[k + 1])
        else:
            break

        i = i + 1

    n = max(key[0] for key in values) + 1
    matrix = np.zeros((n, n))
    for (row, column), value in values.items():
        matrix[row][column] = value

    return matrix

def read_orca_atom_of_ao(lines, n):
    start = None
    for i in range(len(lines)):
        if lines[i].strip() == 'BASIS SET INFORMATION':
            start = i
            break
    if start is None:
        sys.exit("'BASIS SET INFORMATION' not found in the ORCA output.")

    group_pattern = {}
    atom_group = {}
    atom_symbol = {}
    i = start + 1
    while i < len(lines) and 'AUXILIARY' not in lines[i]:
        group_match = re.match(r'\s*Group\s+(\d+)\s+Type\s+\w+\s*:.*contracted to\s+(\S+)\s+pattern', lines[i])
        if group_match:
            group_pattern[int(group_match.group(1))] = group_match.group(2)

        atom_match = re.match(r'\s*Atom\s+(\d+)([A-Za-z]+)\s+basis set group\s*=>\s*(\d+)', lines[i])
        if atom_match:
            atom_group[int(atom_match.group(1))] = int(atom_match.group(3))
            atom_symbol[int(atom_match.group(1))] = atom_match.group(2)

        i = i + 1

    natoms = len(atom_group)
    symbols = [atom_symbol[A] for A in range(natoms)]

    spherical = {'s': 1, 'p': 3, 'd': 5, 'f': 7, 'g': 9, 'h': 11}
    cartesian = {'s': 1, 'p': 3, 'd': 6, 'f': 10, 'g': 15, 'h': 21}

    for functions_per_shell in (spherical, cartesian):
        counts = []
        for A in range(natoms):
            pattern = group_pattern[atom_group[A]]
            count = 0
            for number, shell in re.findall(r'(\d+)([spdfgh])', pattern):
                count = count + int(number) * functions_per_shell[shell]
            counts.append(count)
        if sum(counts) == n:
            atom_of_ao = np.repeat(np.arange(natoms), counts)
            return atom_of_ao, symbols

    sys.exit("Could not assign the ORCA basis functions to atoms (count does not match the matrix size).")

def get_orca_charge_multiplicity(lines):
    charge = None
    multiplicity = None
    charge_pattern = re.compile(r'Total Charge\s+Charge\s*\.+\s*(-?\d+)')
    mult_pattern = re.compile(r'Multiplicity\s+Mult\s*\.+\s*(\d+)')

    for line in lines:
        charge_match = charge_pattern.search(line)
        if charge_match:
            charge = int(charge_match.group(1))
        mult_match = mult_pattern.search(line)
        if mult_match:
            multiplicity = int(mult_match.group(1))

    if charge is None or multiplicity is None:
        sys.exit("Could not find the 'Total Charge'/'Multiplicity' lines in the ORCA output "
                  "(expected in the General Settings block, before the SCF iterations).")

    return charge, multiplicity

def get_orca_data(output_path):
    lines = read_lines(output_path)

    S = read_orca_matrix(
        lines, 'OVERLAP MATRIX',
        "OVERLAP MATRIX not found in the ORCA output. Add: %output Print[ P_Overlap ] 1 end")
    D = read_orca_matrix(
        lines, 'DENSITY',
        "DENSITY not found in the ORCA output. Add: %output Print[ P_Density ] 1 end")

    atom_of_ao, symbols = read_orca_atom_of_ao(lines, len(S))
    Z = [element_charge(s) for s in symbols]

    charge, multiplicity = get_orca_charge_multiplicity(lines)
    is_open_shell = multiplicity != 1

    if is_open_shell:
        SD = read_orca_matrix(
            lines, 'SPIN DENSITY',
            "SPIN DENSITY not found in the ORCA output. For open-shell systems add: "
            "%output Print[ P_Density ] 1 Print[ P_SpinDensity ] 1 end")

        total_electrons = sum(Z) - charge
        n_alpha_expected = (total_electrons + multiplicity - 1) / 2.0
        n_beta_expected = (total_electrons - multiplicity + 1) / 2.0
        trace_D = np.trace(D @ S)

        distance_to_alpha = abs(trace_D - n_alpha_expected)
        distance_to_beta = abs(trace_D - n_beta_expected)

        if min(distance_to_alpha, distance_to_beta) > 0.5:
            sys.exit("Tr(DENSITY . S) = " + str(round(trace_D, 4)) +
                      ", which matches neither the expected n_alpha (" + str(n_alpha_expected) +
                      ") nor n_beta (" + str(n_beta_expected) + "). Check CHARGE/MULTIPLICITY "
                      "or the printed matrices (ORCA).")

        if distance_to_alpha < distance_to_beta:
            P_alpha = D
            P_beta = D - SD
        else:
            P_beta = D
            P_alpha = D + SD

        P_list = [P_alpha, P_beta]
    else:
        P_list = [D]

    return {'S': S, 'P_list': P_list, 'atom_of_ao': atom_of_ao, 'Z': Z, 'symbols': symbols,
            'is_open_shell': is_open_shell, 'charge': charge, 'multiplicity': multiplicity}

def check_data(data, charge, program):
    n = len(data['S'])
    for P in data['P_list']:
        if P.shape != (n, n):
            sys.exit("S and P have different sizes (" + program + ").")

    electrons = sum(np.trace(P @ data['S']) for P in data['P_list'])
    expected = sum(data['Z']) - charge
    if abs(electrons - expected) > 0.05:
        sys.exit("Tr(PS) = " + str(round(electrons, 4)) + " electrons but " + str(expected) +
                 " were expected (" + program + "). Check CHARGE, the printed matrices, or ECPs.")

# ============================================================
# Mulliken / Lowdin: populations, charges and DI
# ============================================================

def atom_indicator(atom_of_ao, natoms):
    E = np.zeros((len(atom_of_ao), natoms))
    E[np.arange(len(atom_of_ao)), atom_of_ao] = 1.0
    return E

def symmetric_powers(S):
    eigenvalues, eigenvectors = np.linalg.eigh(S)
    S_half = eigenvectors @ np.diag(np.sqrt(eigenvalues)) @ eigenvectors.T
    S_inv_half = eigenvectors @ np.diag(1.0 / np.sqrt(eigenvalues)) @ eigenvectors.T
    return S_half, S_inv_half

def charges_and_DI(partition, data, E, S_half):
    # Mulliken: Mayer-type matrix PS. Lowdin: P in the S^(1/2) basis.
    S = data['S']

    matrices = []
    for P in data['P_list']:
        if partition == 'mulliken':
            matrices.append(P @ S)
        else:
            matrices.append(S_half @ P @ S_half)

    populations = E.T @ np.diag(sum(matrices))

    DI = np.zeros((E.shape[1], E.shape[1]))
    for M in matrices:
        if partition == 'mulliken':
            DI = DI + E.T @ (M * M.T) @ E
        else:
            DI = DI + E.T @ (M * M) @ E

    if data['is_open_shell']:
        DI = 2 * DI

    charges = np.array(data['Z']) - populations

    return charges, DI

# ============================================================
# AOMs in the occupied MO basis (as in ESIpy)
# ============================================================

def build_aoms(partition, data, S_half, S_inv_half):
    S = data['S']
    atom_of_ao = data['atom_of_ao']
    natoms = len(data['Z'])

    if data['is_open_shell']:
        densities = data['P_list']
    else:
        densities = [data['P_list'][0] / 2.0]      # per-spin density

    aoms_by_spin = []
    for D in densities:
        eigenvalues, eigenvectors = np.linalg.eigh(S_half @ D @ S_half)
        C_lowdin = eigenvectors[:, eigenvalues > 0.5]   # occupied MOs in the S^(1/2) basis
        C = S_inv_half @ C_lowdin

        aom = []
        for A in range(natoms):
            idxA = np.where(atom_of_ao == A)[0]
            if partition == 'mulliken':
                aom.append(C.T @ S[:, idxA] @ C[idxA, :])
            else:
                aom.append(C_lowdin[idxA, :].T @ C_lowdin[idxA, :])
        aoms_by_spin.append(aom)

    return aoms_by_spin

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
        starts = sorted(connec)    # ESIpy uses only atom 1, which may not belong to the graph
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

def nth_root(value, n):
    if value < 0:
        return -(abs(value) ** (1.0 / n))
    return value ** (1.0 / n)

def combine_channel_value(values, is_open_shell):
    if is_open_shell:
        return values[0] + values[1]
    return 2 * values[0]

def combine_channel_list(lists, is_open_shell):
    if is_open_shell:
        return list(np.add(lists[0], lists[1]))
    return [2 * v for v in lists[0]]

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

def write_info_txt(molecule_name, program, input_data, data, charge, multiplicity,
                    output_files, skipped, rings_by_partition):
    if data['is_open_shell']:
        system_type = 'open-shell'
        equations = 'UHF/UKS'
    else:
        system_type = 'closed-shell'
        equations = 'RHF/RKS'

    text = 'The ' + molecule_name + ' is ' + system_type + ', therefore, we use ' + equations + '.\n\n'
    text = text + 'Program: ' + program + '\n'

    if program == 'PYSCF':
        text = text + 'Functional/basis-set: ' + input_data['FUNCTIONAL/BASIS-SET'] + '\n'
    else:
        text = text + 'Output file read: ' + input_data['OUTPUT_FILE'] + '\n'

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

def analyse_partition(partition, input_data, data, E, S_half, S_inv_half, molecule_name,
                      want_di, want_iring, want_mci, want_pdi, want_flu, want_av1245,
                      output_files, rings_by_partition):
    symbols = data['symbols']
    is_open_shell = data['is_open_shell']
    molinfo = {'symbols': symbols}
    print('=== ' + partition.capitalize() + ' partition ===')

    charges, DI = charges_and_DI(partition, data, E, S_half)

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

    want_aom_based = want_iring or want_mci or want_pdi or want_flu or want_av1245
    if not want_aom_based:
        return

    aoms_by_spin = build_aoms(partition, data, S_half, S_inv_half)
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

    program = get_program(input_data)
    partitions = [PARTITION_NAMES[key] for key in PARTITION_NAMES if is_requested(input_data, key)]
    skipped = [key for key in NOT_IMPLEMENTED_PARTITIONS if is_requested(input_data, key)]

    if not partitions:
        print('Neither MULLIKEN nor LOWDIN is set to YES in INPUT.txt; nothing to do.')
        return
    if skipped:
        print('Not implemented in mulliken_lowdin.py, skipping: ' + ', '.join(skipped))

    molecule_name = input_data['NAME_MOLECULE']

    want_di = is_requested(input_data, 'DI')
    want_iring = is_requested(input_data, 'Iring')
    want_mci = is_requested(input_data, 'MCI')
    want_pdi = is_requested(input_data, 'PDI')
    want_flu = is_requested(input_data, 'FLU')
    want_av1245 = is_requested(input_data, 'AV1245')

    if program == 'PYSCF':
        if 'MOLECULAR_COORDINARES' not in input_data:
            sys.exit('MOLECULAR_COORDINARES is missing in INPUT.txt (required for PROGRAM = PySCF).')
        if 'FUNCTIONAL/BASIS-SET' not in input_data:
            sys.exit('FUNCTIONAL/BASIS-SET is missing in INPUT.txt (required for PROGRAM = PySCF).')
        if 'CHARGE' not in input_data:
            sys.exit('CHARGE is missing in INPUT.txt (required for PROGRAM = PySCF).')
        if 'MULTIPLICITY' not in input_data:
            sys.exit('MULTIPLICITY is missing in INPUT.txt (required for PROGRAM = PySCF).')

        charge = int(input_data['CHARGE'])
        multiplicity = int(input_data['MULTIPLICITY'])
        symbols, coordinates = read_xyz(input_data['MOLECULAR_COORDINARES'])
        data = get_pyscf_data(input_data, symbols, coordinates, charge, multiplicity)
    else:
        if 'OUTPUT_FILE' not in input_data:
            sys.exit('OUTPUT_FILE is missing in INPUT.txt (required for PROGRAM = GAUSSIAN or ORCA).')

        if program == 'GAUSSIAN':
            data = get_gaussian_data(input_data['OUTPUT_FILE'])
        else:
            data = get_orca_data(input_data['OUTPUT_FILE'])

    charge = data['charge']
    multiplicity = data['multiplicity']

    check_data(data, charge, program)

    natoms = len(data['Z'])
    E = atom_indicator(data['atom_of_ao'], natoms)
    S_half, S_inv_half = symmetric_powers(data['S'])

    output_files = []
    rings_by_partition = {}
    for partition in partitions:
        analyse_partition(partition, input_data, data, E, S_half, S_inv_half, molecule_name,
                          want_di, want_iring, want_mci, want_pdi, want_flu, want_av1245,
                          output_files, rings_by_partition)

    write_info_txt(molecule_name, program, input_data, data, charge, multiplicity,
                    output_files, skipped, rings_by_partition)


if __name__ == '__main__':
    main()




