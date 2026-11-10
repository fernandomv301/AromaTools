import os
import re
import sys
import numpy as np

# Insert information!
ATOMICFILES_PATH = 'benzene (s0)/AIM-benzene/benzene_atomicfiles'
# ATOMICFILES_PATH = 'benzene (t1)/AIM-benzene_tr-camb3lyp-AQ/benzene_tr-camb3lyp-AQ_atomicfiles'

def read_lines(path):
    file = open(path, 'r')
    lines = file.readlines()
    file.close()
    return lines

def list_int_files(path):
    names = []
    for name in os.listdir(path):
        if name.endswith('.int'): # .int files!
            names.append(name)
    return names

def read_atom_label(lines, file_name):
    for line in lines:
        match = re.search(r'Integration is over atom\s+([A-Za-z]+)(\d+)', line)
        if match:
            return match.group(1).capitalize(), int(match.group(2))

def read_aom(lines, file_name): # AOM: lower triangle, 6 values per line, over all the orbitals (occupied and virtual)
    start = None
    for i in range(len(lines)):
        if 'The Atomic Overlap Matrix' in lines[i]:
            start = i
            break
    wavefunction_line = ''
    i = start + 1
    while i < len(lines) and not re.match(r'^\s*-?\d', lines[i]):
        if 'wavefunction' in lines[i]:
            wavefunction_line = lines[i].strip()
        i = i + 1

    values = []
    while i < len(lines) and re.match(r'^\s*-?\d', lines[i]):
        for word in lines[i].split():
            values.append(float(word))
        i = i + 1

    n = int(round((-1 + np.sqrt(1 + 8 * len(values))) / 2))
    if n * (n + 1) // 2 != len(values):
        sys.exit(file_name + ": " + str(len(values)) + " values do not form a triangular matrix.")

    lower = np.zeros((n, n))
    lower[np.tril_indices(n)] = values
    matrix = lower + lower.T - np.diag(np.diag(lower))

    return matrix, wavefunction_line

def find_wfx(lines, path): # The .wfx can be in the atomic files folder or in a folder above it
    wfx_name = None
    for line in lines:
        if 'Wfx File:' in line:
            wfx_name = os.path.basename(line.split('Wfx File:')[1].strip())
            break
   
    folder = os.path.abspath(path)
    for level in range(4):
        candidate = os.path.join(folder, wfx_name)
        if os.path.isfile(candidate):
            return candidate
        folder = os.path.dirname(folder)

def read_wfx_tag(text, tag):
    match = re.search('<' + re.escape(tag) + '>\n(.*?)\n</' + re.escape(tag) + '>', text, re.S)
    return match.group(1)

def read_wfx(wfx_path):
    file = open(wfx_path, 'r')
    text = file.read()
    file.close()

    occupations = [float(word) for word in read_wfx_tag(text, 'Molecular Orbital Occupation Numbers').split()]
    spin_types = [line.strip() for line in read_wfx_tag(text, 'Molecular Orbital Spin Types').split('\n')]

    return {'occupations': occupations, 'spin_types': spin_types,
            'n_nuclei': int(float(read_wfx_tag(text, 'Number of Nuclei'))),
            'n_electrons': int(float(read_wfx_tag(text, 'Number of Electrons'))),
            'charge': int(float(read_wfx_tag(text, 'Net Charge'))),
            'multiplicity': int(float(read_wfx_tag(text, 'Electronic Spin Multiplicity')))}

def occupied_orbitals(wfx, is_open_shell): # Indices of the occupied orbitals of each spin
    occupations = wfx['occupations']
    spin_types = wfx['spin_types']

    if is_open_shell:
        sets = [[], []]
        expected = 1.0
    else:
        sets = [[]]
        expected = 2.0

    for i in range(len(occupations)):
        if occupations[i] < 1e-6:
            continue
        if abs(occupations[i] - expected) > 1e-3:
            sys.exit("Orbital " + str(i + 1) + " has occupation " + str(occupations[i]) + " (expected " + str(expected) +
                     "). Fractional occupations (natural orbitals) are not supported.")
        if not is_open_shell:
            sets[0].append(i)
        elif spin_types[i] == 'Alpha':
            sets[0].append(i)
        elif spin_types[i] == 'Beta':
            sets[1].append(i)

    return sets

def read_qtaim(path):
    file_names = list_int_files(path)

    atoms = {} # atom number -> (symbol, AOM)
    wavefunction_lines = []
    first_lines = None
    for file_name in file_names:
        lines = read_lines(os.path.join(path, file_name))
        symbol, number = read_atom_label(lines, file_name)
        matrix, wavefunction_line = read_aom(lines, file_name)

        atoms[number] = (symbol, matrix)
        wavefunction_lines.append(wavefunction_line)
        if first_lines is None:
            first_lines = lines

    natoms = len(atoms)
    missing = [str(number) for number in range(1, natoms + 1) if number not in atoms]

    n = atoms[1][1].shape[0]
    for number in atoms:
        is_open_shell = 'Unrestricted' in wavefunction_lines[0]
        wfx = read_wfx(find_wfx(first_lines, path))

    symbols = []
    for number in range(1, natoms + 1):
        symbols.append(atoms[number][0])

    occupied = occupied_orbitals(wfx, is_open_shell)

    aoms_by_spin = []
    for indices in occupied:
        aom = []
        for number in range(1, natoms + 1):
            aom.append(atoms[number][1][np.ix_(indices, indices)])
        aoms_by_spin.append(aom)

    return {'aoms_by_spin': aoms_by_spin, 'symbols': symbols, 'is_open_shell': is_open_shell,
            'charge': wfx['charge'], 'multiplicity': wfx['multiplicity'], 'n_electrons': wfx['n_electrons'],
            'occupied': occupied, 'n_orbitals': n}

def check_aoms(data): # Integration quality: the AOMs of all the atoms must add up to the identity
    aoms_by_spin = data['aoms_by_spin']

    electrons = 0.0
    for spin in range(len(aoms_by_spin)):
        aom = aoms_by_spin[spin]
        deviation = np.abs(sum(aom) - np.identity(aom[0].shape[0])).max()
        print('  spin ' + str(spin + 1) + ': ' + str(aom[0].shape[0]) + ' occupied orbitals, max|sum of AOMs - I| = ' + '{:.1e}'.format(deviation))
        for atom_aom in aom:
            electrons = electrons + np.trace(atom_aom)

    if not data['is_open_shell']:
        electrons = 2.0 * electrons
    print('  electrons from the AOMs = ' + str(round(electrons, 4)) + ' (wfx: ' + str(data['n_electrons']) + ')')

def format_orbitals(indices): # Orbital numbers (starting at 1) as ranges, for example 1-21
    numbers = [i + 1 for i in indices]
    ranges = []
    start = numbers[0]
    for k in range(1, len(numbers) + 1):
        if k == len(numbers) or numbers[k] != numbers[k - 1] + 1:
            if start == numbers[k - 1]:
                ranges.append(str(start))
            else:
                ranges.append(str(start) + '-' + str(numbers[k - 1]))
            if k < len(numbers):
                start = numbers[k]
    return ', '.join(ranges)

def main():
    data = read_qtaim(ATOMICFILES_PATH)
    aoms_by_spin = data['aoms_by_spin']
    symbols = data['symbols']

    if data['is_open_shell']:
        wavefunction = 'unrestricted'
        spin_names = ['alpha', 'beta']
    else:
        wavefunction = 'restricted'
        spin_names = ['alpha = beta']

    print('Atomic files: ' + ATOMICFILES_PATH)
    print('Wavefunction: ' + wavefunction)
    print('Charge: ' + str(data['charge']) + '   Multiplicity: ' + str(data['multiplicity']) + '   Electrons: ' + str(data['n_electrons']))
    print()

    print('Atoms (atom number = position in aom[spin] + 1), with N(A) = electrons in the atom from the AOM:')
    for A in range(len(symbols)):
        population = 0.0
        for aom in aoms_by_spin:
            population = population + np.trace(aom[A])
        if not data['is_open_shell']:
            population = 2.0 * population
        print('  ' + str(A + 1) + '  ' + symbols[A] + '   N = ' + str(round(population, 6)))
    print()

    print('Orbitals: ' + str(data['n_orbitals']) + ' in each AOM (occupied and virtual), only the occupied ones are kept')
    for spin in range(len(aoms_by_spin)):
        print('  ' + spin_names[spin] + ': ' + str(len(data['occupied'][spin])) + ' occupied, orbitals ' + format_orbitals(data['occupied'][spin]))
    print()

    shape = aoms_by_spin[0][0].shape
    print('What the aromaticity indices need:')
    print("  data['aoms_by_spin'][spin][atom - 1]: " + str(len(aoms_by_spin)) + ' spin set(s) x ' + str(len(symbols)) + ' atoms, each AOM is ' + str(shape[0]) + 'x' + str(shape[1]))
    print("  data['symbols'], data['is_open_shell']")
    print()

    print('Integration quality:')
    check_aoms(data)


if __name__ == '__main__':
    main()
