import sys
import numpy as np
from itertools import permutations
from pyscf import gto, dft, lo
import warnings
warnings.filterwarnings('ignore', category=RuntimeWarning, message='.*matmul.*')
warnings.filterwarnings('ignore', category=UserWarning, message='.*B3LYP.*')

# Insert information!
MOLECULE_NAME = 'benzene'
XYZ_PATH = 'molecules/benzene.xyz'
FUNCTIONAL = 'b3lyp'
BASIS = '6-311G**'
# mol.basis = '6-311G**'
CHARGE = 0
MULTIPLICITY = 1
RING = [1,2,6,10,7,3] # Benzene
# RING = [1,2,3,4,5,6] # Cyclohexane 
# RING = [10,5,6,7,8,9] # Naphtalene
# RING = [1,2,3,4,12,11] # Anthracene R1
# RING = [10,12,11,9,14,13] # Anthracene R2
# RING = [4,3,7,11,12,8] # Phenantrene R1
# RING = [1,2,5,9,8,4] # Phenantrene R2

PARTITIONS = ['mulliken', 'lowdin', 'meta-lowdin', 'nao', 'iao']

FLU_REFERENCES = {
    'mulliken':    {'CC': 1.4536, 'CN': 1.4051, 'BN': 1.1042, 'NN': 1.2831, 'CS': 1.1013, 'CO': 0.9840},
    'lowdin':      {'CC': 1.4910, 'CN': 1.5864, 'BN': 1.3165, 'NN': 1.5277, 'CS': 1.2213, 'CO': 1.2368},
    'meta-lowdin': {'CC': 1.4392, 'CN': 1.4521, 'BN': 1.2048, 'NN': 1.3700, 'CS': 1.1455, 'CO': 1.0802},
    'nao':         {'CC': 1.4347, 'CN': 1.4115, 'BN': 0.9239, 'NN': 1.3708, 'CS': 1.1625, 'CO': 0.9912},
    'iao':         {'CC': 1.4378, 'CN': 1.4386, 'BN': 1.1631, 'NN': 1.3606, 'CS': 1.1433, 'CO': 1.0708},
}

def read_xyz(xyz_path):
    file = open(xyz_path, 'r')
    lines = file.readlines()
    file.close()
    number_of_atoms = int(lines[0].strip())

    symbols = []
    coordinates = []

    for i in range(2, 2 + number_of_atoms):
        parts = lines[i].split()
        symbol = parts[0]
        x = float(parts[1])
        y = float(parts[2])
        z = float(parts[3])

        symbols.append(symbol)
        coordinates.append([x, y, z])

    return symbols, coordinates

def build_atom_text(symbols, coordinates): # Build pyscf input
    atom_text = ''
    for i in range(len(symbols)):
        atom_text = atom_text + symbols[i] + ' '
        atom_text = atom_text + str(coordinates[i][0]) + ' '
        atom_text = atom_text + str(coordinates[i][1]) + ' '
        atom_text = atom_text + str(coordinates[i][2]) + '\n'
    return atom_text

def ao_atom_assignment(mol):
    ao_labels = mol.ao_labels(fmt=None, base=0)
    atom_of_basis_list = []
    for label in ao_labels:
        atom_of_basis_list.append(label[0])
    return np.array(atom_of_basis_list)

# ---------------- PySCF ----------------

def run_scf():
    symbols, coordinates = read_xyz(XYZ_PATH)
    atom_text = build_atom_text(symbols, coordinates)

    mol = gto.Mole()
    mol.atom = atom_text
    mol.basis = BASIS
    mol.charge = CHARGE
    mol.spin = MULTIPLICITY - 1
    mol.build()

    is_open_shell = (MULTIPLICITY != 1)
    if is_open_shell:
        mf = dft.UKS(mol).density_fit()
    else:
        mf = dft.RKS(mol).density_fit()
    mf.xc = FUNCTIONAL
    mf.kernel()

    return mol, mf, is_open_shell

# ======================= BEGIN =======================

def get_orthogonalizer(partition, mol, mf, S): # Make orthogonalization! Lowdin, meta-Lowdin and NAO
    if partition == 'lowdin':
        U = lo.orth.lowdin(S)
    elif partition == 'meta-lowdin':
        U = lo.orth_ao(mf, 'meta_lowdin')
    elif partition == 'nao':
        U = lo.orth_ao(mf, 'nao')
    atom_of_basis = ao_atom_assignment(mol)
    return U, atom_of_basis
 
def get_iao_transform(mol, S, mo_occ_one_spin): # IAO partition
    a = lo.iao.iao(mol, mo_occ_one_spin)
    a = lo.vec_lowdin(a, S)
    atom_of_basis = ao_atom_assignment(lo.iao.reference_mol(mol))
    return a, atom_of_basis

def build_aoms(partition, mol, mf, S, is_open_shell): # Build AO matrix
    natoms = mol.natm

    if is_open_shell:
        mo_occ_sets = [mf.mo_coeff[0][:, mf.mo_occ[0] > 0], mf.mo_coeff[1][:, mf.mo_occ[1] > 0]]
    else:
        mo_occ_sets = [mf.mo_coeff[:, mf.mo_occ > 0]]

    aoms_by_spin = []
    for mo_occ in mo_occ_sets:
        aom = []

        if partition == 'mulliken': # No basis: the AOs of atom A are used as they are
            atom_of_basis = ao_atom_assignment(mol)
            for A in range(natoms):
                idxA = np.where(atom_of_basis == A)[0]
                aom.append(mo_occ.T @ S[:, idxA] @ mo_occ[idxA, :])
        else:
            if partition == 'iao':
                U, atom_of_basis = get_iao_transform(mol, S, mo_occ)
            else:
                U, atom_of_basis = get_orthogonalizer(partition, mol, mf, S)

            C_new = U.T @ S @ mo_occ
            for A in range(natoms):
                idxA = np.where(atom_of_basis == A)[0]
                aom.append(C_new[idxA, :].T @ C_new[idxA, :])

        aoms_by_spin.append(aom)

    return aoms_by_spin

# ======================= END =======================

def compute_iring_os(ring, aom_one_spin):
    product = np.identity(aom_one_spin[0].shape[0])
    for atom in ring:
        product = product @ aom_one_spin[atom - 1]
    return (2 ** (len(ring) - 1)) * np.trace(product)

def compute_iring(ring, aom, is_open_shell):
    irings = [compute_iring_os(ring, spin_aom) for spin_aom in aom]
    if is_open_shell:
        return irings[0] + irings[1]
    return 2 * irings[0]

def compute_di(aom, A, B):
    return compute_iring([A, B], aom, is_open_shell=(len(aom) == 2))

def compute_mci_os(ring, aom_one_spin):
    mats = [aom_one_spin[atom - 1] for atom in ring]
    n = len(mats)
    if n < 3:
        return 0.0

    total = 0.0
    for orden in permutations(range(1, n)):
        producto = mats[0]
        for i in orden:
            producto = producto @ mats[i]
        total = total + np.trace(producto)

    return (2 ** (n - 2)) * total

def compute_mci(ring, aom, is_open_shell):
    mcis = [compute_mci_os(ring, spin_aom) for spin_aom in aom]
    if is_open_shell:
        return mcis[0] + mcis[1]
    return 2 * mcis[0]

def nth_root(value, n): # Normalization!
    if value < 0:
        return -(abs(value) ** (1.0 / n))
    return value ** (1.0 / n)

def compute_pdi(ring, aom):
    if len(ring) != 6:
        return None
    pdi_a = compute_di(aom, ring[0], ring[3])
    pdi_b = compute_di(aom, ring[1], ring[4])
    pdi_c = compute_di(aom, ring[2], ring[5])
    return (pdi_a + pdi_b + pdi_c) / 3.0

def av1245_quadruplets(ring):
    n = len(ring)
    quadruplets = []
    for i in range(n):
        quadruplets.append((ring[i % n], ring[(i + 1) % n], ring[(i + 3) % n], ring[(i + 4) % n]))
    return quadruplets

def compute_av1245(ring, aom, is_open_shell):
    if len(ring) < 6:
        return None, None, []

    terms = []
    for quadruplet in av1245_quadruplets(ring):
        mci_4center = compute_mci(list(quadruplet), aom, is_open_shell)
        terms.append(1000.0 * mci_4center / 3.0)

    av1245_value = sum(terms) / len(terms)
    avmin_value = min(terms, key=abs)
    return av1245_value, avmin_value, terms

def atom_pop(aom, A, is_open_shell): # Electrons in atom A
    traces = [np.trace(spin_aom[A - 1]) for spin_aom in aom]
    if is_open_shell:
        return traces[0] + traces[1]
    return 2.0 * traces[0]

def atom_loc(aom, A, is_open_shell): # Electrons localized in atom A (localization index)
    traces = [np.trace(spin_aom[A - 1] @ spin_aom[A - 1]) for spin_aom in aom]
    if is_open_shell:
        return traces[0] + traces[1]
    return 2.0 * traces[0]

def compute_flu_os(ring, aom_one_spin, symbols, flu_refs):
    n = len(ring)

    bond_types = []
    for i in range(n):
        symbol_a = symbols[ring[i] - 1].capitalize()
        symbol_b = symbols[ring[(i + 1) % n] - 1].capitalize()
        bond_types.append(''.join(sorted([symbol_a, symbol_b])))

    dis = []
    ns = []
    lis = []
    for i in range(n):
        A = aom_one_spin[ring[i] - 1]
        B = aom_one_spin[ring[(i + 1) % n] - 1]
        dis.append(4.0 * np.trace(A @ B))
        ns.append(2.0 * np.trace(A))
        lis.append(2.0 * np.trace(A @ A))

    flu_value = 0.0
    for i in range(n):
        reference = flu_refs[bond_types[i]]
        flu_deloc = (dis[i] - reference) / reference

        a_to_b = dis[i] / 2.0 * (ns[i] - lis[i])
        b_to_a = dis[i] / 2.0 * (ns[(i + 1) % n] - lis[(i + 1) % n])
        flu_polar = a_to_b / b_to_a
        if flu_polar < 1.0:
            flu_polar = 1.0 / flu_polar

        flu_value = flu_value + (flu_deloc * flu_polar) ** 2

    return flu_value / n

def compute_flu(ring, aom, symbols, flu_refs, is_open_shell):
    values = [compute_flu_os(ring, spin_aom, symbols, flu_refs) for spin_aom in aom]
    if None in values:
        return None
    if is_open_shell:
        return values[0] + values[1]
    return values[0]

def compute_boa_os(ring, aom_one_spin):
    n = len(ring)

    dis = [] # DI of the bond i-(i+1) in one spin
    for i in range(n):
        A = aom_one_spin[ring[i] - 1]
        B = aom_one_spin[ring[(i + 1) % n] - 1]
        dis.append(2.0 * np.trace(A @ B))

    even_bonds = dis[0::2]
    odd_bonds = dis[1::2]
    boa = abs(sum(even_bonds) / len(even_bonds) - sum(odd_bonds) / len(odd_bonds))

    boa_c = 0.0 # Difference between neighbouring bonds, around the ring
    for i in range(n):
        boa_c = boa_c + abs(dis[i] - dis[(i + 1) % n]) / n

    return boa, boa_c

def compute_boa(ring, aom, is_open_shell):
    if len(ring) % 2 != 0:
        return None, None

    values = [compute_boa_os(ring, spin_aom) for spin_aom in aom]
    if is_open_shell:
        return values[0][0] + values[1][0], values[0][1] + values[1][1]
    return 2 * values[0][0], 2 * values[0][1]

# MAIN PROGRAM

def main():
    mol, mf, is_open_shell = run_scf()
    S = mol.intor('int1e_ovlp')
    symbols = [mol.atom_symbol(A) for A in range(mol.natm)]
    natoms = mol.natm

    print('Molecule: ' + MOLECULE_NAME)
    print('Functional/basis: ' + FUNCTIONAL + '/' + BASIS)
    print('Charge: ' + str(CHARGE) + '   Multiplicity: ' + str(MULTIPLICITY))
    print('Ring: ' + '-'.join(str(atom) for atom in RING))
    print()

    for partition in PARTITIONS:
        print('=== ' + partition + ' partition ===')

        aom = build_aoms(partition, mol, mf, S, is_open_shell)

        print('Atomic populations (ring atoms):')
        for A in RING:
            population = atom_pop(aom, A, is_open_shell)
            localization = atom_loc(aom, A, is_open_shell)
            charge = mol.atom_charge(A - 1) - population
            print(symbols[A - 1] + str(A) + ' -> delocalized = ' + str(round(population - localization, 6)))
        print("\n")
        print('DI (ring bonds):')
        for i in range(len(RING)):
            A = RING[i]
            B = RING[(i + 1) % len(RING)]
            label = symbols[A - 1] + str(A) + '-' + symbols[B - 1] + str(B)
            print(label + ' = ' + str(round(compute_di(aom, A, B), 6)))

        print("\n")
        print(" === Aromaticity index === ")
        print("\n")
        
        iring_value = compute_iring(RING, aom, is_open_shell)
        print('Iring        = ' + str(round(iring_value, 6)))
        print('Iring^(1/n)  = ' + str(round(nth_root(iring_value, len(RING)), 6)))
        
        mci_value = compute_mci(RING, aom, is_open_shell) # Just 14 rings with atoms
        print('MCI          = ' + str(round(mci_value, 6)))
        print('MCI^(1/n)    = ' + str(round(nth_root(mci_value, len(RING)), 6)))

        pdi_value = compute_pdi(RING, aom)
        if pdi_value is None:
            print('PDI          = not defined (ring size != 6)')
        else:
            print('PDI          = ' + str(round(pdi_value, 6)))

        flu_value = compute_flu(RING, aom, symbols, FLU_REFERENCES[partition], is_open_shell)
        print('FLU          = ' + str(round(flu_value, 6)))

        av1245_value, avmin_value, av_terms = compute_av1245(RING, aom, is_open_shell)
        if av1245_value is None:
            print('AV1245/AVmin = not defined (ring size < 6)')
        else:
            print('AV1245       = ' + str(round(av1245_value, 6)))
            print('AVmin        = ' + str(round(avmin_value, 6)))

        boa_value, boac_value = compute_boa(RING, aom, is_open_shell)
        if boa_value is None:
            print('BOA          = not defined (odd ring size)')
        else:
            print('BOA          = ' + str(round(boa_value, 6)))
            print('BOA_c        = ' + str(round(boac_value, 6)))
        print()


if __name__ == '__main__':
    main()
