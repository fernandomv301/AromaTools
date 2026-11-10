# Methods and references

AromaTools brings together methodologies that were developed by different research groups. This page lists each method implemented in the program together with the work in which it was proposed. **If you use a method through AromaTools, please cite AromaTools and the works of the methods you used.**

Numbers in brackets refer to the [reference list](#reference-list) at the end of this page.

---

## AroGeometric

| Method | Description | References |
|--------|-------------|------------|
| HOMA93 | Harmonic oscillator model of aromaticity, based on bond-length equalization | [1,2,3] |
| HOMAc | HOMA reparameterized with aromatic **and antiaromatic** | [4] |
| HOMER | HOMA reparameterized for the excited state (T₁) | [5] |
| BLA | Bond length alternation: difference between the lengths of consecutive bonds of the ring | [6] |

---

## AroMagnetic

| Method | Description | References |
|--------|-------------|------------|
| NICS | Nucleus-independent chemical shift, calculated with ghost atoms | [7,8] |
| NICS-scan | NICS profiles along the axis perpendicular to the ring | [9] |
| INICS | Integral of the NICS profile along the normal to the ring, independent of the height of the probe | [10,11] |
| NICS-XY-scan | NICS probes placed along the X (and Y) axis, to identify local and global ring currents in multi-ring systems | [12] |
| Predi-XY | Automated generation of NICS-XY-scans for polycyclic aromatic systems using an additivity scheme | [13] |
| Induced magnetic field (B<sup>ind</sup>) | Graphical representation of the induced magnetic field under an external field | [14] |
| σ and π contributions to B<sup>ind</sup> | Separation of the σ and π contributions to the induced magnetic field | [15] |
| Ring current strength | Integration of the shielding tensor along a profile of dummy atoms, using the Ampère–Maxwell law | [16] |
| Pseudo-π model | π-electron magnetic response from a hydrogen-atom model of the carbon framework | [17] |

!!! note
    The induced magnetic field is evaluated from the shielding tensor [18] obtained with Gaussian [19] or ORCA [20,21,22], which are the programs that perform the magnetic calculations.

---

## AroElectronics

AroElectronics **does not reimplement the electronic indices**. It imports from [ESIpy](https://github.com/jgrebol/ESIpy) the functions needed to build the atomic overlap matrices and to compute the indices (`ESI`, the wavefunction reader and the reader of AIMAll files), and adds the AromaTools workflow around them: the `INPUT.txt` file, the detection of rings, the selection of indices and partitions, and the output files. The values obtained with AroElectronics are therefore those of ESIpy [23], and the works in which the indices and partitions were introduced must be cited.

**Indices**

| Method | Description | References |
|--------|-------------|------------|
| Delocalization index (DI) | Number of electron pairs shared between two atoms; the basis of the electronic indices | [24] |
| Iring | Multicenter bond index of a ring | [25] |
| MCI | Multicenter index (sum of Iring over the permutations of the ring atoms) | [26] |
| PDI | Para-delocalization index (six-membered rings) | [27] |
| FLU | Aromatic fluctuation index | [28] |
| AV1245 and AVmin | Average and minimum of the four-center indices of atoms 1-2-4-5 of the ring, for large rings | [29] |
| BOA | Bond-order alternation, the electronic counterpart of BLA (uses the DI instead of the bond length). | [30] |

**Atomic partitions**

| `PARTITION` | References |
|-------------|------------|
| `mulliken` | [31] |
| `lowdin` | [32] |
| `meta-lowdin` | [33] |
| `nao` | [34] |
| `iao` | [35] |
| `qtaim` | Bader's Quantum Theory of Atoms in Molecules [36]; the AOMs are obtained with AIMAll [37] |

PySCF [38] is used to run the calculations from scratch, and wavefunction files are read with IOData [39].

---

## Reference list

1. Kruszewski, J.; Krygowski, T. M. Definition of aromaticity basing on the harmonic oscillator model. *Tetrahedron Lett.* **1972**, *13*, 3839–3842.
2. Krygowski, T. M. Crystallographic studies of inter- and intramolecular interactions reflected in aromatic character of π-electron systems. *J. Chem. Inf. Comput. Sci.* **1993**, *33*, 70–78.
3. Krygowski, T. M.; Cyrański, M. K. Structural aspects of aromaticity. *Chem. Rev.* **2001**, *101*, 1385–1419.
4. Arpa, E. M.; Stafström, S.; Durbeej, B. HOMAc: A Parameterization of the Harmonic Oscillator Model of Aromaticity (HOMA) That Includes Antiaromaticity. *J. Org. Chem.* **2025**, *90*, 1297–1308.
5. Arpa, E. M.; Durbeej, B. HOMER: a reparameterization of the harmonic oscillator model of aromaticity (HOMA) for excited states. *Phys. Chem. Chem. Phys.* **2023**, *25*, 16763.
6. Casademont-Reig, I.; Ramos-Cordoba, E.; Torrent-Sucarrat, M.; Matito, E. Aromaticity descriptors based on electron delocalization. In *Aromaticity: Modern Computational Methods and Applications*; Elsevier, 2021; Chapter 7. DOI: 10.1016/C2019-0-04193-3.
7. Schleyer, P. v. R.; Maerker, C.; Dransfeld, A.; Jiao, H.; van Eikema Hommes, N. J. R. Nucleus-Independent Chemical Shifts: A Simple and Efficient Aromaticity Probe. *J. Am. Chem. Soc.* **1996**, *118*, 6317–6318.
8. Chen, Z.; Wannere, C. S.; Corminboeuf, C.; Puchta, R.; Schleyer, P. v. R. Nucleus-Independent Chemical Shifts (NICS) as an Aromaticity Criterion. *Chem. Rev.* **2005**, *105*, 3842–3888.
9. Stanger, A. Nucleus-Independent Chemical Shifts (NICS): Distance Dependence and Revised Criteria for Aromaticity and Antiaromaticity. *J. Org. Chem.* **2006**, *71*, 883–893.
10. Stanger, A. Reexamination of NICSπ,zz: Height Dependence, Off-Center Values, and Integration. *J. Phys. Chem. A* **2019**, *123*, 3922–3927.
11. Dudek, W. M.; Ostrowski, S.; Dobrowolski, J. Cz. On Aromaticity of the Aromatic α-Amino Acids and Tuning of the NICS Indices to Find the Aromaticity Order. *J. Phys. Chem. A* **2022**, *126*, 3433–3444.
12. Gershoni-Poranne, R.; Stanger, A. The NICS-XY-Scan: Identification of Local and Global Ring Currents in Multi-Ring Systems. *Chem. Eur. J.* **2014**, *20*, 5673–5688. DOI: 10.1002/chem.201304307.
13. Wahab, A.; Fleckenstein, F.; Feusi, S.; Gershoni-Poranne, R. Predi-XY: a python program for automated generation of NICS-XY-scans based on an additivity scheme. *Electron. Struct.* **2020**, *2* (4). DOI: 10.1088/2516-1075/abd081.
14. Merino, G.; Heine, T.; Seifert, G. The Induced Magnetic Field in Cyclic Molecules. *Chem. Eur. J.* **2004**, *10*, 4367–4371.
15. Heine, T.; Islas, R.; Merino, G. σ and π Contributions to the Induced Magnetic Field: Indicators for the Mobility of Electrons in Molecules. *J. Comput. Chem.* **2007**, *28*, 302–309.
16. Berger, R. J. F.; Dimitrova, M.; Nasibullin, R. T.; Valiev, R. R.; Sundholm, D. Integration of global ring currents using the Ampère–Maxwell law. *Phys. Chem. Chem. Phys.* **2022**, *24*, 624–628.
17. Fowler, P. W.; Steiner, E. Pseudo-π currents: rapid and accurate visualisation of ring currents in conjugated hydrocarbons. *Chem. Phys. Lett.* **2002**, *364*, 259–266.
18. Anet, F. A. L.; O'Leary, D. J. The Shielding Tensor. Part I: Understanding its Symmetry Properties. *Concepts Magn. Reson.* **1991**, *3*, 193–214.
19. Frisch, M. J.; Trucks, G. W.; Schlegel, H. B.; *et al.* *Gaussian 16*, Revision C.01; Gaussian, Inc.: Wallingford, CT, 2016.
20. Neese, F. The ORCA program system. *WIREs Comput. Mol. Sci.* **2012**, *2*, 73–78.
21. Neese, F. Software update: The ORCA program system—Version 5.0. *WIREs Comput. Mol. Sci.* **2022**, *12*, e1606.
22. Neese, F. Software update: The ORCA program system—Version 6. *WIREs Comput. Mol. Sci.* **2025**, *15*, e70019.
23. Grèbol-Tomàs, J.; Matito, E.; Salvador, P. Can Aromaticity Be Evaluated Using Atomic Partitions Based on the Hilbert-Space? *Chem. Eur. J.* **2024**, *30*, e202401282.
24. Fradera, X.; Austen, M. A.; Bader, R. F. W. The Lewis Model and Beyond. *J. Phys. Chem. A* **1999**, *103*, 304–314.
25. Giambiagi, M.; de Giambiagi, M. S.; dos Santos Silva, C. D.; de Figueiredo, A. P. Multicenter bond indices as a measure of aromaticity. *Phys. Chem. Chem. Phys.* **2000**, *2*, 3381–3392.
26. Bultinck, P.; Ponec, R.; Van Damme, S. Multicenter bond indices as a new measure of aromaticity in polycyclic aromatic hydrocarbons. *J. Phys. Org. Chem.* **2005**, *18*, 706–718.
27. Poater, J.; Fradera, X.; Duran, M.; Solà, M. The Delocalization Index as an Electronic Aromaticity Criterion: Application to a Series of Planar Polycyclic Aromatic Hydrocarbons. *Chem. Eur. J.* **2003**, *9*, 400–406.
28. Matito, E.; Duran, M.; Solà, M. The aromatic fluctuation index (FLU): A new aromaticity index based on electron delocalization. *J. Chem. Phys.* **2005**, *122*, 014109.
29. Matito, E. An electronic aromaticity index for large rings. *Phys. Chem. Chem. Phys.* **2016**, *18*, 11839–11846.
30. Casademont-Reig, I.; Ramos-Cordoba, E.; Torrent-Sucarrat, M.; Matito, E. How do the Hückel and Baird Rules Fade Away in Annulenes? *Molecules* **2020**, *25*, 711.
31. Mulliken, R. S. Electronic Population Analysis on LCAO–MO Molecular Wave Functions. I. *J. Chem. Phys.* **1955**, *23*, 1833–1840.
32. Löwdin, P.-O. On the Non-Orthogonality Problem Connected with the Use of Atomic Wave Functions in the Theory of Molecules and Crystals. *J. Chem. Phys.* **1950**, *18*, 365–375.
33. Sun, Q.; Chan, G. K.-L. Exact and Optimal Quantum Mechanics/Molecular Mechanics Boundaries. *J. Chem. Theory Comput.* **2014**, *10*, 3784–3790.
34. Reed, A. E.; Weinstock, R. B.; Weinhold, F. Natural population analysis. *J. Chem. Phys.* **1985**, *83*, 735–746.
35. Knizia, G. Intrinsic Atomic Orbitals: An Unbiased Bridge between Quantum Theory and Chemical Concepts. *J. Chem. Theory Comput.* **2013**, *9*, 4834–4843.
36. Bader, R. F. W. *Atoms in Molecules: A Quantum Theory*; Oxford University Press: Oxford, 1990.
37. Keith, T. A. *AIMAll*, TK Gristmill Software: Overland Park, KS, USA. aim.tkgristmill.com.
38. Sun, Q.; Berkelbach, T. C.; Blunt, N. S.; *et al.* PySCF: the Python-based simulations of chemistry framework. *WIREs Comput. Mol. Sci.* **2018**, *8*, e1340.
39. Verstraelen, T.; *et al.* IOData: A python library for reading, writing, and converting computational chemistry file formats and generating input files. *J. Comput. Chem.* **2021**, *42*, 458–464.
