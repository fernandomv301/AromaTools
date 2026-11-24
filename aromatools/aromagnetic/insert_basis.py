#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob

BASIS_FILE = "basis.txt"

with open(BASIS_FILE, "r", encoding="utf-8") as f:
    basis_content = f.read()

for out_file in glob.glob("*.inp"):
    with open(out_file, "r", encoding="utf-8") as f:
        text = f.read()

    # With POP=NBO6read the $NBO keylist must be the last section of the Gaussian input,
    # so the basis set goes before it.
    pos = text.find("\n$NBO")
    if pos != -1:
        pos += 1
        new_text = (text[:pos] + "\n" + basis_content.rstrip("\n") + "\n\n" + text[pos:])
        with open(out_file, "w", encoding="utf-8") as f:
            f.write(new_text)
    else:
        with open(out_file, "a", encoding="utf-8") as f:
            f.write("\n")          # separador seguro
            f.write(basis_content)

print("Basis inserted!!")
