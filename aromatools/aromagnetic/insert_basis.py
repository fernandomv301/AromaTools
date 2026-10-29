#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import glob

BASIS_FILE = "basis.txt"

with open(BASIS_FILE, "r", encoding="utf-8") as f:
    basis_content = f.read()

for out_file in glob.glob("*.inp"):
    with open(out_file, "a", encoding="utf-8") as f:
        f.write("\n")          # separador seguro
        f.write(basis_content)

print("Basis inserted!!")
