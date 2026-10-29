#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import glob


def main():
    patterns = [
        "*.inp",
        "*in_calcgrid.xyz",
        "*in_symgrid.xyz",
        "*_oriented.xyz",
        "grid_info.txt",
        "*_vprofiles.xyz",
        "*_hprofiles.xyz",
        "*_prediXY*.xyz",
        "*_prediXY_path.txt",
        "*_rings.png",
        "*.csv",
        "*.png",
        "int_out.txt",
        "*.out",
        "*.cube",
        "*.vtk",
    ]

    files_to_delete = []
    for pattern in patterns:
        files_to_delete.extend(glob.glob(pattern))

    files_to_delete = sorted(set(files_to_delete))

    if not files_to_delete:
        print("No files found to delete.")
        return

    for f in files_to_delete:
        try:
            os.remove(f)
            print(f"Deleted: {f}")
        except Exception as e:
            print(f"Could not delete {f}: {e}")


if __name__ == "__main__":
    main()