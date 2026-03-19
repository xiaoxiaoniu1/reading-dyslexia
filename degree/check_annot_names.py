import os
import pandas as pd
import numpy as np
from nibabel.freesurfer.io import read_annot

LH_ANNOT = r"e:\阅读障碍\DK-318\lh.DK318.annot"
RH_ANNOT = r"e:\阅读障碍\DK-318\rh.DK318.annot"
ROI_FILE = r"e:\阅读障碍\DK-318\DK318_roi_names.csv"

def check():
    # Load ROI names
    roi_df = pd.read_csv(ROI_FILE)
    roi_names = roi_df["region"].tolist()
    print(f"Total ROIs in CSV: {len(roi_names)}")
    print(f"First 5 ROIs: {roi_names[:5]}")

    # Load Annotations
    lh_labels, lh_ctab, lh_names = read_annot(LH_ANNOT)
    rh_labels, rh_ctab, rh_names = read_annot(RH_ANNOT)
    
    # Decode bytes to strings
    lh_names_str = [n.decode('utf-8') for n in lh_names]
    rh_names_str = [n.decode('utf-8') for n in rh_names]
    
    print(f"\nLH Annot Names (first 5): {lh_names_str[:5]}")
    print(f"RH Annot Names (first 5): {rh_names_str[:5]}")
    
    # Check matching
    # Construct expected names from annot (usually we prepend lh_ or rh_)
    
    # Note: 'unknown' or 'corpuscallosum' might be in annot but not in ROI list
    
    print(f"\nLH Annot Names count: {len(lh_names_str)}")
    print(f"RH Annot Names count: {len(rh_names_str)}")

if __name__ == "__main__":
    check()