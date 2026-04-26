#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Wrapper for DK68 stratified DGLM degree brainmaps.
It reuses the DK318 plotting implementation and swaps DK68-specific
paths, file patterns, and annotation locations.
"""

from pathlib import Path

TEMPLATE = Path("/data/home/tqi/data1/share/after_freesurfer/CODE/degree/plot_significant_DK318_stratified_brainmaps.py")
script_txt = TEMPLATE.read_text()

replacements = [
    ("Plot DK318 stratified DGLM degree brainmaps.", "Plot DK68 stratified DGLM degree brainmaps."),
    ("Plot stratified DK318 DGLM degree brainmaps.", "Plot stratified DK68 DGLM degree brainmaps."),
    ("FILE/DK-318/lh.DK318.annot", "FILE/fsaverage/label/lh.aparc.annot"),
    ("FILE/DK-318/rh.DK318.annot", "FILE/fsaverage/label/rh.aparc.annot"),
    ("MIND_DK318_DGLM_stratified", "MIND_DK68_DGLM_stratified"),
    ("DGLM_DK318_", "DGLM_DK68_"),
    ("DK318", "DK68"),
]

for old, new in replacements:
    script_txt = script_txt.replace(old, new)

namespace = {"__name__": "__main__", "__file__": str(TEMPLATE)}
exec(compile(script_txt, str(TEMPLATE), "exec"), namespace)
