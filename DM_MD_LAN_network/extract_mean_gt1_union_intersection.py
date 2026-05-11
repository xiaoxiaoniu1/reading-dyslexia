#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path
from typing import Set


BASE_DIR = Path("/data/home/tqi/data1/share/after_freesurfer")
INPUT_DIR = BASE_DIR / "FILE" / "test_mean_1.5" / "DM_MD_LAN_network"
OUTPUT_DIR = INPUT_DIR

INPUT_FILES = [
    INPUT_DIR / "reading_association-test_z_FDR_0.01_DK318.csv",
    INPUT_DIR / "v4-topics-400_355_reading_readers_dyslexia_association-test_z_FDR_0.01_DK318.csv",
    INPUT_DIR / "v5-topics-200_150_reading_phonological_readers_association-test_z_FDR_0.01_DK318.csv",
]

INTERSECTION_OUTPUT_TEMPLATE = "mean_gt{threshold_label}_intersection_reading_v4_v5.csv"
UNION_OUTPUT_TEMPLATE = "mean_gt{threshold_label}_union_reading_v4_v5.csv"
OUTPUT_HEADER = "barin_region"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract DK318 brain regions with mean above a threshold and compute union/intersection across three networks."
    )
    parser.add_argument(
        "--threshold",
        type=float,
        default=1.0,
        help="Keep regions with mean strictly greater than this value. Default: 1.0",
    )
    return parser.parse_args()


def format_threshold_label(threshold: float) -> str:
    return str(threshold).replace(".", "p")


def load_regions_with_mean_gt_threshold(csv_path: Path, threshold: float) -> Set[str]:
    regions = set()

    with csv_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        required_columns = {"hemi", "struct_name", "mean"}
        missing_columns = required_columns - set(reader.fieldnames or [])
        if missing_columns:
            raise ValueError(
                f"File {csv_path} is missing required columns: {sorted(missing_columns)}"
            )

        for row in reader:
            mean_value = float(row["mean"])
            if mean_value > threshold:
                region_name = f'{row["hemi"]}_{row["struct_name"]}'
                regions.add(region_name)

    return regions


def write_region_list(output_path: Path, regions: Set[str], header: str) -> None:
    sorted_regions = sorted(regions)
    with output_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow([header])
        for region in sorted_regions:
            writer.writerow([region])


def main() -> None:
    args = parse_args()
    threshold = args.threshold
    threshold_label = format_threshold_label(threshold)
    intersection_output = OUTPUT_DIR / INTERSECTION_OUTPUT_TEMPLATE.format(
        threshold_label=threshold_label
    )
    union_output = OUTPUT_DIR / UNION_OUTPUT_TEMPLATE.format(
        threshold_label=threshold_label
    )

    region_sets = [load_regions_with_mean_gt_threshold(path, threshold) for path in INPUT_FILES]

    intersection_regions = set.intersection(*region_sets)
    union_regions = set.union(*region_sets)

    write_region_list(intersection_output, intersection_regions, OUTPUT_HEADER)
    write_region_list(union_output, union_regions, OUTPUT_HEADER)

    print(f"Threshold: mean > {threshold}")
    for path, regions in zip(INPUT_FILES, region_sets):
        print(f"{path.name}: {len(regions)} regions")
    print(f"Intersection: {len(intersection_regions)} regions -> {intersection_output}")
    print(f"Union: {len(union_regions)} regions -> {union_output}")


if __name__ == "__main__":
    main()
