#!/bin/bash
#SBATCH --job-name=DK318_NEUROSYNTH
#SBATCH --output=/data/home/tqi/data1/share/after_freesurfer/CODE/DM_MD_LAN_network/logs/DK318_NEUROSYNTH_%j.out
#SBATCH --error=/data/home/tqi/data1/share/after_freesurfer/CODE/DM_MD_LAN_network/logs/DK318_NEUROSYNTH_%j.err
#SBATCH --partition=partition_1
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=72:00:00
#SBATCH --account=""

set -euo pipefail

echo "[INFO] Job script started"

fail() {
  echo "ERROR: $*" >&2
  exit 1
}

# FreeSurfer setup
export FREESURFER_HOME="${FREESURFER_HOME:-/data/software/freesurfer}"
export SUBJECTS_DIR="${SUBJECTS_DIR:-$FREESURFER_HOME/subjects}"
export FSFAST_HOME="${FSFAST_HOME:-$FREESURFER_HOME/fsfast}"
export FUNCTIONALS_DIR="${FUNCTIONALS_DIR:-$FREESURFER_HOME/sessions}"
export FSF_OUTPUT_FORMAT="${FSF_OUTPUT_FORMAT:-nii.gz}"
export MNI_DIR="${MNI_DIR:-$FREESURFER_HOME/mni}"
export MINC_BIN_DIR="${MINC_BIN_DIR:-$MNI_DIR/bin}"
export MINC_LIB_DIR="${MINC_LIB_DIR:-$MNI_DIR/lib}"
export MNI_DATAPATH="${MNI_DATAPATH:-$MNI_DIR/data}"
export LOCAL_DIR="${LOCAL_DIR:-$FREESURFER_HOME/local}"

[[ -d "$FREESURFER_HOME" ]] || fail "FREESURFER_HOME not found: $FREESURFER_HOME"
[[ -d "$FREESURFER_HOME/bin" ]] || fail "FreeSurfer bin dir not found: $FREESURFER_HOME/bin"

export PATH="$FREESURFER_HOME/bin:$FSFAST_HOME/bin:$MINC_BIN_DIR:${PATH:-/usr/bin:/bin}"

if [[ -d "$MINC_LIB_DIR" ]]; then
  export LD_LIBRARY_PATH="$MINC_LIB_DIR:${LD_LIBRARY_PATH:-}"
fi

if [[ -d "$FREESURFER_HOME/lib" ]]; then
  export LD_LIBRARY_PATH="$FREESURFER_HOME/lib:${LD_LIBRARY_PATH:-}"
fi

echo "[INFO] FreeSurfer environment initialized without SetUpFreeSurfer.sh"
echo "[INFO] Checking FreeSurfer commands in PATH"

# Input/output paths
INPUT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network"
OUTPUT_DIR="/data/home/tqi/data1/share/after_freesurfer/FILE/test_mean_1.5/DM_MD_LAN_network"
TRG_SUBJECT="fsaverage"
ANNOT_DIR="$SUBJECTS_DIR/$TRG_SUBJECT/label"
LH_ANNOT="$ANNOT_DIR/lh.DK318.annot"
RH_ANNOT="$ANNOT_DIR/rh.DK318.annot"

PROJECT_FRACTION="0.5"
INTERP="nearest"

MAPS=(
  "$INPUT_DIR/phonological_association-test_z_FDR_0.01.nii.gz"
  "$INPUT_DIR/reading_association-test_z_FDR_0.01.nii.gz"
  "$INPUT_DIR/sentence comprehension_association-test_z_FDR_0.01.nii.gz"
  "$INPUT_DIR/word form_association-test_z_FDR_0.01.nii.gz"
)

command -v mri_vol2surf >/dev/null 2>&1 || fail "mri_vol2surf not found"
command -v mri_segstats >/dev/null 2>&1 || fail "mri_segstats not found"
echo "[INFO] FreeSurfer commands available"
echo "[INFO] mri_vol2surf: $(command -v mri_vol2surf)"
echo "[INFO] mri_segstats: $(command -v mri_segstats)"

[[ -d "$INPUT_DIR" ]] || fail "Input dir not found: $INPUT_DIR"
[[ -d "$OUTPUT_DIR" ]] || fail "Output dir not found: $OUTPUT_DIR"
[[ -d "$ANNOT_DIR" ]] || fail "Annot dir not found: $ANNOT_DIR"
[[ -d "$SUBJECTS_DIR/$TRG_SUBJECT" ]] || fail "Target subject not found: $SUBJECTS_DIR/$TRG_SUBJECT"

[[ -f "$LH_ANNOT" ]] || fail "Left annot not found: $LH_ANNOT"
[[ -f "$RH_ANNOT" ]] || fail "Right annot not found: $RH_ANNOT"

TMP_DIR="$(mktemp -d)"
trap 'rm -rf "$TMP_DIR"' EXIT

for MAP in "${MAPS[@]}"; do
  [[ -f "$MAP" ]] || fail "Input map not found: $MAP"

  base="$(basename "$MAP")"
  stem="${base%.nii.gz}"
  csv_stem="${stem// /_}"
  csv_out="$OUTPUT_DIR/${csv_stem}_DK318.csv"

  echo "network,hemi,index,segid,nvertices,area_mm2,struct_name,mean,stddev,min,max,range" > "$csv_out"

  for hemi in lh rh; do
    if [[ "$hemi" == "lh" ]]; then
      annot="$LH_ANNOT"
    else
      annot="$RH_ANNOT"
    fi

    surf_out="$TMP_DIR/${csv_stem}.${hemi}.mgh"
    sum_out="$TMP_DIR/${csv_stem}.${hemi}.sum"

    echo "[INFO] Projecting: $base -> $hemi"
    mri_vol2surf \
      --mov "$MAP" \
      --mni152reg \
      --trgsubject "$TRG_SUBJECT" \
      --hemi "$hemi" \
      --projfrac "$PROJECT_FRACTION" \
      --interp "$INTERP" \
      --out_type mgh \
      --o "$surf_out"

    echo "[INFO] Parcel summary: $base -> $hemi"
    mri_segstats \
      --annot "$TRG_SUBJECT" "$hemi" "$annot" \
      --i "$surf_out" \
      --sum "$sum_out" \
      --excludeid 0 \
      --no-global-stats >/dev/null

    awk -v network="$stem" -v hemi="$hemi" '
      BEGIN { OFS = "," }
      /^#/ { next }
      NF < 10 { next }
      {
        index_col = $1
        segid = $2
        nvertices = $3
        area = $4
        mean = $(NF-4)
        stddev = $(NF-3)
        minv = $(NF-2)
        maxv = $(NF-1)
        rangev = $NF

        struct_name = ""
        for (i = 5; i <= NF-5; i++) {
          struct_name = struct_name (i == 5 ? "" : " ") $i
        }

        if (tolower(struct_name) == "unknown") {
          next
        }

        gsub(/"/, "\"\"", network)
        gsub(/"/, "\"\"", hemi)
        gsub(/"/, "\"\"", struct_name)

        print "\"" network "\"", "\"" hemi "\"", index_col, segid, nvertices, area, "\"" struct_name "\"", mean, stddev, minv, maxv, rangev
      }
    ' "$sum_out" >> "$csv_out"
  done

  echo "[DONE] $csv_out"
done

echo "All networks finished."
