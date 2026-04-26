import os
import sys


def load_template(path: str) -> str:
    with open(path, "r", encoding="utf-8") as f:
        return f.read()


def main() -> None:
    template_path = "/data/home/tqi/data1/share/after_freesurfer/CODE/degree/degree_analysis_DK318_DGLM.py"
    script_txt = load_template(template_path)

    replacements = [
        ("DK-318", "DK-68"),
        ("DK318", "DK68"),
        ("degree_analysis_DK318_DGLM.py", "degree_analysis_DK68_DGLM.py"),
        ("MIND_DK318/DGLM", "MIND_DK68_DGLM"),
        ("MIND_DK318_DGLM", "MIND_DK68_DGLM"),
        ("DGLM_DK318_", "DGLM_DK68_"),
    ]
    for old, new in replacements:
        script_txt = script_txt.replace(old, new)

    script_txt = script_txt.replace(
        'LH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-318", "lh.DK318.annot")',
        'LH_ANNOT = None',
    )
    script_txt = script_txt.replace(
        'RH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-318", "rh.DK318.annot")',
        'RH_ANNOT = None',
    )
    script_txt = script_txt.replace(
        'LH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-68", "lh.DK68.annot")',
        'LH_ANNOT = None',
    )
    script_txt = script_txt.replace(
        'RH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK-68", "rh.DK68.annot")',
        'RH_ANNOT = None',
    )
    script_txt = script_txt.replace(
        'LH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK68", "lh.DK68.annot")',
        'LH_ANNOT = None',
    )
    script_txt = script_txt.replace(
        'RH_ANNOT = os.path.join(BASE_DIR, r"FILE/DK68", "rh.DK68.annot")',
        'RH_ANNOT = None',
    )
    script_txt = script_txt.replace(
        '    lh_labels, _, lh_names = read_annot(LH_ANNOT)',
        '    lh_labels, _, lh_names = read_annot(os.path.join(fs_dir, "label", "lh.aparc.annot"))',
    )
    script_txt = script_txt.replace(
        '    rh_labels, _, rh_names = read_annot(RH_ANNOT)',
        '    rh_labels, _, rh_names = read_annot(os.path.join(fs_dir, "label", "rh.aparc.annot"))',
    )

    dk68_globals = {
        "__name__": "__main__",
        "__file__": os.path.abspath(__file__),
        "sys": sys,
        "os": os,
    }
    exec(compile(script_txt, template_path, "exec"), dk68_globals)


if __name__ == "__main__":
    main()
