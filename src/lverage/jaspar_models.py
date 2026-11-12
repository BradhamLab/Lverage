"""Cis-BP similarity-regression support adapted for Lverage.

Adapted from ``__init__.py`` in the JASPAR profile inference tool at
https://github.com/amathelier/JASPAR-inference-tool, source revision
ac4024148f1f140c7370d7b13a70556b4d13f7a3 (2026.1.0).

Copyright (c) 2019 Oriol Fornes and Xi Zhang. The original MIT terms are
preserved in ``jaspar_inference_license.txt``. Upstream attribution also names
Anthony Mathelier and the JASPAR Consortium.
"""

import json
from pathlib import Path

import numpy as np


PFAM_TO_CISBP = {
    "AP2": "AP2",
    "KilA-N": "APSES",
    "ARID": "ARID/BRIGHT",
    "AT_hook": "AT hook",
    "HLH": "bHLH",
    "bZIP_1": "bZIP",
    "zf-C2H2": "C2H2 ZF",
    "DM": "DM",
    "zf-Dof": "Dof",
    "E2F_TDP": "E2F",
    "Ets": "Ets",
    "Forkhead": "Forkhead",
    "GATA": "GATA",
    "GCM": "GCM",
    "Homeobox": "Homeodomain",
    "Pou": "Homeodomain,POU",
    "HSF_DNA-bind": "HSF",
    "MADF_DNA_bdg": "MADF",
    "Myb_DNA-binding": "Myb/SANT",
    "NAM": "NAC/NAM",
    "zf-C4": "Nuclear receptor",
    "HTH_psq": "Pipsqueak",
    "RFX_DNA_binding": "RFX",
    "SAND": "SAND",
    "SBP": "SBP",
    "HMG_box": "Sox",
    "T-box": "T-box",
    "TCP": "TCP",
    "TCR": "TCR/CxC",
    "WRKY": "WRKY",
    "Zn_clus": "Zinc cluster",
    None: "NO_THRESHOLD",
}
CISBP_TO_PFAM = {value: key for key, value in PFAM_TO_CISBP.items()}


def read_sr_model(path):
    """Load a Cis-BP SR model with the same array/NaN conversions upstream uses."""
    with Path(path).open(encoding="utf-8") as handle:
        model = json.load(handle)
    if "SR.Weights" in model:
        model["SR.Weights"] = np.asarray(model["SR.Weights"])
        model["SR.FeatureScales.mean"] = np.asarray(model["SR.FeatureScales.mean"])
        standard_deviation = np.asarray(model["SR.FeatureScales.sd"])
        standard_deviation[standard_deviation == 0] = np.nan
        model["SR.FeatureScales.sd"] = standard_deviation
    if np.isnan(model["Threshold.Dis"]):
        model["Threshold.Dis"] = None
    return model


def load_cisbp_models(data_dir):
    """Load family-keyed Cis-BP models exactly as the upstream loader does."""
    models = {}
    model_dir = Path(data_dir) / "cisbp"
    for filename in model_dir.iterdir():
        if filename.suffix != ".json":
            continue
        model = read_sr_model(filename)
        models.setdefault(CISBP_TO_PFAM[model["Family_Name"]], model)
    return models


def score_alignment_result(result, model, apply_identical_rule=True):
    """Apply the upstream HSim identity/SR thresholds without approximation."""
    if apply_identical_rule and result["PctID_L"] == 1:
        return result["PctID_L"], "HSim"

    if model["Model.Class"] == "SequenceIdentity":
        score = result[model["Model.Name"]]
        if score >= model["Threshold.HSim"]:
            return score, "HSim"
    else:
        score = result[model["Baseline"]["Name"]]
        if score >= model["Baseline"]["Threshold.HSim"]:
            return score, "HSim"

        feature_name = "ByPos." + model["SR.Features"].replace("_", ".")
        values = np.asarray(result[feature_name])
        if len(values) == len(model["SR.FeatureScales.mean"]):
            scaled = (values - model["SR.FeatureScales.mean"]) / model["SR.FeatureScales.sd"]
            scaled[np.isnan(scaled)] = 0
            score = model["SR.Intercept"] + np.dot(model["SR.Weights"], scaled)
            if model["SR.LogisticTransform"]:
                score = 1 / (1 + np.exp(-score))
            if score >= model["Threshold.HSim"]:
                return score, "HSim"

    return np.nan, np.nan
