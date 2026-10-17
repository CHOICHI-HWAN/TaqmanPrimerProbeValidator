import copy
import json
from pathlib import Path

import pytest

from virus_profile import ProfileError, chrom_name, load_profile, segment_keys, validate

ROOT = Path(__file__).resolve().parent.parent


def test_shipped_lassa_profile_is_valid_and_matches_historic_values():
    p = load_profile(ROOT / "profiles" / "lassa.json")
    assert segment_keys(p) == ["s", "l"]
    assert (p["segments"]["s"]["min_len"], p["segments"]["l"]["min_len"]) == (3000, 6000)
    assert (p["segments"]["s"]["reference_accession"], p["segments"]["l"]["reference_accession"]) == ("AY628203", "AY628202")
    assert p["population"]["pooled_label"] == "ALL_LASV" and chrom_name(p, "s") == "LASV_S"
    assert p["diversity"] == {"window_size": 1, "min_coverage": 0.95, "min_identity": 0.65, "min_lineage_n": 10}


def _good():
    return load_profile(ROOT / "profiles" / "lassa.json")


@pytest.mark.parametrize("mutate,msg", [
    (lambda p: p.pop("name"), "'name' is required"),
    (lambda p: p["source"].update(type="ftp"), "source.type"),
    (lambda p: p["source"].update(base_url="data/x"), "base_url"),
    (lambda p: p["segments"].clear(), "at least one segment"),
    (lambda p: p["segments"]["s"].update(min_len=0), "min_len"),
    (lambda p: p["segments"]["s"].pop("reference_accession"), "reference_accession"),
    (lambda p: p["segments"].update({"S": p["segments"]["s"]}), "lowercase"),
    (lambda p: p["metadata"].pop("country"), "metadata.country"),
    (lambda p: p["diversity"].update(min_identity=1.5), "min_coverage / min_identity"),
    (lambda p: p.update(profile_version=2), "profile_version"),
])
def test_invalid_profiles_are_rejected_with_a_clear_message(mutate, msg):
    p = copy.deepcopy(_good())
    mutate(p)
    with pytest.raises(ProfileError, match=msg):
        validate(p)


def test_missing_and_broken_files(tmp_path):
    with pytest.raises(ProfileError, match="not found"):
        load_profile(tmp_path / "nope.json")
    (tmp_path / "bad.json").write_text("{not json")
    with pytest.raises(ProfileError, match="valid JSON"):
        load_profile(tmp_path / "bad.json")


def test_diversity_block_is_optional(tmp_path):
    p = copy.deepcopy(_good()); p.pop("diversity")
    (tmp_path / "p.json").write_text(json.dumps(p))
    assert load_profile(tmp_path / "p.json")["diversity"]["min_lineage_n"] == 10


def test_profile_is_found_from_config_yaml(tmp_path):
    from virus_profile import find_profile_path, load_for_root
    (tmp_path / "config").mkdir()
    (tmp_path / "profiles").mkdir()
    (tmp_path / "profiles" / "x.json").write_text((ROOT / "profiles" / "lassa.json").read_text())
    (tmp_path / "config" / "config.yaml").write_text('mfeprimer_cpu: 4\nprofile: "profiles/x.json"   # which virus\n')
    assert find_profile_path(tmp_path) == tmp_path / "profiles" / "x.json"
    assert load_for_root(tmp_path)["id"] == "lassa"
    (tmp_path / "config" / "config.yaml").write_text("mfeprimer_cpu: 4\n")          # no key -> historic default
    assert find_profile_path(tmp_path) == tmp_path / "profiles" / "lassa.json"


def test_nipah_profile_flat_source_and_clades_from_auspice():
    from virus_profile import source_urls
    p = load_profile(ROOT / "profiles" / "nipah.json")
    assert segment_keys(p) == ["genome"]
    assert source_urls(p, "genome") == ("https://data.nextstrain.org/files/workflows/nipah/sequences.fasta.zst",
                                         "https://data.nextstrain.org/files/workflows/nipah/metadata.tsv.zst")
    assert p["lineage"]["mode"] == "auspice_clades" and p["lineage"]["attr"] == "clade_membership" and p["lineage"]["label"] == "Clade"
    assert chrom_name(p, "genome") == "NIV_GENOME" and p["qc"]["max_n_fraction"] == 0.01
    # the LASV profile keeps the segment sub-folder layout and the metadata-column lineage
    q = load_profile(ROOT / "profiles" / "lassa.json")
    assert source_urls(q, "s")[0].endswith("/lassa/s/sequences.fasta.zst") and q["lineage"]["mode"] == "metadata_column"


def test_auspice_mode_needs_a_path():
    p = copy.deepcopy(_good())
    p["lineage"] = {"mode": "auspice_clades"}
    with pytest.raises(ProfileError, match="lineage.auspice"):
        validate(p)
    p["lineage"] = {"mode": "nope"}
    with pytest.raises(ProfileError, match="lineage.mode"):
        validate(p)


def test_ui_wording_block_has_defaults_and_shipped_values():
    p = copy.deepcopy(_good())
    p.pop("ui", None)
    ui = validate(p)["ui"]
    assert ui["name_ko"] == p["short_name"] and ui["lineage_ko"] == "리니지" and ui["web_url"] == p["source"]["base_url"]
    n = load_profile(ROOT / "profiles" / "nipah.json")["ui"]
    assert n["lineage_ko"] == "클레이드" and n["web_url"].startswith("https://nextstrain.org/") and not n["orientation_note"]
    assert load_profile(ROOT / "profiles" / "lassa.json")["ui"]["lineage_ko"] == "리니지"
