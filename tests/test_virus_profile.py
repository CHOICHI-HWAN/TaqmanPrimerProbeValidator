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
