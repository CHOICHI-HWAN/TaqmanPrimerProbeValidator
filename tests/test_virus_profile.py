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


def test_shipped_profile_pathoplexus_source_and_country_lineage():
    from virus_profile import source_urls
    p = load_profile(ROOT / "profiles" / "machv.json")
    assert segment_keys(p) == ['genome']
    assert source_urls(p, "genome") == ("https://api.ncbi.nlm.nih.gov/datasets/v2", "3052317")
    assert p["lineage"]["mode"] == "metadata_column" and p["segments"]["genome"]["lineage_col"] == "country_lineage"
    assert chrom_name(p, "genome") == "MACHV_GENOME" and p["qc"]["max_n_fraction"] == 0.01
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
    n = load_profile(ROOT / "profiles" / "machv.json")["ui"]
    assert n["lineage_ko"] == "발생 국가" and n["web_url"].startswith("https://") and not n["orientation_note"]
    assert load_profile(ROOT / "profiles" / "lassa.json")["ui"]["lineage_ko"] == "리니지"


def test_pathoplexus_source_needs_an_organism_and_yields_lapis_target():
    from virus_profile import source_urls
    p = copy.deepcopy(_good())
    p["source"] = {"type": "pathoplexus", "base_url": "https://lapis.pathoplexus.org"}
    with pytest.raises(ProfileError, match="source.organism"):
        validate(p)
    p["source"]["organism"] = "ebola-sudan"
    assert source_urls(validate(p), "genome" if "genome" in p["segments"] else next(iter(p["segments"]))) == ("https://lapis.pathoplexus.org", "ebola-sudan")


def test_ncbi_source_needs_a_taxon():
    from virus_profile import source_urls
    p = copy.deepcopy(_good())
    p["source"] = {"type": "ncbi", "base_url": "https://api.ncbi.nlm.nih.gov/datasets/v2"}
    with pytest.raises(ProfileError, match="source.taxon"):
        validate(p)
    p["source"]["taxon"] = 186539
    assert source_urls(validate(p), next(iter(p["segments"]))) == ("https://api.ncbi.nlm.nih.gov/datasets/v2", "186539")


def test_pathoplexus_segmented_organism_gets_the_segment_in_its_target():
    from virus_profile import source_urls
    p = copy.deepcopy(_good())
    p["source"] = {"type": "pathoplexus", "base_url": "https://lapis.pathoplexus.org", "organism": "cchf", "segment_subdirs": True}
    assert source_urls(validate(p), next(iter(p["segments"])))[1].startswith("cchf/")


def test_pathoplexus_lineage_field_is_part_of_the_target():
    from virus_profile import source_urls
    p = copy.deepcopy(_good())
    p["source"] = {"type": "pathoplexus", "base_url": "https://lapis.pathoplexus.org", "organism": "mpox", "lineage_field": "clade"}
    assert source_urls(validate(p), next(iter(p["segments"])))[1] == "mpox@clade"


def test_pathoplexus_sublineage_field_is_appended_to_the_target():
    from virus_profile import source_urls
    p = copy.deepcopy(_good())
    p["source"] = {"type": "pathoplexus", "base_url": "https://lapis.pathoplexus.org", "organism": "mpox", "lineage_field": "clade", "sublineage_field": "lineage"}
    assert source_urls(validate(p), next(iter(p["segments"])))[1] == "mpox@clade:lineage"


def test_ncbi_segmented_virus_gets_the_segment_in_its_target():
    from virus_profile import source_urls
    p = copy.deepcopy(_good())
    p["source"] = {"type": "ncbi", "base_url": "https://api.ncbi.nlm.nih.gov/datasets/v2", "taxon": 11588, "segment_subdirs": True}
    assert source_urls(validate(p), "s")[1] == "11588/S"


def test_ncbi_several_taxa_with_filters_build_one_spec_per_segment():
    from virus_profile import source_urls
    p = copy.deepcopy(_good())
    seg = next(iter(p["segments"]))
    p["segments"][seg].update(ncbi_segment="4", ncbi_subtype="H3")
    p["source"] = {"type": "ncbi", "base_url": "https://api.ncbi.nlm.nih.gov/datasets/v2", "taxon": [114727, 119210], "since": 2015, "host_filter": {"114727": 9606},
                   "group_by": "subtype", "segment_subdirs": True}
    assert source_urls(validate(p), seg)[1] == "114727,119210/4?since=2015&host=114727:9606&sub=H3&group=subtype"
