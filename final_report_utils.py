from collections import defaultdict
import ast
from unittest import result
from typing import Any, Callable, Iterable, Optional
from structures import (
    geo_syns,
    taxa_syns,
    muts_interest,
    seg_lens,
    muts_loci_meaning,
    iav_segments,
    segment_to_number_iav,
)
import pandas as pd
import os
import csv
import zipfile
from metadata_tree_utils import load_reference_dicts
from flu_utils import parse_clstr
from html import escape
import re


def _norm_key(value: Any) -> Optional[str]:
    """
    Normalize a metadata/reference key for case-insensitive lookup.

    Uses lower() consistently. This function is used both for raw metadata
    values and for canonical reference dictionary keys.

    Returns
    -------
    str or None
        Lowercase stripped string, or None for missing/empty values.
    """
    if value is None:
        return None

    try:
        if pd.isna(value):
            return None
    except Exception:
        pass

    text = str(value).strip()

    if not text:
        return None

    return text.lower()


def _build_lower_lookup(values: Iterable[str]) -> dict[str, str]:
    """
    Build a lowercase lookup index from canonical reference names.

    Example
    -------
    {"homo sapiens": "Homo sapiens"}

    If duplicate lowercase keys occur, the first one is kept.
    """
    lookup = {}

    for value in values:
        key = _norm_key(value)

        if key is None:
            continue

        if key not in lookup:
            lookup[key] = str(value).strip()

    return lookup


def _resolve_with_synonyms(
    value: Any,
    synonyms: Optional[dict[str, Optional[str]]],
) -> Optional[str]:
    """
    Resolve a raw value through a synonym dictionary using lowercase matching.

    Both the raw input and synonym keys are compared as lowercase strings.

    Behaviour
    ---------
    - None / NaN / empty -> None
    - if lower(value) matches lower(synonym_key), return synonym value
    - synonym value may be None; downstream code must treat None as Unknown
    - if no synonym matches, return stripped original value

    Notes
    -----
    This function does not require synonym dictionaries to already have
    lowercase keys. It normalizes both sides.
    """
    value_key = _norm_key(value)

    if value_key is None:
        return None

    raw_value = str(value).strip()

    if not synonyms:
        return raw_value

    for syn_key, canonical_value in synonyms.items():
        syn_key_norm = _norm_key(syn_key)

        if syn_key_norm is None:
            continue

        if value_key == syn_key_norm:
            if canonical_value is None:
                return None

            canonical_text = str(canonical_value).strip()

            if not canonical_text:
                return None

            return canonical_text

    return raw_value


def _resolve_against_reference(
    value: Any,
    reference_keys: Iterable[str],
    synonyms: Optional[dict[str, Optional[str]]] = None,
) -> Optional[str]:
    """
    Resolve a value through synonyms and then against canonical reference keys.

    Resolution order
    ----------------
    1. raw value -> synonyms
    2. None from synonyms -> None
    3. lowercase comparison against canonical reference keys
    4. if found, return canonical reference spelling
    5. if not found, return resolved value as-is

    This keeps the canonical spelling from geo_dict/taxa_dict whenever possible.
    """
    resolved = _resolve_with_synonyms(value, synonyms)

    if resolved is None:
        return None

    reference_lookup = _build_lower_lookup(reference_keys)
    resolved_key = _norm_key(resolved)

    if resolved_key is None:
        return None

    if resolved_key in reference_lookup:
        return reference_lookup[resolved_key]

    return str(resolved).strip()


def _unknown_if_none(value: Optional[str]) -> str:
    """
    Convert None/empty values to the report-level Unknown label.
    """
    if value is None:
        return "Unknown"

    text = str(value).strip()

    if not text:
        return "Unknown"

    return text


def _normalize_collection_date(value: Any) -> str:
    """Return a stable display string for collection dates."""
    try:
        if pd.isna(value):
            return "Unknown"
    except Exception:
        pass

    text = str(value).strip()
    if not text or text.lower() == 'nan':
        return "Unknown"

    return text

def _empty_card(message: str) -> str:
    """Return a consistent empty-state card for report sections."""
    return f"""
    <div class="card">
      <p class="muted">{escape(str(message))}</p>
    </div>
    """


def _safe_read_table(path: str, required_cols: Optional[Iterable[str]] = None) -> tuple[Optional[pd.DataFrame], Optional[str]]:
    """Read a TSV defensively and return a dataframe or a user-facing error message."""
    if not path or not os.path.exists(path):
        return None, f"Required report input is missing: {path}"

    try:
        df = pd.read_table(path, index_col=False)
    except Exception as exc:
        return None, f"Could not read report input {path}: {exc}"

    if df is None or df.empty:
        return None, f"Report input is empty: {path}"

    if required_cols:
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            return None, f"Report input is missing required columns ({', '.join(missing)}): {path}"

    return df, None


def _safe_read_metadata_csv(path: str, required_cols: Optional[Iterable[str]] = None) -> tuple[Optional[pd.DataFrame], Optional[str]]:
    """Read the semicolon-delimited metadata CSV defensively."""
    if not path or not os.path.exists(path):
        return None, f"Required metadata input is missing: {path}"

    try:
        df = pd.read_csv(path, sep=';', index_col=False, low_memory=False)
    except Exception as exc:
        return None, f"Could not read metadata input {path}: {exc}"

    if df is None or df.empty:
        return None, f"Metadata input is empty: {path}"

    df.columns = [str(col).strip().upper() for col in df.columns]

    if required_cols:
        missing = [col for col in required_cols if col not in df.columns]
        if missing:
            return None, f"Metadata input is missing required columns ({', '.join(missing)}): {path}"

    return df, None


def _coerce_metadata_counter(cell: Any) -> dict[str, float]:
    """Parse a metadata distribution cell into a numeric dictionary."""
    if isinstance(cell, dict):
        source = dict(cell)
    elif pd.isna(cell):
        return {}
    else:
        text = str(cell).strip().replace('%', '')
        if not text:
            return {}
        try:
            source = ast.literal_eval(text)
        except Exception:
            return {}

    if not isinstance(source, dict):
        return {}

    cleaned = {}
    for key, value in source.items():
        try:
            cleaned[str(key).strip()] = float(value)
        except Exception:
            continue

    return cleaned


def _normalize_distribution_for_export(counter: dict[str, float], ndigits: int = 2) -> dict[str, float]:
    """Round a distribution deterministically while forcing the total to sum to 1.0."""
    if not counter:
        return {}

    cleaned = {}
    total = 0.0
    for key, value in counter.items():
        key_text = str(key).strip()
        try:
            value_num = float(value)
        except Exception:
            continue
        if not key_text or value_num <= 0:
            continue
        cleaned[key_text] = value_num
        total += value_num

    if not cleaned or total <= 0:
        return {}

    scale = 10 ** int(ndigits)
    ordered = sorted(cleaned.items(), key=lambda item: (-float(item[1]), str(item[0])))
    normalized = [(key, value / total) for key, value in ordered]
    base_units = []
    used_units = 0

    for key, value in normalized:
        exact_units = value * scale
        floor_units = int(exact_units)
        remainder = exact_units - floor_units
        base_units.append([key, floor_units, remainder])
        used_units += floor_units

    remaining_units = scale - used_units
    if remaining_units > 0:
        for key, floor_units, remainder in sorted(base_units, key=lambda item: (-item[2], str(item[0]))):
            if remaining_units <= 0:
                break
            for entry in base_units:
                if entry[0] == key:
                    entry[1] += 1
                    remaining_units -= 1
                    break

    rounded = {}
    for key, units, _ in base_units:
        rounded[key] = round(units / scale, ndigits)

    return rounded


def _serialize_distribution(counter: dict[str, float], ndigits: int = 2) -> str:
    """Serialize a distribution as a stable dict-style string."""
    normalized = _normalize_distribution_for_export(counter, ndigits=ndigits)
    if not normalized:
        return 'Unknown'

    parts = []
    for key, value in sorted(normalized.items(), key=lambda item: (-float(item[1]), str(item[0]))):
        parts.append(f"{repr(str(key))}: {value:.{ndigits}f}")
    return '{' + ', '.join(parts) + '}'


def _static_rollup_distribution(cell: Any, mapper: Callable[[str], str]) -> str:
    """Convert a metadata distribution into a serialized rolled-up distribution."""
    counter = _coerce_metadata_counter(cell)
    if not counter:
        return 'Unknown'

    rolled = _rollup_counts(counter, mapper, drop_unknown=True)
    if not rolled:
        return 'Unknown'

    return _serialize_distribution(rolled, ndigits=2)


def _collect_report_clusters(id_report_fps: Iterable[str]) -> list[str]:
    """Return the unique real cluster labels present in one or more ID reports."""
    seen = set()
    ordered = []

    for id_report_fp in id_report_fps:
        df, error = _safe_read_table(id_report_fp, required_cols=['CLUSTER'])
        if error or df is None:
            continue

        for value in df['CLUSTER'].tolist():
            cluster = str(value).strip()
            if not cluster.startswith('>Cluster '):
                continue
            if cluster in seen:
                continue
            seen.add(cluster)
            ordered.append(cluster)

    return ordered


def export_id_report_rollup(
    id_report_fp: str,
    output_fp: str,
    geo_dict: dict[str, tuple[str, str]],
    taxa_dict: dict[str, list[str]],
) -> None:
    """Write a per-sample ID report TSV with static taxonomy/geography roll-ups."""
    df, error = _safe_read_table(
        id_report_fp,
        required_cols=['SAMPLE_NAME', 'HOST', 'COUNTRY', 'ASSIGNED_BY'],
    )
    if error or df is None:
        raise ValueError(error or f'Could not load ID report: {id_report_fp}')

    for col in ['HOST', 'COUNTRY', 'ASSIGNED_BY']:
        if col not in df.columns:
            df[col] = 'Unknown'

    host_family = []
    host_class = []
    country_region = []
    country_continent = []

    for row in df.itertuples(index=False):
        assigned_by = str(getattr(row, 'ASSIGNED_BY', '')).strip()
        host_value = getattr(row, 'HOST', 'Unknown')
        country_value = getattr(row, 'COUNTRY', 'Unknown')

        if assigned_by == 'L-BLAST':
            host_family.append(rollup_taxa(host_value, taxa_dict, level='Family'))
            host_class.append(rollup_taxa(host_value, taxa_dict, level='Class'))
            country_region.append(rollup_geo(country_value, geo_dict, mode='region'))
            country_continent.append(rollup_geo(country_value, geo_dict, mode='continent'))
        else:
            host_family.append(
                _static_rollup_distribution(host_value, lambda term: rollup_taxa(term, taxa_dict, level='Family'))
            )
            host_class.append(
                _static_rollup_distribution(host_value, lambda term: rollup_taxa(term, taxa_dict, level='Class'))
            )
            country_region.append(
                _static_rollup_distribution(country_value, lambda loc: rollup_geo(loc, geo_dict, mode='region'))
            )
            country_continent.append(
                _static_rollup_distribution(country_value, lambda loc: rollup_geo(loc, geo_dict, mode='continent'))
            )

    df['HOST_FAMILY'] = host_family
    df['HOST_CLASS'] = host_class
    df['COUNTRY_REGION'] = country_region
    df['COUNTRY_CONTINENT'] = country_continent
    df.to_csv(output_fp, sep='\t', index=False)


def export_cluster_composition(
    id_report_fps: Iterable[str],
    metadata_csv_fp: str,
    cluster_clstr_fp: str,
    output_fp: str,
) -> None:
    """Write a cluster composition TSV for the clusters used in an analysis."""
    clusters = _collect_report_clusters(id_report_fps)

    metadata_df, error = _safe_read_metadata_csv(
        metadata_csv_fp,
        required_cols=['ACCESSION', 'SEGMENT', 'GENOTYPE', 'HOST', 'COUNTRY'],
    )
    if error or metadata_df is None:
        raise ValueError(error or f'Could not load metadata CSV: {metadata_csv_fp}')

    if not os.path.exists(cluster_clstr_fp):
        raise ValueError(f'Required cluster membership input is missing: {cluster_clstr_fp}')

    collection_date_col = 'COLLECTION_DATE' if 'COLLECTION_DATE' in metadata_df.columns else None
    metadata_df['ACCESSION'] = metadata_df['ACCESSION'].astype(str).str.strip()
    metadata_df = metadata_df.dropna(subset=['ACCESSION'])
    metadata_df = metadata_df[metadata_df['ACCESSION'] != '']
    metadata_df = metadata_df.drop_duplicates(subset=['ACCESSION'], keep='first')
    metadata_df = metadata_df.set_index('ACCESSION', drop=False)

    cluster_members = parse_clstr(cluster_clstr_fp, access_only=True)
    rows = []
    seen = set()

    for cluster in clusters:
        members = cluster_members.get(cluster, [])
        for accession in members:
            accession_text = str(accession).strip()
            key = (cluster, accession_text)
            if not accession_text or key in seen:
                continue
            seen.add(key)

            if accession_text not in metadata_df.index:
                continue

            meta_row = metadata_df.loc[accession_text]
            rows.append({
                'CLUSTER': cluster,
                'ACCESSION': accession_text,
                'SEGMENT': str(meta_row.get('SEGMENT', '')).strip(),
                'GENOTYPE': str(meta_row.get('GENOTYPE', '')).strip(),
                'HOST': str(meta_row.get('HOST', '')).strip(),
                'COUNTRY': str(meta_row.get('COUNTRY', '')).strip(),
                'COLLECTION_DATE': _normalize_collection_date(meta_row.get(collection_date_col, 'Unknown')) if collection_date_col else 'Unknown',
            })

    out_df = pd.DataFrame(
        rows,
        columns=['CLUSTER', 'ACCESSION', 'SEGMENT', 'GENOTYPE', 'HOST', 'COUNTRY', 'COLLECTION_DATE'],
    )
    out_df.to_csv(output_fp, sep='\t', index=False)

def rollup_geo(
    location: Any,
    geo_dict: dict[str, tuple[str, str]],
    mode: str = 'region'
) -> str:
    """
    Roll up a country name to a broader geographic label.

    Applies geo_syns before consulting the canonical geo_dict.

    Unknown/missing values mapped to None by geo_syns are returned as
    "Unknown".
    """
    mode = str(mode).strip().lower()

    canonical_location = _resolve_against_reference(
        value=location,
        reference_keys=geo_dict.keys(),
        synonyms=geo_syns,
    )

    if canonical_location is None:
        return "Unknown"

    if mode == "country":
        if canonical_location in geo_dict:
            return canonical_location
        return "Unknown"

    if canonical_location in geo_dict:
        if mode == "region":
            return geo_dict[canonical_location][0]
        if mode == "continent":
            return geo_dict[canonical_location][1]

    return "Unknown"

def search_tax_level(term: Any, taxa_dict: dict[str, list[str]]) -> str | tuple[str, str]:
    """
    Find the taxonomy level represented by a host term.

    Applies taxa_syns before consulting the canonical taxa_dict.

    Returns
    -------
    tuple[str, str] or str
        (level_name, canonical_key) when found, otherwise "Unknown".
    """
    canonical_term = _resolve_against_reference(
        value=term,
        reference_keys=taxa_dict.keys(),
        synonyms=taxa_syns,
    )

    if canonical_term is None:
        return "Unknown"

    canonical_term = str(canonical_term).strip()

    if not canonical_term:
        return "Unknown"

    if "sp." in canonical_term or "spp." in canonical_term:
        canonical_term = (
            canonical_term
            .replace("sp.", "")
            .replace("spp.", "")
            .strip()
        )

    if not canonical_term:
        return "Unknown"

    tax_level = {
        0: 'Species',
        1: 'Genus',
        2: 'Subfamily',
        3: 'Family',
        4: 'Division',
        5: 'Order',
        6: 'Class'
    }

    canonical_key_lookup = _build_lower_lookup(taxa_dict.keys())
    canonical_term_key = _norm_key(canonical_term)

    if canonical_term_key is None:
        return "Unknown"

    # Fast path: term matches a canonical taxa_dict key, ignoring case.
    if canonical_term_key in canonical_key_lookup:
        canonical_key = canonical_key_lookup[canonical_term_key]
        lineage = taxa_dict.get(canonical_key, [])

        for level, value in enumerate(lineage):
            if value == "X":
                continue

            value_key = _norm_key(value)

            if value_key == canonical_term_key:
                return tax_level[level], canonical_key

    # Fallback: term may match a lineage value rather than a taxa_dict key.
    for key, lineage in taxa_dict.items():
        for level, value in enumerate(lineage):
            if value == "X":
                continue

            if _norm_key(value) == canonical_term_key:
                return tax_level[level], key

    return "Unknown"

def rollup_taxa(term: Any, taxa_dict: dict[str, list[str]], level: str = 'Genus') -> str:
    """
    Roll up a host term to a requested taxonomy level.

    Applies taxa_syns before consulting the canonical taxa_dict.

    Unknown/missing values mapped to None by taxa_syns are returned as
    "Unknown".
    """
    level_tax = {
        'SPECIES': 0,
        'GENUS': 1,
        'SUBFAMILY': 2,
        'FAMILY': 3,
        'DIVISION': 4,
        'ORDER': 5,
        'CLASS': 6
    }

    result = search_tax_level(term, taxa_dict)

    if result == "Unknown":
        return "Unknown"

    t_lev, key = result

    if key not in taxa_dict:
        return "Unknown"

    target_level = str(level).upper()
    current_level = str(t_lev).upper()

    if target_level not in level_tax or current_level not in level_tax:
        return "Unknown"

    target_idx = level_tax[target_level]
    current_idx = level_tax[current_level]

    if current_idx <= target_idx:
        value = taxa_dict[key][target_idx]
    else:
        value = taxa_dict[key][current_idx]

    if value == "X" or not str(value).strip():
        return "Unknown"

    return str(value).strip()
    
def _rollup_counts(counter: Optional[dict[str, int | float]], mapper: Callable[[str], str], drop_unknown: bool = False) -> dict[str, int | float]:
    """
    Generic aggregator: maps each key via `mapper(key) -> new_key` and sums values.
    If `drop_unknown` is True, entries that map to 'Unknown' are discarded.
    """
    if counter is None:
        return {}

    if not isinstance(counter, dict):
        return {}

    out = defaultdict(float)

    for k, v in counter.items():
        try:
            new_k = mapper(k)
        except Exception:
            new_k = 'Unknown'

        if new_k == 'Unknown' and drop_unknown:
            continue

        try:
            out[new_k] += float(v)
        except Exception:
            continue

    if all(isinstance(v, int) for v in counter.values()):
        return {k: int(v) for k, v in out.items()}

    return dict(out)

def rollup_counts_geo(counter: dict[str, int | float], geo_dict: dict[str, tuple[str, str]], mode: str = 'region', drop_unknown: bool = False) -> dict[str, int | float]:
    """
    Aggregate geographic counts at a broader location level.

    Parameters
    ----------
    counter : dict
        Mapping of location labels to counts or proportions.
    geo_dict : dict
        Geographic lookup table.
    mode : {"country", "region", "continent"}, default "region"
        Geographic level to aggregate to.
    drop_unknown : bool, default False
        Whether entries that map to ``"Unknown"`` should be discarded.

    Returns
    -------
    dict
        Aggregated counts or proportions keyed by the selected geographic
        level.
    """
    return _rollup_counts(
        counter,
        lambda loc: rollup_geo(loc, geo_dict, mode=mode),
        drop_unknown=drop_unknown
    )
def rollup_counts_taxa(counter: dict[str, int | float], taxa_dict: dict[str, list[str]], level: str = 'Genus', drop_unknown: bool = False) -> dict[str, int | float]:
    """
    Aggregate host counts at a broader taxonomy level.

    Parameters
    ----------
    counter : dict
        Mapping of host labels to counts or proportions.
    taxa_dict : dict
        Taxonomy lookup table.
    level : str, default "Genus"
        Taxonomy level to aggregate to.
    drop_unknown : bool, default False
        Whether entries that map to ``"Unknown"`` should be discarded.

    Returns
    -------
    dict
        Aggregated counts or proportions keyed by the selected taxonomy level.
    """
    return _rollup_counts(
        counter,
        lambda term: rollup_taxa(term, taxa_dict, level=level),
        drop_unknown=drop_unknown
    )
def normalize(counter: Optional[dict[str, int | float]]) -> dict[str, float]:
    """
    Normalize numeric values so they sum to 1.

    Parameters
    ----------
    counter : dict
        Mapping of labels to numeric values.

    Returns
    -------
    dict
        Normalized values. Invalid inputs return an empty dictionary.
    """
    if counter is None or not isinstance(counter, dict):
        return {}

    total = sum(counter.values())
    if total == 0:
        return {k: 0.0 for k in counter}

    return {k: v / total for k, v in counter.items()}

def clean_host(host_cell: Any, taxa_dict: dict[str, list[str]], level: str = 'Species', drop_unknown: bool = True) -> dict[str, float]:
    """
    Parse and normalize host distributions from a report cell.

    Parameters
    ----------
    host_cell : dict or str or Any
        Host distribution stored as a dictionary or a string representation of
        one.
    taxa_dict : dict
        Taxonomy lookup table.
    level : str, default "Species"
        Taxonomy level used for roll-up.
    drop_unknown : bool, default True
        Whether unresolved host labels should be removed.

    Returns
    -------
    dict
        Cleaned host proportions keyed by taxonomy label.
    """
    if isinstance(host_cell, dict):
        host_dict = dict(host_cell)
    elif pd.isna(host_cell):
        return {}
    else:
        host_cell = str(host_cell).strip().replace('%', '')
        if not host_cell:
            return {}
        try:
            host_dict = ast.literal_eval(host_cell)
        except Exception:
            return {}

    if not isinstance(host_dict, dict):
        return {}

    cleaned = {}
    for key, value in host_dict.items():
        try:
            cleaned[str(key).strip()] = float(value)
        except Exception:
            continue

    cleaned = rollup_counts_taxa(cleaned, taxa_dict, level=level, drop_unknown=drop_unknown)

    if round(sum(cleaned.values()), 3) > 1.1:
        for key in cleaned:
            cleaned[key] = round(cleaned[key] / 100, 3)
    else:
        for key in cleaned:
            cleaned[key] = round(cleaned[key], 3)

    return cleaned

def clean_country(country_cell: Any, geo_dict: dict[str, tuple[str, str]], mode: str = 'country', drop_unknown: bool = True) -> dict[str, float]:
    """
    Parse and normalize geographic distributions from a report cell.

    Parameters
    ----------
    country_cell : dict or str or Any
        Geographic distribution stored as a dictionary or a string
        representation of one.
    geo_dict : dict
        Geographic lookup table.
    mode : {"country", "region", "continent"}, default "country"
        Geographic level to return.
    drop_unknown : bool, default True
        Whether unresolved locations should be removed during roll-up.

    Returns
    -------
    dict
        Cleaned geographic proportions keyed by the requested level.
    """
    if isinstance(country_cell, dict):
        country_dict = dict(country_cell)
    elif pd.isna(country_cell):
        return {}
    else:
        country_cell = str(country_cell).strip().replace('%', '')
        if not country_cell:
            return {}
        try:
            country_dict = ast.literal_eval(country_cell)
        except Exception:
            return {}

    if not isinstance(country_dict, dict):
        return {}

    cleaned = {}
    for key, value in country_dict.items():
        try:
            cleaned[str(key).strip()] = float(value)
        except Exception:
            continue

    cleaned = rollup_counts_geo(
            cleaned,
            geo_dict,
            mode=mode,
            drop_unknown=drop_unknown,
            )

    if sum(cleaned.values()) > 1.0:
        for key in cleaned:
            cleaned[key] = round(cleaned[key] / 100, 3)
    else:
        for key in cleaned:
            cleaned[key] = round(cleaned[key], 3)

    return cleaned

def clean_genotype(genotype_cell: Any) -> dict[str, float]:
    """
    Parse and normalize genotype distributions from a report cell.

    Parameters
    ----------
    genotype_cell : dict or str or Any
        Genotype distribution stored as a dictionary or a string
        representation of one.

    Returns
    -------
    dict
        Cleaned genotype proportions.
    """
    if isinstance(genotype_cell, dict):
        genotype_dict = dict(genotype_cell)
    elif pd.isna(genotype_cell):
        return {}
    else:
        genotype_cell = str(genotype_cell).strip().replace('%', '')
        if not genotype_cell:
            return {}
        try:
            genotype_dict = ast.literal_eval(genotype_cell)
        except Exception:
            return {}

    if not isinstance(genotype_dict, dict):
        return {}

    cleaned = {}
    for key, value in genotype_dict.items():
        try:
            cleaned[str(key).strip()] = float(value)
        except Exception:
            continue

    if sum(cleaned.values()) > 1.0:
        for key in cleaned:
            cleaned[key] = round(cleaned[key] / 100, 3)
    else:
        for key in cleaned:
            cleaned[key] = round(cleaned[key], 3)

    return cleaned

def convert_muts(flags: dict, muts_loci_meaning: dict) -> dict[str, str]:
    """
    Convert recorded mutation hits into report-ready labels and meanings.

    Parameters
    ----------
    flags : dict
        Pipeline state dictionary containing per-segment mutation hits.
    muts_loci_meaning : dict
        Mapping from mutation identifiers to display label and biological
        meaning.

    Returns
    -------
    dict
        Mapping of display mutation labels to their biological effects.
    """
    muts = {}

    for segment in ('PB2', 'PB1', 'PA', 'HA', 'NP', 'NA', 'MP', 'NS'):
        segment_muts = flags['Sample'].get(f'{segment}_muts', [])
        if segment_muts:
            for mut in segment_muts:
                if mut in muts_loci_meaning:
                    muts[muts_loci_meaning[mut][0]] = muts_loci_meaning[mut][1]

    return muts

def _derive_batch_artifact_names(batch_summary_fp: str) -> tuple[str, str]:
    """
    Derive output artifact names from a batch summary filename.

    Example
    -------
    batch_summary_batch_root_20260420_153012.tsv
      -> batch_index_batch_root_20260420_153012.html
      -> batch_reports_batch_root_20260420_153012.zip
    """
    base = os.path.splitext(os.path.basename(batch_summary_fp))[0]

    if base.startswith("batch_summary_"):
        suffix = base[len("batch_summary_"):]
    else:
        suffix = base

    index_name = f"batch_index_{suffix}.html"
    zip_name = f"batch_reports_{suffix}.zip"
    return index_name, zip_name


def _get_row_reports_dir_rel(row: dict[str, str]) -> str:
    """Resolve the per-sample reports directory relative to the batch root."""
    reports_dir_rel = str(row.get("reports_dir_rel", "")).strip()
    if reports_dir_rel:
        return reports_dir_rel

    sample_dir = str(row.get("sample_dir", "")).strip()
    if sample_dir:
        return sample_dir

    return ""


def _build_batch_row_paths(batch_root: str, row: dict[str, str]) -> tuple[str, str]:
    """Return the absolute and relative reports directory for a batch row."""
    reports_dir_rel = _get_row_reports_dir_rel(row)
    reports_dir_abs = os.path.join(batch_root, reports_dir_rel) if reports_dir_rel else batch_root
    return reports_dir_abs, reports_dir_rel


def _join_rel_path(*parts: str) -> str:
    """Join relative path fragments using forward slashes for HTML and ZIP output."""
    clean_parts = [part for part in parts if part]
    if not clean_parts:
        return ""
    return os.path.join(*clean_parts).replace(os.sep, "/")


def _html_link(href: str, label: str) -> str:
    """Render an anchor that always opens in a new tab."""
    return (
        f'<a href="{escape(str(href).strip())}" target="_blank" '
        f'rel="noopener noreferrer">{escape(str(label).strip())}</a>'
    )


_GITHUB_REPOSITORY_LINK = _html_link(
    'https://github.com/jopereira88/AFluID',
    'AFluID GitHub repository',
)


def _representative_link(accession: Any) -> str:
    """Render a representative accession as an NCBI Nuccore link."""
    accession_text = str(accession).strip()
    if not accession_text or accession_text == 'Not Available':
        return 'Not Available'
    return _html_link(
        f'https://www.ncbi.nlm.nih.gov/nuccore/{accession_text}',
        accession_text,
    )


def _segment_sort_key(segment_name: Any) -> tuple[int, str]:
    """Return a stable sort key for canonical IAV segment ordering."""
    segment_text = str(segment_name).strip()
    return (segment_to_number_iav.get(segment_text, 999), segment_text)


_MUTATION_LOCUS_ORDER = (
    'PB2',
    'PB1',
    'PB1-F2',
    'PA',
    'PA-X',
    'HA',
    'NP',
    'NA',
    'M1',
    'M2',
)


def _mutation_locus_sort_key(locus_name: Any) -> tuple[int, str]:
    """Return a stable sort key for report mutation loci."""
    locus_text = str(locus_name).strip()
    try:
        order_idx = _MUTATION_LOCUS_ORDER.index(locus_text)
    except ValueError:
        order_idx = len(_MUTATION_LOCUS_ORDER)
    return (order_idx, locus_text)


def _mutation_label_sort_key(label: Any) -> tuple[tuple[int, str], str]:
    """Sort mutation display labels by locus first, then full label."""
    label_text = str(label).strip()
    locus_text = label_text.split(':', 1)[0] if ':' in label_text else label_text
    return (_mutation_locus_sort_key(locus_text), label_text)


def order_mutation_labels(labels: Iterable[Any]) -> list[str]:
    """Return unique mutation display labels in report order."""
    cleaned = []
    seen = set()

    for label in labels:
        label_text = str(label).strip()
        if not label_text or label_text in seen:
            continue
        seen.add(label_text)
        cleaned.append(label_text)

    return sorted(cleaned, key=_mutation_label_sort_key)


def _display_mutation_label(segment: str, mutation: str, mut_lookup: dict[str, tuple[str, str]]) -> str:
    """Return the locus-prefixed mutation label used in the mutations tab."""
    mutation_text = str(mutation).strip()
    if not mutation_text:
        return ''
    if mutation_text in mut_lookup:
        return mut_lookup[mutation_text][0]
    return f'{segment}:{mutation_text}'


def _tool_version_placeholder_map() -> dict[str, str]:
    """Return the supported per-tool HTML placeholder tags."""
    return {
        '<!--TOOL_VERSION_BLASTN-->': 'blastn',
        '<!--TOOL_VERSION_CD_HIT_EST_2D-->': 'cd-hit-est-2d',
        '<!--TOOL_VERSION_FLUMUT-->': 'flumut',
        '<!--TOOL_VERSION_FLUMUTDB-->': 'FluMutDB',
        '<!--TOOL_VERSION_GENIN2-->': 'genin2',
        '<!--TOOL_VERSION_NEXTCLADE-->': 'nextclade',
    }


def _tool_version_text(value: Any) -> str:
    """Normalize a tool-version value into display text."""
    if isinstance(value, dict):
        version = str(value.get('version', '')).strip()
        build = str(value.get('build', '')).strip()
        detail = ' '.join(part for part in (version, build) if part)
    else:
        detail = str(value).strip()
    return detail or 'Unknown'


def _apply_named_tool_version_placeholders(
    html_text: str,
    tool_versions: Optional[dict[str, Any]] = None,
) -> str:
    """Replace all supported per-tool placeholders in the supplied HTML."""
    if not tool_versions:
        tool_versions = {}

    for placeholder_tag, tool_key in _tool_version_placeholder_map().items():
        replacement = escape(_tool_version_text(tool_versions.get(tool_key, '')))
        html_text = html_text.replace(placeholder_tag, replacement)

    return html_text


def apply_tool_version_placeholder(
    html_text: str,
    tool_versions: Optional[dict[str, Any]] = None,
    placeholder_tag: Optional[str] = None,
) -> str:
    """Replace a caller-provided placeholder tag with tool version text.

    This helper supports both a caller-provided aggregate placeholder and the
    built-in per-tool placeholders returned by ``_tool_version_placeholder_map()``.
    To use individual tool tags in a template, add markers such as
    ``<!--TOOL_VERSION_BLASTN-->`` or ``<!--TOOL_VERSION_NEXTCLADE-->``.

    Parameters
    ----------
    html_text : str
        HTML content to transform.
    tool_versions : dict, optional
        Mapping such as ``flags['Tool Versions']``. Values may be strings or
        dictionaries containing fields like ``version`` and ``build``.
    placeholder_tag : str, optional
        Exact aggregate marker text to replace with the full multi-line tool
        version block. When omitted, only the per-tool placeholders are
        processed.

    Returns
    -------
    str
        Updated HTML when the placeholder is present, otherwise the original
        HTML.
    """
    html_text = _apply_named_tool_version_placeholders(html_text, tool_versions)

    if not placeholder_tag or placeholder_tag not in html_text:
        return html_text

    if not tool_versions:
        return html_text.replace(placeholder_tag, '')

    lines = []
    for tool_name in sorted(tool_versions):
        value = tool_versions[tool_name]
        detail = _tool_version_text(value)
        lines.append(f'{escape(str(tool_name))}: {escape(detail)}')

    replacement = '<br>'.join(lines)
    return html_text.replace(placeholder_tag, replacement)


def create_batch_index(reports_p: str, ok_rows: list[dict[str, str]], all_rows: Optional[list[dict[str, str]]] = None, output_name: str = "batch_index.html") -> str:
    """
    Create an HTML landing page for batch outputs.

    Parameters
    ----------
    reports_p : str
        Batch root containing generated report files.
    ok_rows : list[dict]
        Summary rows for successful samples.
    all_rows : list[dict], optional
        Full batch manifest rows, including failures.
    output_name : str, default "batch_index.html"
        Filename for the generated HTML index.

    Returns
    -------
    str
        Path to the generated HTML file.
    """
    out_fp = os.path.join(reports_p, output_name)

    html = []
    html.append("<!DOCTYPE html>")
    html.append('<html lang="en">')
    html.append("<head>")
    html.append('<meta charset="utf-8">')
    html.append('<meta name="viewport" content="width=device-width, initial-scale=1">')
    html.append("<title>AFluID Batch Index</title>")
    html.append("""
    <style>
      body{font-family:system-ui,-apple-system,Segoe UI,Roboto,Arial,sans-serif;margin:0;padding:2rem;background:#fff;color:#111}
      .container{max-width:1200px;margin:0 auto}
      h1{margin-top:0}
      table{border-collapse:collapse;width:100%}
      th,td{border:1px solid #ddd;padding:.6rem;text-align:left;vertical-align:top}
      th{background:#f5f5f5}
      .muted{color:#666}
      .missing{color:#999}
      .failed{color:#b91c1c}
      a{text-decoration:none}
      a:hover{text-decoration:underline}
      code{background:#f3f4f6;padding:.1rem .3rem;border-radius:.25rem}
    </style>
    """)
    html.append("</head>")
    html.append("<body>")
    html.append('<div class="container">')
    html.append("<h1>AFluID Batch Index</h1>")

    if all_rows is not None:
        total = len(all_rows)
        oks = len(ok_rows)
        fails = total - oks
        html.append(f'<p class="muted">Samples in manifest: {total} | Successful: {oks} | Failed: {fails}</p>')

    html.append("<table>")
    html.append("""
    <thead>
      <tr>
        <th>Input file</th>
        <th>File tag</th>
        <th>Status</th>
        <th>Final report</th>
        <th>Extras</th>
        <th>Error</th>
      </tr>
    </thead>
    <tbody>
    """)

    ok_lookup = {(str(r.get("input_file", "")).strip(), str(r.get("file_tag", "")).strip()) for r in ok_rows}

    rows_to_render = all_rows if all_rows is not None else ok_rows

    for row in rows_to_render:
        input_file = str(row.get("input_file", "")).strip()
        file_tag = str(row.get("file_tag", "")).strip()
        status = str(row.get("status", "")).strip().lower()
        error = str(row.get("error", "")).strip()

        final_report = f"{file_tag}_final_report.html"
        id_report = f"{file_tag}_ID_Report.txt"
        flags_json = f"{file_tag}_flags.json"
        reports_dir_abs, reports_dir_rel = _build_batch_row_paths(reports_p, row)

        is_ok = (input_file, file_tag) in ok_lookup and status == "ok"

        final_report_rel = _join_rel_path(reports_dir_rel, final_report)
        id_report_rel = _join_rel_path(reports_dir_rel, id_report)
        flags_rel = _join_rel_path(reports_dir_rel, flags_json)

        if is_ok and os.path.exists(os.path.join(reports_dir_abs, final_report)):
            final_report_link = _html_link(final_report_rel, final_report)
        else:
            final_report_link = '<span class="missing">missing</span>'

        extras = []
        if is_ok and os.path.exists(os.path.join(reports_dir_abs, id_report)):
            extras.append(_html_link(id_report_rel, id_report))
        if is_ok and os.path.exists(os.path.join(reports_dir_abs, flags_json)):
            extras.append(_html_link(flags_rel, flags_json))

        if file_tag and os.path.isdir(reports_dir_abs):
            for fname in sorted(os.listdir(reports_dir_abs)):
                if fname.startswith(f"{file_tag}_") and fname not in {final_report, id_report, flags_json}:
                    full = os.path.join(reports_dir_abs, fname)
                    if os.path.isfile(full):
                        fname_rel = _join_rel_path(reports_dir_rel, fname)
                        extras.append(_html_link(fname_rel, fname))

        extras_html = "<br>".join(extras) if extras else '<span class="missing">none</span>'
        status_html = escape(status) if status != "failed" else f'<span class="failed">{escape(status)}</span>'
        error_html = escape(error) if error else ""

        html.append(f"""
        <tr>
          <td>{escape(input_file)}</td>
          <td><code>{escape(file_tag)}</code></td>
          <td>{status_html}</td>
          <td>{final_report_link}</td>
          <td>{extras_html}</td>
          <td>{error_html}</td>
        </tr>
        """)

    html.append("</tbody></table>")
    html.append("</div></body></html>")

    with open(out_fp, "w", encoding="utf-8") as f:
        f.write("\n".join(html))

    return out_fp


def create_batch_zip(reports_p: str, ok_rows: list[dict[str, str]], zip_name: str = "batch_reports.zip", extra_files: Optional[list[str]] = None) -> Optional[str]:
    """
    Bundle batch report artifacts into a ZIP archive.

    Parameters
    ----------
    reports_p : str
        Batch root containing generated report files.
    ok_rows : list[dict]
        Summary rows for successful samples.
    zip_name : str, default "batch_reports.zip"
        Filename for the generated archive.
    extra_files : list[str], optional
        Additional filenames inside ``reports_p`` to include.

    Returns
    -------
    str or None
        Path to the generated archive, or None when no files were eligible.
    """
    if extra_files is None:
        extra_files = []

    zip_fp = os.path.join(reports_p, zip_name)
    files_to_add: dict[str, str] = {}

    for row in ok_rows:
        file_tag = str(row.get("file_tag", "")).strip()
        if not file_tag:
            continue

        reports_dir_abs, reports_dir_rel = _build_batch_row_paths(reports_p, row)
        if not os.path.isdir(reports_dir_abs):
            continue

        for fname in os.listdir(reports_dir_abs):
            if fname.startswith(f"{file_tag}_"):
                full = os.path.join(reports_dir_abs, fname)
                if os.path.isfile(full):
                    arcname = _join_rel_path(reports_dir_rel, fname)
                    files_to_add[arcname] = full

    for fname in extra_files:
        full = os.path.join(reports_p, fname)
        if os.path.isfile(full):
            files_to_add[fname] = full

    if not files_to_add:
        return None

    with zipfile.ZipFile(zip_fp, "w", compression=zipfile.ZIP_DEFLATED) as zf:
        for arcname in sorted(files_to_add):
            zf.write(files_to_add[arcname], arcname=arcname)

    return zip_fp


def maybe_create_batch_artifacts(reports_p: str, batch_summary_fp: str, extra_files: Optional[list[str]] = None) -> bool:
    """
    Create batch index + zip only if the provided batch summary exists, is valid,
    and contains at least one successful sample.

    Parameters
    ----------
    reports_p : str
        Reports directory.
    batch_summary_fp : str
        Path to a specific batch summary TSV.

    Returns
    -------
    bool
        True if artifacts were generated, False otherwise.
    """
    if not batch_summary_fp:
        return False

    if extra_files is None:
        extra_files = []

    if not os.path.exists(batch_summary_fp):
        return False

    rows = []
    try:
        with open(batch_summary_fp, "r", encoding="utf-8") as f:
            reader = csv.DictReader(f, delimiter="\t")

            required = {"input_file", "file_tag", "status", "error"}

            if reader.fieldnames is None:
                return False

            if not required.issubset(set(reader.fieldnames)):
                return False

            for row in reader:
                rows.append(row)
    except Exception:
        return False

    if not rows:
        return False

    ok_rows = [
        row for row in rows
        if str(row.get("status", "")).strip().lower() == "ok"
           and str(row.get("file_tag", "")).strip() != ""
    ]

    if not ok_rows:
        return False

    index_name, zip_name = _derive_batch_artifact_names(batch_summary_fp)

    create_batch_index(
        reports_p=reports_p,
        ok_rows=ok_rows,
        all_rows=rows,
        output_name=index_name
    )

    create_batch_zip(
        reports_p=reports_p,
        ok_rows=ok_rows,
        zip_name=zip_name,
        extra_files=[os.path.basename(batch_summary_fp), index_name, *extra_files]
    )

    return True

def create_sample_card(flags: dict) -> str:
    """
    Create the summary card for the Sample tab.

    Parameters
    ----------
    flags : dict
        Pipeline state dictionary containing inferred sample metadata and tool
        statuses.

    Returns
    -------
    str
        HTML fragment for the Sample summary card.
    """
    GENOTYPE=flags['Sample']['Genotype']
    H=flags['Sample']['H_gen']
    N=flags['Sample']['N_gen']
    CLADE=flags['Sample']['clade'] if flags['Sample']['clade'] else 'N/A'
    C_BLAST=flags['Master']['C_BLAST']
    L_BLAST=flags['Master']['L_BLAST']
    NEXTCLADE=flags['Master']['nextclade']
    FLUMUT=flags['Master']['flumut']
    GENIN=flags['Master']['genin']
    SINGLE=flags['Master']['single']
    card=f"""
    <div class="card" id="sample-card">
  <p>
    <strong>Genotype:</strong> {GENOTYPE}
    &nbsp;·&nbsp; <strong>H:</strong> {H}
    &nbsp;·&nbsp; <strong>N:</strong> {N}
    &nbsp;·&nbsp; <strong>Clade:</strong> {CLADE}
  </p>

  <div class="chips" id="pipeline-chips">
    <!-- Preencha com chips de estado dos módulos da pipeline -->
    <!-- Exemplo (repita/ajuste conforme preciso): -->
    <span class="chip">C_BLAST: <b>{C_BLAST}</b></span>
    <span class="chip">L_BLAST: <b>{L_BLAST}</b></span>
    <span class="chip">Nextclade: <b>{NEXTCLADE}</b></span>
    <span class="chip">FluMut: <b>{FLUMUT}</b></span>
    <span class="chip">GenIn: <b>{GENIN}</b></span>
    <span class="chip">Single: <b>{SINGLE}</b></span>
  </div>
</div>"""
    return card

def create_sample_table(
    df_segments: str,
    flags: dict,
    seg_lens: dict,
    mappings: dict[str, str],
    muts_loci_meaning: dict[str, tuple[str, str]],
) -> str:
    """
    Create the per-segment summary table for the Sample tab.

    Parameters
    ----------
    df_segments : str
        Path to the ID report TSV used to determine row ordering.
    flags : dict
        Pipeline state dictionary containing segment assignments, references,
        sequence lengths, and mutation hits.
    seg_lens : dict
        Expected minimum and maximum length thresholds per segment.
    mappings : dict
        Mapping from original sequence headers to internal identifiers.

    Returns
    -------
    str
        HTML table fragment.
    """
    len_thres = {
        'PB2': (seg_lens['PB2_min'], seg_lens['PB2_max']),
        'PB1': (seg_lens['PB1_min'], seg_lens['PB1_max']),
        'PA': (seg_lens['PA_min'], seg_lens['PA_max']),
        'HA': (seg_lens['HA_min'], seg_lens['HA_max']),
        'NP': (seg_lens['NP_min'], seg_lens['NP_max']),
        'NA': (seg_lens['NA_min'], seg_lens['NA_max']),
        'MP': (seg_lens['MP_min'], seg_lens['MP_max']),
        'NS': (seg_lens['NS_min'], seg_lens['NS_max'])
    }

    seq_names = {}
    for segment in iav_segments:
        for name in flags['Sample'].get(segment, []):
            if name in mappings:
                seq_names[mappings[name]] = segment

    df, error = _safe_read_table(df_segments, required_cols=['SAMPLE_NAME'])
    if error:
        return _empty_card(f'Sample segment summary unavailable. {error}')

    table = '''
    <table class="tabela" id="sample-segments">
    <thead>
    <tr>
        <th>Segment</th>
        <th>Length</th>
        <th>Reference</th>
        <th>Mutations</th>
        <th>Note</th>
    </tr>
    </thead>
    <tbody>
    '''

    for segment in iav_segments:
        if segment in len_thres:
            length = flags['Sample'].get(f'{segment}_len', 0)
            ref_list = flags['Sample'].get(f'{segment}_ref', [])
            reference = _representative_link(ref_list[0] if ref_list else "Not Available")
            muts_list = flags['Sample'].get(f'{segment}_muts', [])
            mutations = ', '.join(
                _display_mutation_label(segment, mut, muts_loci_meaning)
                for mut in muts_list
                if str(mut).strip()
            ) if muts_list else "None"

            min_len, max_len = len_thres[segment]
            if length < min_len:
                length_class = "len-below"
                note = 'Length outside of expected range'
            elif length > max_len:
                length_class = "len-above"
                note = 'Length outside of expected range'
            else:
                length_class = "len-ok"
                note = ''

            length_html = f'<span class="{length_class}">{length}</span>'
        else:
            length_html = "0"
            reference = "Not Available"
            mutations = "None"
            note = "Segment not found"

        table += f"""
          <tr>
            <td>{segment}</td>
            <td>{length_html}</td>
            <td>{reference}</td>
            <td>{mutations}</td>
            <td>{note}</td>
          </tr>
        """

    table += """
      </tbody>
    </table>
    """

    return table


def create_segment_background_table(
    df_segments: str,
    flags: dict,
    seg_lens: dict,
    mappings: dict[str, str],
    geo_dict: dict[str, tuple[str, str]],
    taxa_dict: dict[str, list[str]],
) -> str:
    """
    Create the Segment Background table for single-sample reports.

    Parameters
    ----------
    df_segments : str
        Path to the ID report TSV.
    flags : dict
        Pipeline state dictionary containing segment assignments, references,
        and sequence lengths.
    seg_lens : dict
        Expected minimum and maximum length thresholds per segment.
    mappings : dict
        Mapping from original sequence headers to internal identifiers.

    Returns
    -------
    str
        HTML table fragment.

    Notes
    -----
    Cluster- and C-BLAST-derived rows are normalized from dict-like metadata
    distributions, while L-BLAST-derived rows are rendered from flat string
    metadata.
    """
    len_thres = {
        'PB2': (seg_lens['PB2_min'], seg_lens['PB2_max']),
        'PB1': (seg_lens['PB1_min'], seg_lens['PB1_max']),
        'PA': (seg_lens['PA_min'], seg_lens['PA_max']),
        'HA': (seg_lens['HA_min'], seg_lens['HA_max']),
        'NP': (seg_lens['NP_min'], seg_lens['NP_max']),
        'NA': (seg_lens['NA_min'], seg_lens['NA_max']),
        'MP': (seg_lens['MP_min'], seg_lens['MP_max']),
        'NS': (seg_lens['NS_min'], seg_lens['NS_max'])
    }

    seq_names = {}
    for segment in iav_segments:
        for name in flags['Sample'].get(segment, []):
            if name in mappings:
                seq_names[mappings[name]] = segment

    df, error = _safe_read_table(df_segments, required_cols=['SAMPLE_NAME'])
    if error:
        return _empty_card(f'Segment background unavailable. {error}')

    for col, default in {
        'ASSIGNED_BY': 'NA',
        '%ID': 'NA',
        'CLUSTER': 'NA',
        'HOST': 'Unknown',
        'GENOTYPE': 'Unknown',
        'COUNTRY': 'Unknown',
        'COLLECTION_DATE': 'Unknown',
    }.items():
        if col not in df.columns:
            df[col] = default

    table = """
    <div class="sb-table-wrap">
    <table class="tabela" id="sb-flat-table">
    <thead>
      <tr>
        <th>Segment</th>
        <th>Length</th>
        <th>Note</th>
        <th>Cluster</th>
        <th>Reference</th>
        <th>%ID</th>
        <th>Assigned_By</th>
        <th>Genotype</th>
        <th>Host</th>
        <th>Country</th>
        <th>Collection_date</th>
      </tr>
    </thead>
    <tbody>"""

    mask = df['ASSIGNED_BY'] != 'L-BLAST'

    df['HOST_CLEAN'] = None
    df['GENOTYPE_CLEAN'] = None
    df['COUNTRY_CLEAN'] = None

    df.loc[mask, 'HOST_CLEAN'] = df.loc[mask, 'HOST'].apply(
        lambda x: clean_host(x, taxa_dict, level='Species', drop_unknown=True)
    )

    df.loc[mask, 'GENOTYPE_CLEAN'] = df.loc[mask, 'GENOTYPE'].apply(
        lambda x: clean_genotype(x)
    )

    df.loc[mask, 'COUNTRY_CLEAN'] = df.loc[mask, 'COUNTRY'].apply(
        lambda x: clean_country(x, geo_dict, mode='country', drop_unknown=True)
    )

    df['SEGMENT_NAME'] = df['SAMPLE_NAME'].apply(lambda value: seq_names.get(value, 'Not Available'))
    df = df.sort_values(
        by='SEGMENT_NAME',
        key=lambda series: series.map(lambda value: _segment_sort_key(value)),
        kind='stable',
    ).reset_index(drop=True)

    for row in range(len(df)):
        segment = df.loc[row, 'SEGMENT_NAME']

        if segment in len_thres:
            length = flags['Sample'].get(f'{segment}_len', 0)
            ref_list = flags['Sample'].get(f'{segment}_ref', [])
            reference = _representative_link(ref_list[0] if ref_list else "Not Available")

            min_len, max_len = len_thres[segment]
            if length < min_len:
                length_class = "len-below"
                note = 'Length outside of expected range'
            elif length > max_len:
                length_class = "len-above"
                note = 'Length outside of expected range'
            else:
                length_class = "len-ok"
                note = ''

            length_html = f'<span class="{length_class}">{length}</span>'
        else:
            length_html = "0"
            reference = "Not Available"
            note = "Segment not found"

        assigned_by = df.loc[row, 'ASSIGNED_BY']
        perc_id = df.loc[row, '%ID']
        cluster = df.loc[row, 'CLUSTER']
        collection_date = _normalize_collection_date(df.loc[row, 'COLLECTION_DATE'])

        if assigned_by != 'L-BLAST':
            host_species = df.loc[row, 'HOST_CLEAN']
            genotype = df.loc[row, 'GENOTYPE_CLEAN']
            geo_countries = df.loc[row, 'COUNTRY_CLEAN']

            host_species = host_species if isinstance(host_species, dict) else {}
            genotype = genotype if isinstance(genotype, dict) else {}
            geo_countries = geo_countries if isinstance(geo_countries, dict) else {}

            host_order = rollup_counts_taxa(host_species, taxa_dict, level='Order', drop_unknown=True)
            host_class = rollup_counts_taxa(host_species, taxa_dict, level='Class', drop_unknown=True)
            geo_regions = rollup_counts_geo(geo_countries, geo_dict, mode='region', drop_unknown=True)
            geo_continents = rollup_counts_geo(geo_countries, geo_dict, mode='continent', drop_unknown=True)

        else:
            host_species = df.loc[row, 'HOST']
            genotype = df.loc[row, 'GENOTYPE']
            geo_countries = df.loc[row, 'COUNTRY']

            host_order = rollup_taxa(host_species, taxa_dict, level='Order')
            host_class = rollup_taxa(host_species, taxa_dict, level='Class')
            geo_regions = rollup_geo(geo_countries, geo_dict, mode='region')
            geo_continents = rollup_geo(geo_countries, geo_dict, mode='continent')

        if note != "":
            note += "; Segment not found" if segment == "Not Available" else ""
        else:
            note = "Segment not found" if segment == "Not Available" else ""

        table += f"""
        <tr>
          <td>{segment}</td>
          <td>{length_html}</td>
          <td>{note}</td>
          <td>{cluster}</td>
          <td>{reference}</td>
          <td>{perc_id}</td>
          <td>{assigned_by}</td>
          <td>{genotype}</td>
          <td>
            <span class="host host--species">{host_species}</span>
            <span class="host host--order">{host_order}</span>
            <span class="host host--class">{host_class}</span>
          </td>
          <td>
            <span class="geo geo--country">{geo_countries}</span>
            <span class="geo geo--region">{geo_regions}</span>
            <span class="geo geo--continent">{geo_continents}</span>
          </td>
          <td>{collection_date}</td>
        </tr>
        """

    table += """
        </tbody>
        </table>
        </div>
    """

    return table


def create_segment_mutations_table(flags: dict, muts_loci_meaning: dict) -> str:
    """
    Create the Segment Mutations table for the final report.

    Parameters
    ----------
    flags : dict
        Pipeline state dictionary containing mutation hits per segment.
    muts_loci_meaning : dict
        Mapping from mutation identifiers to display label and biological
        meaning.

    Returns
    -------
    str
        HTML table fragment.
    """
    muts = convert_muts(flags, muts_loci_meaning)

    table = """
    <table class="tabela" id="mutations-table">
    <thead>
      <tr>
        <th>Mutation</th>
        <th>Biological effects</th>
      </tr>
    </thead>
    <tbody>"""

    for mut in order_mutation_labels(muts):
        meaning = muts[mut]
        table += f"""
        <tr>
          <td>{mut}</td>
          <td>{meaning}</td>
        </tr>
        """

    table += """
        </tbody>
        </table>
    """

    return table

def create_genin_const_table(flags: dict) -> str:
    """
    Create the GenIn2 constellation table for the final report.

    Parameters
    ----------
    flags : dict
        Pipeline state dictionary containing parsed GenIn2 results.

    Returns
    -------
    str
        HTML fragment containing either a results table or an empty-state card.
    """
    genin_const = flags['Sample'].get('Genin_genotypes', [])

    if not genin_const:
        return """
        <div class="card">
          <p class="muted">No GenIn2 data available.</p>
        </div>
        """

    table = """
    <table class="tabela" id="genin-const-table">
      <thead>
        <tr>
          <th>Sample Name</th>
          <th>Genotype</th>
          <th>Sub-genotype</th>
          <th>PB2</th>
          <th>PB1</th>
          <th>PA</th>
          <th>NP</th>
          <th>NA</th>
          <th>MP</th>
          <th>NS</th>
          <th>Notes</th>
        </tr>
      </thead>
      <tbody>
    """

    if isinstance(genin_const, dict):
        # Backward compatibility for previously serialized flags.
        rows = []
        for sample_name, data in genin_const.items():
            row = {
                'Sample Name': sample_name,
                'Genotype': data.get('Genotype', ''),
                'Sub-genotype': data.get('Sub-genotype', ''),
                'PB2': '',
                'PB1': '',
                'PA': '',
                'NP': '',
                'NA': '',
                'MP': '',
                'NS': '',
                'Notes': '',
            }
            segment = data.get('segment', '')
            value = data.get('value', '')
            if segment in ('PB2', 'PB1', 'PA', 'NP', 'NA', 'MP', 'NS'):
                row[segment] = value
            rows.append(row)
    else:
        rows = genin_const

    for data in rows:
        sample_name = data.get('Sample Name', '')
        genotype = data.get('Genotype', '')
        sub_genotype = data.get('Sub-genotype', '')
        pb2 = data.get('PB2', '')
        pb1 = data.get('PB1', '')
        pa = data.get('PA', '')
        np = data.get('NP', '')
        na = data.get('NA', '')
        mp = data.get('MP', '')
        ns = data.get('NS', '')
        notes = data.get('Notes', '')

        table += f"""
        <tr>
          <td>{sample_name}</td>
          <td>{genotype}</td>
          <td>{sub_genotype}</td>
          <td>{pb2}</td>
          <td>{pb1}</td>
          <td>{pa}</td>
          <td>{np}</td>
          <td>{na}</td>
          <td>{mp}</td>
          <td>{ns}</td>
          <td>{notes}</td>
        </tr>
        """

    table += """
      </tbody>
    </table>
    """

    return table


html_skeleton = r"""<!DOCTYPE html>
<html lang="EN">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title><!--FILENAME--> Final Report </title>

  <style>
    :root{
      --bg:#fff; --fg:#111827; --muted:#6b7280;
      --accent:#2563eb; --accent-weak:#e5edff;
      --border:#e5e7eb; --card:#f9fafb;
      --sans: system-ui,-apple-system, Segoe UI, Roboto, Noto Sans, Ubuntu, Cantarell, Arial;
    }
    html,body{height:100%}
    body{margin:0;font-family:var(--sans);color:var(--fg);background:var(--bg);line-height:1.55}

    header{padding:1rem;border-bottom:1px solid var(--border)}
    .container{max-width:1100px;margin:0 auto;padding:1rem}

    .tabs{border:1px solid var(--border); border-radius:.75rem; background:#fff; overflow:hidden}
    .tab-controls{display:flex; gap:.5rem; padding:.5rem; border-bottom:1px solid var(--border); background:var(--card); flex-wrap:wrap}
    .tab-controls label{cursor:pointer; padding:.5rem .75rem; border:1px solid var(--border); border-radius:.5rem; background:#fff; font-weight:600; font-size:.95rem}
    input[name="tab"]{position:absolute; left:-9999px}

    input#tab1:checked ~ .tab-controls label[for="tab1"],
    input#tab2:checked ~ .tab-controls label[for="tab2"],
    input#tab3:checked ~ .tab-controls label[for="tab3"]
    <!--TAB4_ACTIVE_CSS-->{
      background:var(--accent-weak); border-color:var(--accent); color:var(--accent)
    }

    .tab-panel{display:none;padding:1rem}
    input#tab1:checked ~ .tab-panels #panel1,
    input#tab2:checked ~ .tab-panels #panel2,
    input#tab3:checked ~ .tab-panels #panel3
    <!--TAB4_DISPLAY_CSS-->{display:block}

    .card{background:#fff;border:1px solid var(--border);border-radius:.75rem;padding:1rem}
    table{border-collapse:collapse;width:100%}
    th, td{border:1px solid var(--border);padding:.5rem;text-align:left;font-size:.95rem;vertical-align:top}
    th{background:#f3f4f6}
    .muted{color:var(--muted)}

    .chips .chip{display:inline-block;margin:.15rem .3rem;padding:.15rem .45rem;border:1px solid var(--border);border-radius:.4rem;background:#f9fafb}
    .len-ok{background:#ecfdf5; border:1px solid #bbf7d0; padding:.1rem .3rem; border-radius:.3rem;}
    .len-below{background:#fef2f2; border:1px solid #fecaca; padding:.1rem .3rem; border-radius:.3rem;}
    .len-above{background:#fff7ed; border:1px solid #fed7aa; padding:.1rem .3rem; border-radius:.3rem;}

    .vh{
      position:absolute; width:1px; height:1px;
      margin:-1px; padding:0; border:0; overflow:hidden;
      clip:rect(0 0 0 0); clip-path:inset(50%); white-space:nowrap;
    }

    #panel2 .sb-card{display:grid; grid-template-rows:auto 1fr; gap:.25rem}
    #panel2 .toolbars{display:flex; gap:.75rem; justify-content:flex-end; align-items:center; margin:0 0 .25rem 0; padding-right:.5rem}
    .controls{display:inline-flex; gap:.25rem; align-items:center}
    .controls .group-label{color:var(--muted); font-weight:600; font-size:.85rem; margin-right:.25rem}
    .controls label{
      display:inline-block; white-space:nowrap;
      padding:.15rem .5rem; font-size:.8rem; line-height:1;
      border:1px solid var(--border); border-radius:.4rem; background:#fff; cursor:pointer;
    }

    #panel2 .sb-table-wrap .host{display:none}
    #panel2 .sb-table-wrap .geo{display:none}

    #panel2 #hv-species:checked ~ .sb-card .sb-table-wrap .host--species{display:inline}
    #panel2 #hv-order:checked   ~ .sb-card .sb-table-wrap .host--order  {display:inline}
    #panel2 #hv-class:checked   ~ .sb-card .sb-table-wrap .host--class  {display:inline}

    #panel2 #gv-country:checked   ~ .sb-card .sb-table-wrap .geo--country  {display:inline}
    #panel2 #gv-region:checked    ~ .sb-card .sb-table-wrap .geo--region   {display:inline}
    #panel2 #gv-continent:checked ~ .sb-card .sb-table-wrap .geo--continent{display:inline}

    #panel2 #hv-species:checked ~ .sb-card .controls-host label[for="hv-species"],
    #panel2 #hv-order:checked   ~ .sb-card .controls-host label[for="hv-order"],
    #panel2 #hv-class:checked   ~ .sb-card .controls-host label[for="hv-class"],
    #panel2 #gv-country:checked   ~ .sb-card .controls-geo label[for="gv-country"],
    #panel2 #gv-region:checked    ~ .sb-card .controls-geo label[for="gv-region"],
    #panel2 #gv-continent:checked ~ .sb-card .controls-geo label[for="gv-continent"]{
      background:var(--accent-weak); color:var(--accent); border-color:var(--accent);
    }

    @media print{
      .tab-controls{display:none}
      .tab-panel{display:block}
      .tabs{border:none}
    }
  </style>
</head>
<body>
  <header class="container">
    <h1><!--FILENAME--> Final Report</h1>
    <p class="muted"> HTML interactive report for Single-Sample mode</p>
  </header>

  <main class="container">
    <section class="tabs card">
      <input type="radio" id="tab1" name="tab" checked>
      <input type="radio" id="tab2" name="tab">
      <input type="radio" id="tab3" name="tab">
      <!--TAB4_INPUT-->

      <div class="tab-controls">
        <label for="tab1">Sample</label>
        <label for="tab2">Segment Background</label>
        <label for="tab3">Segment Mutations</label>
        <!--TAB4_LABEL-->
      </div>

      <div class="tab-panels">
        <article id="panel1" class="tab-panel">
          <h2>Sample Information</h2>
          <p class="muted">
            This page displays summarised information from cd-hit (<!--TOOL_VERSION_CD_HIT_EST_2D-->)
            and BLAST (<!--TOOL_VERSION_BLASTN-->). Clade assignments, when available, were
            performed with Nextclade (<!--TOOL_VERSION_NEXTCLADE-->), and mutation analysis was
            provided by FluMut (<!--TOOL_VERSION_FLUMUT-->) using FluMutDB (<!--TOOL_VERSION_FLUMUTDB-->).
          </p>
          <!--SAMPLE_CARD-->
          <!--SAMPLE_TABLE-->
        </article>

        <article id="panel2" class="tab-panel">
          <h2>Segment Background</h2>
          <p class="muted">
            Placements of sequences within the reported diversity. To examine cluster composition,
            please refer to the cluster composition TSV file.
          </p>

          <input class="vh" type="radio" id="hv-species" name="hostview" checked>
          <input class="vh" type="radio" id="hv-order"   name="hostview">
          <input class="vh" type="radio" id="hv-class"   name="hostview">

          <input class="vh" type="radio" id="gv-country"   name="geoview" checked>
          <input class="vh" type="radio" id="gv-region"    name="geoview">
          <input class="vh" type="radio" id="gv-continent" name="geoview">

          <div class="sb-card card">
            <div class="toolbars">
              <div class="controls controls-host" aria-label="Host view">
                <span class="group-label">Host:</span>
                <label for="hv-species">Species</label>
                <label for="hv-order">Order</label>
                <label for="hv-class">Class</label>
              </div>
              <div class="controls controls-geo" aria-label="Geography view">
                <span class="group-label">Geography:</span>
                <label for="gv-country">Country</label>
                <label for="gv-region">Region</label>
                <label for="gv-continent">Continent</label>
              </div>
            </div>

            <div class="sb-table-wrap">
              <!--SEGMENT_BACKGROUND_TABLE-->
            </div>
          </div>
        </article>

        <article id="panel3" class="tab-panel">
          <h2>Segment Mutations</h2>
          <p class="muted">
            Obtained using FluMut (<!--TOOL_VERSION_FLUMUT-->) running FluMutDB (<!--TOOL_VERSION_FLUMUTDB-->).
            These mutations were further filtered for mutations of importance as cited by Alvarez et al. 2025 - 34 key
            mutations and 6 mutations of secondary importance - and Mohaptra et al. 2023 - for all antiviral resistance related mutations. 
            For a list of all mutations found, please refer to the final FluMut spreadsheet report.
            For a list of the mutations of interest parsed in this pipeline, please refer to the """ + _GITHUB_REPOSITORY_LINK + r""".
          </p>
          <!--SEGMENT_MUTATIONS_TABLE-->
        </article>

        <!--TAB4_PANEL-->
      </div>
    </section>
  </main>
</body>
</html>
"""
def to_html_report(
    sample_card: str,
    sample_table: str,
    seg_back_table: str,
    mut_table: str,
    output_file_path: str,
    filename: str,
    genin2: bool = False,
    genin2table: str = '',
    tool_versions: Optional[dict[str, Any]] = None,
    tool_version_placeholder_tag: Optional[str] = None,
    html_skeleton: str = html_skeleton
) -> None:
    """
    Render the single-sample final HTML report.

    Parameters
    ----------
    sample_card : str
        HTML fragment for the Sample summary card.
    sample_table : str
        HTML fragment for the Sample tab table.
    seg_back_table : str
        HTML fragment for the Segment Background tab.
    mut_table : str
        HTML fragment for the Segment Mutations tab.
    output_file_path : str
        Path to the HTML file to write.
    filename : str
        Report name displayed in the HTML title and header.
    genin2 : bool, default False
        Whether the GenIn2 tab should be included.
    genin2table : str, default ''
        HTML fragment for the GenIn2 tab.
    tool_versions : dict, optional
        Optional tool/version metadata, typically ``flags['Tool Versions']``.
        When combined with ``tool_version_placeholder_tag``, the helper
        ``apply_tool_version_placeholder()`` replaces that caller-provided tag.
    tool_version_placeholder_tag : str, optional
        Exact placeholder string to replace with rendered tool versions. No
        built-in placeholder is injected into the skeleton automatically.
    html_skeleton : str, default ``html_skeleton``
        Template HTML document containing replacement markers.

    Returns
    -------
    None
    """
    html_content = html_skeleton
    html_content = html_content.replace('<!--FILENAME-->', filename.replace("_final_report", ""))
    html_content = html_content.replace('<!--SAMPLE_CARD-->', sample_card)
    html_content = html_content.replace('<!--SAMPLE_TABLE-->', sample_table)
    html_content = html_content.replace('<!--SEGMENT_BACKGROUND_TABLE-->', seg_back_table)
    html_content = html_content.replace('<!--SEGMENT_MUTATIONS_TABLE-->', mut_table)

    if genin2:
        tab4_input = '<input type="radio" id="tab4" name="tab">'
        tab4_label = '<label for="tab4">Genin2</label>'
        tab4_active_css = ',\n    input#tab4:checked ~ .tab-controls label[for="tab4"]'
        tab4_display_css = ',\n    input#tab4:checked ~ .tab-panels #panel4'
        tab4_panel = """
        <article id="panel4" class="tab-panel">
          <h2>Genin2</h2>
          <p class="muted">
            GenIn2 performs genotype constellation assignment. This tab reports the
            constellation assignment for this sample using GenIn2 (<!--TOOL_VERSION_GENIN2-->).
          </p>
          <!--GENIN2_TABLE-->
        </article>
        """
    else:
        tab4_input = ''
        tab4_label = ''
        tab4_active_css = ''
        tab4_display_css = ''
        tab4_panel = ''

    html_content = html_content.replace('<!--TAB4_INPUT-->', tab4_input)
    html_content = html_content.replace('<!--TAB4_LABEL-->', tab4_label)
    html_content = html_content.replace('<!--TAB4_ACTIVE_CSS-->', tab4_active_css)
    html_content = html_content.replace('<!--TAB4_DISPLAY_CSS-->', tab4_display_css)
    html_content = html_content.replace('<!--TAB4_PANEL-->', tab4_panel)

    if genin2:
        html_content = html_content.replace('<!--GENIN2_TABLE-->', genin2table)

    html_content = apply_tool_version_placeholder(
        html_content,
        tool_versions=tool_versions,
        placeholder_tag=tool_version_placeholder_tag,
    )

    with open(output_file_path, 'w', encoding='utf-8') as f:
        f.write(html_content)

###### MULTI SAMPLE

html_skeleton_multi = r"""
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width,initial-scale=1" />
  <title><!--FILENAME--> Final Report</title>
  <style>
    :root{
      --bg:#f7f8fa;
      --card:#ffffff;
      --text:#1f2937;
      --muted:#6b7280;
      --line:#e5e7eb;
      --accent:#2563eb;
      --accent-soft:#dbeafe;
      --ok:#065f46;
      --ok-bg:#d1fae5;
      --warn:#92400e;
      --warn-bg:#fef3c7;
      --bad:#991b1b;
      --bad-bg:#fee2e2;
    }

    *{box-sizing:border-box}
    html,body{margin:0;padding:0;background:var(--bg);color:var(--text);font-family:system-ui,-apple-system,Segoe UI,Roboto,Arial,sans-serif}
    body{line-height:1.45}

    .container{max-width:1400px;margin:0 auto;padding:24px}
    .card{
      background:var(--card);
      border:1px solid var(--line);
      border-radius:16px;
      padding:18px;
      box-shadow:0 1px 2px rgba(0,0,0,.04);
    }

    h1,h2,h3{margin:0 0 12px 0}
    h1{font-size:1.9rem}
    h2{font-size:1.2rem}
    p{margin:0 0 10px 0}
    .muted{color:var(--muted)}

    .tabs{margin-top:18px}
    .tabs input[type="radio"]{display:none}

    .tab-controls{
      display:flex;
      gap:8px;
      flex-wrap:wrap;
      margin-bottom:18px;
      border-bottom:1px solid var(--line);
      padding-bottom:12px;
    }

    .tab-controls label{
      display:inline-block;
      padding:10px 14px;
      border:1px solid var(--line);
      border-radius:12px;
      background:#fff;
      cursor:pointer;
      font-weight:600;
      transition:.15s ease-in-out;
    }

    .tab-controls label:hover{
      border-color:var(--accent);
      background:#f8fbff;
    }

    #tab1:checked ~ .tab-controls label[for="tab1"],
    #tab2:checked ~ .tab-controls label[for="tab2"],
    #tab3:checked ~ .tab-controls label[for="tab3"]
    <!--TAB4_ACTIVE_CSS-->
    {
      border-color:var(--accent);
      background:var(--accent-soft);
      color:#0f3ea8;
    }

    .tab-panels .tab-panel{display:none}
    #tab1:checked ~ .tab-panels #panel1{display:block}
    #tab2:checked ~ .tab-panels #panel2{display:block}
    #tab3:checked ~ .tab-panels #panel3{display:block}
    <!--TAB4_DISPLAY_CSS-->

    .toolbar{
      display:flex;
      gap:16px;
      flex-wrap:wrap;
      align-items:center;
      margin:14px 0 16px 0;
      padding:12px 14px;
      border:1px solid var(--line);
      border-radius:12px;
      background:#fbfbfc;
    }

    .toolbar-group{
      display:flex;
      gap:8px;
      align-items:center;
      flex-wrap:wrap;
    }

    .toolbar-group .label{
      font-size:.95rem;
      font-weight:700;
      color:var(--text);
      margin-right:4px;
    }

    .toolbar input[type="radio"]{display:none}

    .toolbar label{
      display:inline-block;
      padding:6px 10px;
      border:1px solid var(--line);
      border-radius:999px;
      cursor:pointer;
      background:#fff;
      font-size:.92rem;
    }

    #host-species:checked ~ .toolbar .host-controls label[for="host-species"],
    #host-order:checked   ~ .toolbar .host-controls label[for="host-order"],
    #host-class:checked   ~ .toolbar .host-controls label[for="host-class"],
    #geo-country:checked  ~ .toolbar .geo-controls label[for="geo-country"],
    #geo-region:checked   ~ .toolbar .geo-controls label[for="geo-region"],
    #geo-continent:checked ~ .toolbar .geo-controls label[for="geo-continent"]{
      border-color:var(--accent);
      background:var(--accent-soft);
      color:#0f3ea8;
      font-weight:700;
    }

    .table-wrap{
      width:100%;
      overflow:auto;
      border:1px solid var(--line);
      border-radius:14px;
      background:#fff;
    }

    table.tabela{
      width:100%;
      border-collapse:collapse;
      min-width:980px;
      background:#fff;
    }

    .tabela thead th{
      position:sticky;
      top:0;
      z-index:1;
      background:#f9fafb;
      border-bottom:1px solid var(--line);
      padding:10px 12px;
      text-align:left;
      font-size:.93rem;
      white-space:nowrap;
    }

    .tabela td{
      border-bottom:1px solid var(--line);
      padding:10px 12px;
      vertical-align:top;
      font-size:.93rem;
    }

    .tabela tbody tr:hover{background:#fafcff}

    .len-ok{
      display:inline-block;
      padding:2px 8px;
      border-radius:999px;
      background:var(--ok-bg);
      color:var(--ok);
      font-weight:700;
    }

    .len-below,.len-above{
      display:inline-block;
      padding:2px 8px;
      border-radius:999px;
      background:var(--warn-bg);
      color:var(--warn);
      font-weight:700;
    }

    .host, .geo{
      display:none;
      white-space:pre-wrap;
      word-break:break-word;
    }

    #host-species:checked ~ .table-wrap .host--species{display:inline}
    #host-order:checked   ~ .table-wrap .host--order{display:inline}
    #host-class:checked   ~ .table-wrap .host--class{display:inline}

    #geo-country:checked   ~ .table-wrap .geo--country{display:inline}
    #geo-region:checked    ~ .table-wrap .geo--region{display:inline}
    #geo-continent:checked ~ .table-wrap .geo--continent{display:inline}

    .empty-state{
      padding:18px;
      border:1px dashed var(--line);
      border-radius:12px;
      background:#fcfcfd;
      color:var(--muted);
    }

    @media print{
      body{background:#fff}
      .container{max-width:none;padding:0}
      .card{box-shadow:none;border:1px solid #ddd}
      .tab-controls,.toolbar{display:none !important}
      .tab-panels .tab-panel{display:block !important;page-break-inside:avoid}
      .table-wrap{overflow:visible}
    }
  </style>
</head>
<body>
  <header class="container">
    <h1><!--FILENAME--> Final Report</h1>
    <p class="muted">HTML interactive report for Consensus Multi-Sample mode</p>
  </header>

  <main class="container">
    <section class="tabs card">

      <input type="radio" id="tab1" name="tab" checked>
      <input type="radio" id="tab2" name="tab">
      <input type="radio" id="tab3" name="tab">
      <!--TAB4_INPUT-->

      <div class="tab-controls">
        <label for="tab1">Segment Background</label>
        <label for="tab2">Nextclade</label>
        <label for="tab3">FluMut</label>
        <!--TAB4_LABEL-->
      </div>

      <div class="tab-panels">

        <article id="panel1" class="tab-panel">
          <h2>Segment Background</h2>
          <p class="muted">
            Segment-level identification summary for all sequences in the multi-FASTA, combining
            cd-hit (<!--TOOL_VERSION_CD_HIT_EST_2D-->) and BLAST (<!--TOOL_VERSION_BLASTN-->)
            assignments. Host and geography views can be toggled below. To examine cluster
            composition, please refer to the cluster composition TSV file.
          </p>

          <input type="radio" id="host-species" name="host-view" checked>
          <input type="radio" id="host-order" name="host-view">
          <input type="radio" id="host-class" name="host-view">

          <input type="radio" id="geo-country" name="geo-view" checked>
          <input type="radio" id="geo-region" name="geo-view">
          <input type="radio" id="geo-continent" name="geo-view">

          <div class="toolbar">
            <div class="toolbar-group host-controls">
              <span class="label">Host:</span>
              <label for="host-species">Species</label>
              <label for="host-order">Order</label>
              <label for="host-class">Class</label>
            </div>

            <div class="toolbar-group geo-controls">
              <span class="label">Geography:</span>
              <label for="geo-country">Country</label>
              <label for="geo-region">Region</label>
              <label for="geo-continent">Continent</label>
            </div>
          </div>

          <!--MULTI_SEGMENT_BACKGROUND_TABLE-->
        </article>

        <article id="panel2" class="tab-panel">
          <h2>Nextclade</h2>
          <p class="muted">
            Clade calls reported by Nextclade (<!--TOOL_VERSION_NEXTCLADE-->) for eligible HA sequences.
          </p>
          <!--MULTI_NEXTCLADE_TABLE-->
        </article>

        <article id="panel3" class="tab-panel">
          <h2>FluMut</h2>
          <p class="muted">
            Obtained using FluMut (<!--TOOL_VERSION_FLUMUT-->) running FluMutDB (<!--TOOL_VERSION_FLUMUTDB-->).
            These mutations were further filtered for mutations of importance as cited by Alvarez et al. 2025 - 34 key
            mutations and 6 mutations of secondary importance - and Mohaptra et al. 2023 - for all antiviral resistance related mutations. 
            For a list of all mutations found, please refer to the final FluMut spreadsheet report.
            For a list of the mutations of interest parsed in this pipeline, please refer to the """ + _GITHUB_REPOSITORY_LINK + r""".
          </p>
          <!--MULTI_FLUMUT_TABLE-->
        </article>

        <!--TAB4_PANEL-->

      </div>
    </section>
  </main>
</body>
</html>
"""
def _read_fasta_lengths(fasta_fp: str) -> dict[str, int]:
    """
    Read a FASTA file and return {header_with_> : sequence_length}.

    Keeps headers with the leading '>' to stay compatible with the existing
    mappings convention used elsewhere in the pipeline.
    """
    lengths = {}

    if not fasta_fp or not os.path.exists(fasta_fp):
        return lengths

    current_header = None
    current_len = 0

    with open(fasta_fp, "r", encoding="utf-8") as fh:
        for raw_line in fh:
            line = raw_line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    lengths[current_header] = current_len
                current_header = line
                current_len = 0
            else:
                current_len += len(line)

    if current_header is not None:
        lengths[current_header] = current_len

    return lengths


def _build_original_seq_len_map(formatted_fasta_fp: str, mappings: dict[str, str]) -> dict[str, int]:
    """
    Build {normalized_original_sample_name: sequence_length} from the formatted
    FASTA and the existing mappings dict {internal_header: original_header}.
    """
    fasta_lengths = _read_fasta_lengths(formatted_fasta_fp)

    out = {}
    for internal_header, original_header in mappings.items():
        if internal_header in fasta_lengths:
            norm_original = str(original_header).strip()
            if norm_original.startswith('>'):
                norm_original = norm_original[1:].strip()
            out[norm_original] = fasta_lengths[internal_header]

    return out

def _segment_num_to_name(segment_cell: Any) -> str:
    """
    Convert an ID report segment value into a canonical segment name.

    Parameters
    ----------
    segment_cell : Any
        Segment cell value, including plain integers/strings or serialized
        dictionary-like values such as ``"{4: 1.0}"``.

    Returns
    -------
    str
        One of ``PB2``, ``PB1``, ``PA``, ``HA``, ``NP``, ``NA``, ``MP``,
        ``NS``, or ``"Not Available"`` when parsing fails.
    """
    seg_map = {value: key for key, value in segment_to_number_iav.items()}

    seg_num = None

    if pd.isna(segment_cell):
        return "Not Available"

    value = segment_cell
    if not isinstance(value, dict):
        s = str(value).strip()
        if not s:
            return "Not Available"
        try:
            value = ast.literal_eval(s) if s.startswith("{") else s
        except Exception:
            value = s

    if isinstance(value, dict) and value:
        key = next(iter(value.keys()))
        try:
            seg_num = int(key)
        except Exception:
            seg_num = None
    else:
        try:
            seg_num = int(value)
        except Exception:
            seg_num = None

    return seg_map.get(seg_num, "Not Available")


def create_segment_background_table_multi(
    df_segments: str,
    seg_lens: dict,
    mappings: dict[str, str],
    formatted_fasta_fp: str,
    geo_dict: dict[str, tuple[str, str]],
    taxa_dict: dict[str, list[str]],
) -> str:
    """
    Create Segment Background table for consensus multi-sample mode.

    Data sources
    ------------
    - df_segments: path to *_ID_Report.txt
    - formatted_fasta_fp: path to format_<file_tag>.fasta
    - mappings: dict {internal_header: original_header}

    Notes
    -----
    - Segment name is mined directly from the SEGMENT column of the ID report.
    - Sequence length is mined directly from the formatted FASTA.
    - Host / Geography / Genotype normalization logic is preserved.
    """
    len_thres = {
        'PB2': (seg_lens['PB2_min'], seg_lens['PB2_max']),
        'PB1': (seg_lens['PB1_min'], seg_lens['PB1_max']),
        'PA':  (seg_lens['PA_min'],  seg_lens['PA_max']),
        'HA':  (seg_lens['HA_min'],  seg_lens['HA_max']),
        'NP':  (seg_lens['NP_min'],  seg_lens['NP_max']),
        'NA':  (seg_lens['NA_min'],  seg_lens['NA_max']),
        'MP':  (seg_lens['MP_min'],  seg_lens['MP_max']),
        'NS':  (seg_lens['NS_min'],  seg_lens['NS_max']),
    }

    df, error = _safe_read_table(df_segments, required_cols=['SAMPLE_NAME', 'SEGMENT'])
    if error:
        return _empty_card(f'No segment background data available. {error}')

    for col, default in {
        'ASSIGNED_BY': 'NA',
        '%ID': 'NA',
        'CLUSTER': 'NA',
        'HOST': 'Unknown',
        'GENOTYPE': 'Unknown',
        'COUNTRY': 'Unknown',
        'REPRESENTATIVE': 'Not Available',
        'COLLECTION_DATE': 'Unknown',
    }.items():
        if col not in df.columns:
            df[col] = default

    seq_len_map = _build_original_seq_len_map(formatted_fasta_fp, mappings)

    # Derive canonical segment name directly from ID_Report
    df['SEGMENT_NAME'] = df['SEGMENT'].apply(_segment_num_to_name)
    df = df.sort_values(
        by='SEGMENT_NAME',
        key=lambda series: series.map(lambda value: _segment_sort_key(value)),
        kind='stable',
    ).reset_index(drop=True)

    # Normalize SAMPLE_NAME before length mapping
    df['SAMPLE_NAME_NORM'] = df['SAMPLE_NAME'].astype(str).str.strip()
    df['SAMPLE_NAME_NORM'] = df['SAMPLE_NAME_NORM'].str.replace(r'^>', '', regex=True)

    #print("DEBUG formatted_fasta_fp:", formatted_fasta_fp)
    #print("DEBUG fasta exists:", os.path.exists(formatted_fasta_fp))
    #print("DEBUG seq_len_map sample:", list(seq_len_map.items())[:5])
    #print("DEBUG SAMPLE_NAME sample:", df['SAMPLE_NAME'].astype(str).tolist()[:5])
    #print("DEBUG SAMPLE_NAME_NORM sample:", df['SAMPLE_NAME_NORM'].tolist()[:5])

    # Mine sequence lengths directly from formatted FASTA
    df['SEQ_LEN'] = df['SAMPLE_NAME_NORM'].map(seq_len_map).fillna(0).astype(int)

    #print("DEBUG SEQ_LEN sample:", df['SEQ_LEN'].tolist()[:5])

    # Prepare cleaned metadata views for clustered / C-BLAST rows
    mask = df['ASSIGNED_BY'] != 'L-BLAST'

    df['HOST_CLEAN'] = None
    df['GENOTYPE_CLEAN'] = None
    df['COUNTRY_CLEAN'] = None

    df.loc[mask, 'HOST_CLEAN'] = df.loc[mask, 'HOST'].apply(
        lambda x: clean_host(x, taxa_dict, level='Species', drop_unknown=True)
    )
    df.loc[mask, 'GENOTYPE_CLEAN'] = df.loc[mask, 'GENOTYPE'].apply(
        lambda x: clean_genotype(x)
    )
    df.loc[mask, 'COUNTRY_CLEAN'] = df.loc[mask, 'COUNTRY'].apply(
        lambda x: clean_country(x, geo_dict, mode='country', drop_unknown=True)
    )

    table = """
    <div class="table-wrap">
    <table class="tabela" id="sb-multi-table">
    <thead>
      <tr>
        <th>Sample Name</th>
        <th>Segment</th>
        <th>SeqLen</th>
        <th>Note</th>
        <th>Cluster</th>
        <th>Reference</th>
        <th>%ID</th>
        <th>Assigned_By</th>
        <th>Genotype</th>
        <th>Host</th>
        <th>Country</th>
        <th>Collection_date</th>
      </tr>
    </thead>
    <tbody>
    """

    for row in range(len(df)):
        sample_name = str(df.loc[row, 'SAMPLE_NAME'])
        segment = df.loc[row, 'SEGMENT_NAME']
        length = int(df.loc[row, 'SEQ_LEN']) if not pd.isna(df.loc[row, 'SEQ_LEN']) else 0
        reference = _representative_link(df.loc[row, 'REPRESENTATIVE']) if 'REPRESENTATIVE' in df.columns else "Not Available"
        perc_id = df.loc[row, '%ID'] if '%ID' in df.columns else "NA"
        cluster = df.loc[row, 'CLUSTER'] if 'CLUSTER' in df.columns else "NA"
        collection_date = _normalize_collection_date(df.loc[row, 'COLLECTION_DATE']) if 'COLLECTION_DATE' in df.columns else "Unknown"
        assigned_by = df.loc[row, 'ASSIGNED_BY'] if 'ASSIGNED_BY' in df.columns else "NA"

        note_parts = []

        if segment in len_thres:
            min_len, max_len = len_thres[segment]

            if length < min_len:
                length_class = "len-below"
                note_parts.append("Length outside of expected range")
            elif length > max_len:
                length_class = "len-above"
                note_parts.append("Length outside of expected range")
            else:
                length_class = "len-ok"

            length_html = f'<span class="{length_class}">{length}</span>'
        else:
            length_html = '<span class="len-below">0</span>'
            note_parts.append("Segment not found")

        if assigned_by != 'L-BLAST':
            host_species = df.loc[row, 'HOST_CLEAN']
            genotype = df.loc[row, 'GENOTYPE_CLEAN']
            geo_countries = df.loc[row, 'COUNTRY_CLEAN']

            host_species = host_species if isinstance(host_species, dict) else {}
            genotype = genotype if isinstance(genotype, dict) else {}
            geo_countries = geo_countries if isinstance(geo_countries, dict) else {}

            host_order = rollup_counts_taxa(host_species, taxa_dict, level='Order', drop_unknown=True)
            host_class = rollup_counts_taxa(host_species, taxa_dict, level='Class', drop_unknown=True)

            geo_regions = rollup_counts_geo(geo_countries, geo_dict, mode='region', drop_unknown=True)
            geo_continents = rollup_counts_geo(geo_countries, geo_dict, mode='continent', drop_unknown=True)
        else:
            host_species = df.loc[row, 'HOST']
            genotype = df.loc[row, 'GENOTYPE']
            geo_countries = df.loc[row, 'COUNTRY']

            host_order = rollup_taxa(host_species, taxa_dict, level='Order')
            host_class = rollup_taxa(host_species, taxa_dict, level='Class')

            geo_regions = rollup_geo(geo_countries, geo_dict, mode='region')
            geo_continents = rollup_geo(geo_countries, geo_dict, mode='continent')

        note = "; ".join(note_parts)

        table += f"""
        <tr>
          <td>{sample_name}</td>
          <td>{segment}</td>
          <td>{length_html}</td>
          <td>{note}</td>
          <td>{cluster}</td>
          <td>{reference}</td>
          <td>{perc_id}</td>
          <td>{assigned_by}</td>
          <td>{genotype}</td>
          <td>
            <span class="host host--species">{host_species}</span>
            <span class="host host--order">{host_order}</span>
            <span class="host host--class">{host_class}</span>
          </td>
          <td>
            <span class="geo geo--country">{geo_countries}</span>
            <span class="geo geo--region">{geo_regions}</span>
            <span class="geo geo--continent">{geo_continents}</span>
          </td>
          <td>{collection_date}</td>
        </tr>
        """

    table += """
    </tbody>
    </table>
    </div>
    """

    return table

def create_nextclade_table_multi(reports_p: str, file_tag: str) -> str:
    """
    Create a compact Nextclade table for consensus multi-sample mode.

    Reads any existing Nextclade TSVs for the current file_tag and reports only:
    - Dataset (derived from filename: H1/H3/H5)
    - seqName
    - clade
    """
    nextclade_files = [
        ("H1", os.path.join(reports_p, f"{file_tag}_H1_nextclade.tsv")),
        ("H3", os.path.join(reports_p, f"{file_tag}_H3_nextclade.tsv")),
        ("H5", os.path.join(reports_p, f"{file_tag}_H5_nextclade.tsv")),
    ]

    frames = []

    for dataset, fp in nextclade_files:
        if not os.path.exists(fp):
            continue

        try:
            df = pd.read_table(fp, index_col=False)
        except Exception:
            continue

        if df is None or df.empty:
            continue

        if "seqName" not in df.columns:
            continue

        if "clade" not in df.columns:
            df["clade"] = "N/A"

        keep = df[["seqName", "clade"]].copy()
        keep.insert(0, "Dataset", dataset)
        keep = keep.fillna("N/A")

        frames.append(keep)

    if not frames:
        return _empty_card('No Nextclade data available.')

    merged = pd.concat(frames, ignore_index=True)

    table = """
    <div class="table-wrap">
    <table class="tabela" id="nextclade-multi-table">
    <thead>
      <tr>
        <th>Dataset</th>
        <th>seqName</th>
        <th>clade</th>
      </tr>
    </thead>
    <tbody>
    """

    for row in merged.itertuples(index=False):
        dataset = str(row.Dataset)
        seq_name = str(row.seqName)
        clade = str(row.clade)

        table += f"""
        <tr>
          <td>{dataset}</td>
          <td>{seq_name}</td>
          <td>{clade}</td>
        </tr>
        """

    table += """
    </tbody>
    </table>
    </div>
    """

    return table
def create_flumut_table_multi(reports_p: str, file_tag: str, muts_interest: dict[str, Iterable[str]], muts_loci_meaning: dict) -> str:
    """
    Create a FluMut table for consensus multi-sample mode using only markers.tsv.

    Behaviour
    ---------
    - Reads <file_tag>_markers.tsv
    - Uses only mutations present in muts_interest
    - Reuses the same locus aliasing / mutation normalization logic used by mut_miner():
        * A123T -> 123T
        * PB1-F2 -> PB1
        * PA-X -> PA
        * HA1-5 / HA2-5 -> HA
        * NA-1 / NA-2 -> NA
        * M1 / M2 -> MP
        * NS-1 / NS-2 -> NS

    Output columns
    --------------
    - Sample
    - Segment/Locus
    - Mutation
    - Biological effects
    """
    markers_fp = os.path.join(reports_p, f"{file_tag}_markers.tsv")

    if not os.path.exists(markers_fp):
        return _empty_card('No FluMut marker data available.')

    df, error = _safe_read_table(markers_fp, required_cols=['Sample', 'Mutations in your sample'])
    if error:
        return _empty_card(f'FluMut marker data unavailable. {error}')

    # Same alias logic as mut_miner()
    loci_seg = {
        'PB1-F2': 'PB1',
        'PA-X': 'PA',
        'HA1-5': 'HA',
        'HA2-5': 'HA',
        'NA-1': 'NA',
        'NA-2': 'NA',
        'M1': 'MP',
        'M2': 'MP',
        'NS-1': 'NS',
        'NS-2': 'NS'
    }

    # Regex used in mut_miner(): A123T -> 123T
    mut_patt = r'([A-Z])(\d+)([A-Z])'

    # Flatten muts_interest to a quick lookup:
    # {segment: set_of_interest_mutations}
    interest_lookup = {}
    for seg, muts in muts_interest.items():
        try:
            interest_lookup[seg] = set(muts)
        except Exception:
            interest_lookup[seg] = set()

    hits = []

    for row in df.itertuples(index=False):
        sample = getattr(row, "Sample", None)
        raw_mut = getattr(row, "_asdict")().get("Mutations in your sample", None) if hasattr(row, "_asdict") else None

        # Fallback for tuple access if needed
        if raw_mut is None:
            try:
                raw_mut = row[df.columns.get_loc("Mutations in your sample")]
            except Exception:
                raw_mut = None

        if pd.isna(sample) or pd.isna(raw_mut):
            continue

        sample = str(sample).strip()
        raw_mut = str(raw_mut).strip()

        if not sample or not raw_mut:
            continue

        # Expected format from FluMut markers.tsv: "Locus:Mutation"
        # Example: "HA:A156T"
        if ":" not in raw_mut:
            continue

        locus_raw, mut_raw = raw_mut.split(":", 1)
        locus_raw = str(locus_raw).strip()
        mut_raw = str(mut_raw).strip()

        if not locus_raw or not mut_raw:
            continue

        # Canonical segment for matching against muts_interest
        segment = loci_seg.get(locus_raw, locus_raw)
        display_locus = locus_raw

        # Normalize mutation like mut_miner(): A123T -> 123T
        norm_mut = mut_raw
        ex = re.match(mut_patt, norm_mut)
        if ex:
            norm_mut = f"{ex.group(2)}{ex.group(3)}"

        # Keep only mutations of interest for the canonical segment
        if segment not in interest_lookup:
            continue

        if norm_mut not in interest_lookup[segment]:
            continue

        # Biological effect lookup follows the same convention already used elsewhere
        effect = muts_loci_meaning.get(norm_mut, ["", ""])
        if isinstance(effect, (list, tuple)) and len(effect) > 1:
            biological_effect = effect[1]
        else:
            biological_effect = ""

        hits.append({
            "Sample": sample,
            "Segment/Locus": display_locus,
            "Mutation": norm_mut,
            "Biological effects": biological_effect
        })

    if not hits:
        return """
        <div class="card">
          <p class="muted">No mutations of interest found in FluMut markers.</p>
        </div>
        """

    hits_df = pd.DataFrame(hits).drop_duplicates()
    hits_df['__sample_sort'] = hits_df['Sample'].astype(str)
    hits_df['__locus_sort'] = hits_df['Segment/Locus'].map(_mutation_locus_sort_key)
    hits_df['__mutation_sort'] = hits_df['Mutation'].astype(str)
    hits_df = hits_df.sort_values(
        by=["__sample_sort", "__locus_sort", "__mutation_sort"],
        kind="stable"
    ).drop(columns=['__sample_sort', '__locus_sort', '__mutation_sort'])

    table = """
    <div class="table-wrap">
    <table class="tabela" id="flumut-multi-table">
    <thead>
      <tr>
        <th>Sample</th>
        <th>Segment/Locus</th>
        <th>Mutation</th>
        <th>Biological effects</th>
      </tr>
    </thead>
    <tbody>
    """

    for row in hits_df.itertuples(index=False):
        sample = str(row[0])
        segment = str(row[1])
        mutation = str(row[2])
        effect = str(row[3]) if row[3] is not None else ""

        table += f"""
        <tr>
          <td>{sample}</td>
          <td>{segment}</td>
          <td>{mutation}</td>
          <td>{effect}</td>
        </tr>
        """

    table += """
    </tbody>
    </table>
    </div>
    """

    return table
def to_html_report_multi(
    filename: str,
    seg_back_table: str,
    nextclade_table: str,
    flumut_table: str,
    output_file_path: str,
    genin2: bool = False,
    genin2table: str = '',
    tool_versions: Optional[dict[str, Any]] = None,
    tool_version_placeholder_tag: Optional[str] = None,
    html_skeleton_multi: str = html_skeleton_multi
) -> None:
    """
    Render the consensus multi-sample final HTML report.

    Parameters
    ----------
    filename : str
        Base filename/tag to display in the report header.
    seg_back_table : str
        Pre-rendered HTML for the Segment Background tab.
    nextclade_table : str
        Pre-rendered HTML for the Nextclade tab.
    flumut_table : str
        Pre-rendered HTML for the FluMut tab.
    output_file_path : str
        Destination HTML file path.
    genin2 : bool, default False
        Whether to enable the optional GenIn2 tab.
    genin2table : str, default ''
        Pre-rendered HTML for the GenIn2 tab.
    tool_versions : dict, optional
        Optional tool/version metadata, typically ``flags['Tool Versions']``.
        When combined with ``tool_version_placeholder_tag``, the helper
        ``apply_tool_version_placeholder()`` replaces that caller-provided tag.
    tool_version_placeholder_tag : str, optional
        Exact placeholder string to replace with rendered tool versions. No
        built-in placeholder is injected into the skeleton automatically.
    html_skeleton_multi : str
        HTML skeleton template for multi-sample mode.
    """
    html = html_skeleton_multi

    html = html.replace('<!--FILENAME-->', str(filename))
    html = html.replace('<!--MULTI_SEGMENT_BACKGROUND_TABLE-->', seg_back_table)
    html = html.replace('<!--MULTI_NEXTCLADE_TABLE-->', nextclade_table)
    html = html.replace('<!--MULTI_FLUMUT_TABLE-->', flumut_table)

    if genin2:
        html = html.replace(
            '<!--TAB4_INPUT-->',
            '<input type="radio" id="tab4" name="tab">'
        )
        html = html.replace(
            '<!--TAB4_LABEL-->',
            '<label for="tab4">Genin2</label>'
        )
        html = html.replace(
            '<!--TAB4_ACTIVE_CSS-->',
            ',\n    #tab4:checked ~ .tab-controls label[for="tab4"]'
        )
        html = html.replace(
            '<!--TAB4_DISPLAY_CSS-->',
            '#tab4:checked ~ .tab-panels #panel4{display:block}'
        )
        html = html.replace(
            '<!--TAB4_PANEL-->',
            f'''
        <article id="panel4" class="tab-panel">
          <h2>Genin2</h2>
          <p class="muted">
            GenIn2 performs genotype constellation assignment. This tab reports
            constellation assignments for eligible samples in this analysis using
            GenIn2 (<!--TOOL_VERSION_GENIN2-->).
          </p>
          {genin2table}
        </article>
        '''
        )
    else:
        html = html.replace('<!--TAB4_INPUT-->', '')
        html = html.replace('<!--TAB4_LABEL-->', '')
        html = html.replace('<!--TAB4_ACTIVE_CSS-->', '')
        html = html.replace('<!--TAB4_DISPLAY_CSS-->', '')
        html = html.replace('<!--TAB4_PANEL-->', '')

    html = apply_tool_version_placeholder(
        html,
        tool_versions=tool_versions,
        placeholder_tag=tool_version_placeholder_tag,
    )

    with open(output_file_path, 'w', encoding='utf-8') as f:
        f.write(html)

def generate_final_report_single(
    df_segments: str,
    flags: dict,
    seg_lens: dict,
    mappings: dict[str, str],
    muts_loci_meaning: dict,
    html_skeleton: str,
    output_file_path: str,
    filename: str,
    geo_dict: dict[str, tuple[str, str]],
    taxa_dict: dict[str, list[str]],
) -> None:
    """
    Generate the existing single-sample final HTML report.
    """
    sample_card = create_sample_card(flags)

    try:
        sample_table = create_sample_table(df_segments, flags, seg_lens, mappings, muts_loci_meaning)
    except Exception as exc:
        sample_table = _empty_card(f'Sample segment summary unavailable. {exc}')

    try:
        seg_back_table = create_segment_background_table(
                            df_segments,
                            flags,
                            seg_lens,
                            mappings,
                            geo_dict=geo_dict,
                            taxa_dict=taxa_dict,
                        )
    except Exception as exc:
        seg_back_table = _empty_card(f'Segment background unavailable. {exc}')

    try:
        mut_table = create_segment_mutations_table(flags, muts_loci_meaning)
    except Exception as exc:
        mut_table = _empty_card(f'Segment mutation summary unavailable. {exc}')

    if flags['Master'].get('genin', False):
        genin2table = create_genin_const_table(flags)
        to_html_report(
            filename=filename,
            sample_card=sample_card,
            sample_table=sample_table,
            seg_back_table=seg_back_table,
            mut_table=mut_table,
            output_file_path=output_file_path,
            genin2=True,
            genin2table=genin2table,
            tool_versions=flags.get('Tool Versions', {}),
        )
    else:
        to_html_report(
            filename=filename,
            sample_card=sample_card,
            sample_table=sample_table,
            seg_back_table=seg_back_table,
            mut_table=mut_table,
            output_file_path=output_file_path,
            genin2=False,
            tool_versions=flags.get('Tool Versions', {}),
        )


def generate_final_report_multi(
    df_segments: str,
    flags: dict,
    seg_lens: dict,
    mappings: dict[str, str],
    muts_loci_meaning: dict,
    html_skeleton_multi: str,
    output_file_path: str,
    filename: str,
    runs_p: str,
    geo_dict: dict[str, tuple[str, str]],
    taxa_dict: dict[str, list[str]],
) -> None:
    """
    Generate the consensus multi-sample final HTML report.

    Expected inputs
    ---------------
    df_segments : str
        Path to <file_tag>_ID_Report.txt
    output_file_path : str
        Path to the destination final HTML report
    filename : str
        File tag / display name for the report header
    """
    reports_p = os.path.dirname(os.path.abspath(df_segments))
    file_tag = os.path.basename(df_segments).replace('_ID_Report.txt', '')
    formatted_fasta_fp = os.path.join(runs_p, f'format_{file_tag}.fasta')

    try:
        seg_back_table = create_segment_background_table_multi(
                        df_segments=df_segments,
                        seg_lens=seg_lens,
                        mappings=mappings,
                        formatted_fasta_fp=formatted_fasta_fp,
                        geo_dict=geo_dict,
                        taxa_dict=taxa_dict,
                        )
    except Exception as exc:
        seg_back_table = _empty_card(f'Segment background unavailable. {exc}')

    try:
        nextclade_table = create_nextclade_table_multi(
            reports_p=reports_p,
            file_tag=file_tag
        )
    except Exception as exc:
        nextclade_table = _empty_card(f'Nextclade summary unavailable. {exc}')

    try:
        flumut_table = create_flumut_table_multi(
            reports_p=reports_p,
            file_tag=file_tag,
            muts_interest=muts_interest,
            muts_loci_meaning=muts_loci_meaning
        )
    except Exception as exc:
        flumut_table = _empty_card(f'FluMut summary unavailable. {exc}')

    if flags['Master'].get('genin', False):
        genin2table = create_genin_const_table(flags)
        to_html_report_multi(
            filename=filename,
            seg_back_table=seg_back_table,
            nextclade_table=nextclade_table,
            flumut_table=flumut_table,
            output_file_path=output_file_path,
            genin2=True,
            genin2table=genin2table,
            tool_versions=flags.get('Tool Versions', {}),
            html_skeleton_multi=html_skeleton_multi
        )
    else:
        to_html_report_multi(
            filename=filename,
            seg_back_table=seg_back_table,
            nextclade_table=nextclade_table,
            flumut_table=flumut_table,
            output_file_path=output_file_path,
            genin2=False,
            tool_versions=flags.get('Tool Versions', {}),
            html_skeleton_multi=html_skeleton_multi
        )

def generate_final_report(
    df_segments: str,
    flags: dict,
    seg_lens: dict,
    mappings: dict[str, str],
    muts_loci_meaning: dict,
    html_skeleton: str,
    output_file_path: str,
    filename: str,
    metadata_p: str,
    geo: str,
    taxa: str,
    runs_p: Optional[str] = None
) -> None:
    """
    Dispatcher for final report generation.

    Single-sample mode:
        uses the original single-sample HTML report

    Multi-sample mode:
        uses the consensus multi-sample HTML report
    """
    geo_dict, taxa_dict= load_reference_dicts(
            data_dir=metadata_p,
            geo_json_relpath=geo,
            tax_json_relpath=taxa
            )
    is_single = flags.get('Master', {}).get('single', True)

    if is_single:
        generate_final_report_single(
            df_segments,
            flags,
            seg_lens,
            mappings,
            muts_loci_meaning,
            html_skeleton,
            output_file_path,
            filename,
            geo_dict=geo_dict,
            taxa_dict=taxa_dict,
        )
    else:
        generate_final_report_multi(
            df_segments,
            flags,
            seg_lens,
            mappings,
            muts_loci_meaning,
            html_skeleton,
            output_file_path,
            filename,
            runs_p=runs_p,
            geo_dict=geo_dict,
            taxa_dict=taxa_dict,
        )
