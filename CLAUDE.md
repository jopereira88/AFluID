# AFluID — orientation

Influenza A segment / subtype / genotype identification pipeline. Assignment is
CD-HIT-est-2d clustering at 99% against a curated cluster DB, then BLAST
fallbacks. Downstream tools are FluMut, GenIn2, Nextclade and GenBank reference
retrieval. Output is TSV + HTML reports. Pure Python scripts, no package
structure, run from the project root inside the `fluid` conda env (`env.yaml`).

Line numbers below are approximate (snapshot 2026-11-04). Use them as a
starting point and `grep -n "def <name>"` to confirm.

## File map

| File | Lines | Role |
|---|---|---|
| `main.py` | ~3550 | CLI entry point and the whole run pipeline (preprocess → cluster → BLAST → report → tools) |
| `final_report_utils.py` | ~3400 | HTML/TSV report generation, roll-ups, batch index/zip, HTML skeletons |
| `metadata_utils.py` | ~1170 | Table classes: `MetadataTable` (ABC), `SequenceMetadata` (CSV), `ClusterMetadata`, `BlastReportTable`, `ClusterReportTable`, `BlastClustReportTable` |
| `flu_utils.py` | ~1000 | FASTA I/O (`seq_get`, `dict_to_fasta`, `concat_fasta`), `parse_clstr`, `sp_blastn`, JSON helpers, H/N genotype mining |
| `update.py` | ~840 | DB update maintenance command (metadata cleaning → QC → AFluID run → cross-check → concat → re-install) |
| `gb_utils.py` | ~530 | GenBank fetch/filter/convert (Biopython Entrez) |
| `metadata_tree_utils.py` | ~490 | Loads `geo_tree.json` / `host_taxonomy_tree.json` into `geo_dict` / `taxa_dict` (`load_reference_dicts`) |
| `structures.py` | ~480 | Static data: `flagdict` template, `muts_interest`, `muts_loci_meaning`, segment maps, length thresholds, `geo_syns`, `taxa_syns` |
| `install.py` | ~240 | Builds dirs, BLAST DB (`makeblastdb`), cluster DB (`cd-hit-est`), `cluster_charact()` cluster metadata |
| `dbSetup.py` | ~40 | Filters a FASTA down to accessions present in the metadata CSV |
| `test_*.py` | | `unittest` suites: HTML reports, TSV exports, metadata trees, collection-date compatibility |
| `config_example.ini` | | Template for `config.ini` (gitignored). Sections: Paths, Sequence_Size, CD-HIT, blast, getref, Filenames, Functions, Versions |

Gitignored data dirs: `blast_db/ cluster_db/ metadata/ references/ reports/ runs/ samples/ logs/`, plus `update/`.

## main.py layout

- **Tool plumbing (~58–190):** `_set_tool_status`, `_write_flags_json`, the tool-version capture helpers (`_capture_tool_version`, `_record_tool_versions`), `_run_command`, `_require_files`.
- **Identification stages:**
  - `fasta_preprocess` (~208) renames headers to `Seq-N` via `seq_enum` and returns the `mappings` dict (internal → original header). It writes `runs/.../format_<tag>.fasta`.
  - `cd_hit_est_2d` (~291) → `cluster_assign` (~346) → `cluster_compile` / `cluster_miner` (~423/491) produce `<tag>_clust_report.txt`.
  - `bclust` (~611) runs BLAST vs. cluster representatives (C_BLAST, ≥90% id). `best_blast` (~581) picks the best hit.
  - `reblast` (~699) runs BLAST vs. the bulk local DB (L_BLAST).
  - `report_compiler` (~772) writes `<tag>_ID_Report.txt`.
  - `compute_best_segments_iterative` (~886) picks the best contig per segment (contig + single mode only).
- **Routing:** `redirector` (~1036) reads the ID report and fills `flags['Final Report']` targets for FluMut, Nextclade (H1/H3/H5), GenIn and GetRef.
- **Downstream tools:**
  - FluMut: `conform_to_flumut`, `run_flumut`, `remap_flumut_report`, `mut_miner` (~1234–1535)
  - GetRef: `get_reference` (~1536)
  - Nextclade: `fasta_to_nextclade`, `run_nextclade`, `remap_nextclade`, `mine_clade` (~1589–1734)
  - GenIn2: `to_genin2`, `run_genin2`, `remap_genin2` (~1735–1895)
- **CLI and paths:** `parser` (~1896) holds all flags and validation rules. `resolve_*` helpers handle config, outdir, single-sample and batch output paths. `make_output_stem`, `make_file_tag` and `make_batch_sample_tag` build names.
- **Batch:**
  - `run_batch_pipeline` (~2081) writes `batch_summary_<name>.tsv`, the batch cluster composition, the index and the zip.
  - `extract_batch_summary_fields` and `_normalize_genin_summary` fill the summary columns.
  - `find_fasta_files` lists sample files.
  - `demultiplex_batch_fasta` (~2658) handles `-bf --bf-regex`.
- **`run_pipeline_for_file` (~2721):** the per-file orchestrator. It deep-copies `flagdict`, preprocesses, then runs either the `--skip_cdhit` path (bclust → reblast) or the normal path (cd-hit → cluster_assign → cluster_miner → bclust → reblast). After that: ID report → redirector → roll-up / cluster-composition exports → FluMut → Nextclade → GetRef → GenIn2 → flags.json → final HTML report (via `generate_final_report`).
- **`main` (~3245):** reads the config and resolves `[Paths]` against the project root. It cleans `runs/` unless `-p` is set (or `remove_previous=False`) and creates the `runs/session_*` workspace (`create_run_workspace`). It then dispatches to batch (`-b`), batch-fasta (`-bf`) or single file.

## final_report_utils.py layout

- **Normalization/lookup (~24–245):** `_resolve_with_synonyms` (uses `structures.geo_syns`/`taxa_syns`), safe table readers.
- **Roll-ups (~245–915):** `export_id_report_rollup`, `export_cluster_composition`, `rollup_geo`, `rollup_taxa`, `clean_host`, `clean_country`, `clean_genotype`. Exported distributions round to 2 dp and sum to 1.0. L-BLAST rows stay scalar.
- **Batch artifacts (~943–1437):** `create_batch_index`, `create_batch_zip`, `maybe_create_batch_artifacts`.
- **Single-sample HTML:** `create_sample_card`, `create_sample_table`, `create_segment_background_table`, `create_segment_mutations_table`, `create_genin_const_table`, then `html_skeleton` (~1949) and `to_html_report` (~2132).
- **Multi-sample HTML (`-ss off`):** `html_skeleton_multi` (~2229), `create_*_multi`, `to_html_report_multi` (~3100).
- **Entry:** `generate_final_report` (~3344) dispatches to `_single` or `_multi`. Each must use its matching skeleton (a past bug passed the wrong one).
- **HTML placeholders:**
  - Content sections use `<!--SAMPLE_CARD-->`, `<!--SEGMENT_BACKGROUND_TABLE-->` and so on.
  - Tool versions use `<!--TOOL_VERSION_BLASTN|CD_HIT_EST_2D|FLUMUT|FLUMUTDB|GENIN2|NEXTCLADE-->`, filled by `apply_tool_version_placeholder`.

## Key concepts / invariants

- **`flags` dict:** a per-file deep copy of `structures.flagdict`. This is the central mutable state. `flags['Master']` holds the tool on/off switches (adjusted by mode, `-ff` and `-off`). `flags['Sample']` holds per-segment results, mutations and genotypes. The dict is dumped as `<tag>_flags.json`.
- **Tool Status / error flow:**
  - `flags['Tool Status'][name] = {status: completed|failed|skipped, detail}`.
  - Fatal stages run through `_run_stage` / the local `stage()` closure: `preprocess`, `cd_hit`, `cluster_parse`, `c_blast`, `l_blast`, `id_report`, `best_segments`, `redirector`, `rollup_export`, plus `final_report`. These mark the stage failed, save flags.json, and re-raise with `exc.afluid_stage` / `exc.afluid_flags`.
  - Non-fatal tools catch their own errors: `flumut`, `flumut_db_update`, `nextclade`, `getref`, `genin`.
  - `_run_command` appends a stderr tail to its error message.
  - Batch summary: `*_status` columns plus `tool_errors` (from `_format_tool_errors`). Batch-wide failures go to `batch_issues_<name>.tsv`. New tools need a status entry, and `_BATCH_TOOL_STATUS_FIELDS` should be updated when the tool should get its own column.
- **Routing vs. forcing:** `flags['Master'][tool]` means a tool is *enabled*. `forced_*` (from `-ff`) means it is *forced*, and that is what gets passed to `redirector`. Consensus + `-ss on` sends sequences to FluMut only when `H_gen == 'H5'`. Consensus + `-ss off` deliberately forces FluMut for every sequence.
- **Header handling:** samples are renamed internally to `Seq-N`, and every report gets remapped back to the original headers through `mappings`. Any new output that shows sample names needs that remap step. Remaps must tolerate segment suffixes like `Seq-1_HA` and a leading `>`.
- **Duplicate FASTA headers:** `flu_utils.seq_get` gives every duplicated header a suffix (`_A`, `_B`, … `_AA`), *including the first occurrence*. Unique headers stay unchanged.
- **Segments:** canonical order is `structures.iav_segments` (PB2, PB1, PA, HA, NP, NA, MP, NS), and reports sort by it. Number ↔ name mapping comes from `iupac_to_int` / `int_to_iupac`.
- **Metadata:** columns are looked up by header name, not by position. The geo/host hierarchy comes from the JSON trees, not from `structures.py`. Missing host ranks fall back to the nearest higher rank.
- **Mutations of interest:** `structures.muts_interest` and `muts_loci_meaning`. Labels are shown locus-prefixed (`PB1-F2:66S`). The README table must be kept in sync with these.
- **Paths:**
  - Relative `[Paths]` in the config resolve against the project root, not the cwd.
  - `--outdir` resolves against the cwd.
  - Update reports go to `[Paths] update_reports`, which falls back to `update/reports`.
- **Outputs:**
  - Single sample: `<reports>/<stem>/` or `<outdir>/<stem>/{reports,references}/`.
  - Batch: `<reports>/<batch_name>/<sample>/`, with the summary, index and zip at the batch root.
  - Existing targets raise an error unless `--replace` is given. `-o` names are used verbatim, with no timestamp or hex suffix.
- **Run workspaces:** each run gets its own `runs/session_*` directory. It is removed on success and kept on failure.

## Report files (per output tag)

- `<tag>_ID_Report.txt`. Columns: SAMPLE_NAME, SEGMENT, CLUSTER, REPRESENTATIVE, %ID, ASSIGNED_BY, GENOTYPE, HOST, COUNTRY, COLLECTION_DATE.
- `<tag>_ID_Report_rollup.tsv`, which adds HOST_FAMILY, HOST_CLASS, COUNTRY_REGION and COUNTRY_CONTINENT.
- `<tag>_clust_report.txt` and `<tag>_cluster_composition.tsv`.
- `<tag>_flags.json`, plus `<tag>_tool_versions.tsv` with `-tv`.
- FluMut: `<tag>_markers.tsv`, `_mutations.tsv`, `_literature.tsv`. Nextclade: `<tag>_H{1,3,5}_nextclade.tsv`. GenIn2 output.
- `<tag>_final_report.html`.

## Commands

```bash
conda activate fluid
python3 main.py -f SAMPLE.fasta -m consensus            # single sample (file in samples/)
python3 main.py -b -bd SUBDIR -m consensus -bn NAME     # batch over samples/SUBDIR
python3 main.py -f MULTI.fasta -m consensus -bf --bf-regex '(?P<sample>EPI_ISL_\d+)'
python3 install.py config.ini                           # build DBs
python3 update.py -c config.ini -ff UPDATE.fasta -mf UPDATE.csv
python3 -m unittest discover -p 'test_*.py'             # tests; need metadata/ + JSON trees present
```

## Local-only docs (gitignored, root folder)

- `CHANGELOG_2026-05-28.md`: roll-up TSVs, cluster composition, tool-version placeholders, column reordering, JSON trees.
- `CHANGELOG_2026-07-22.md`: `-rm` replaced by `-p/--parallel`, multi-sample skeleton fix, FluMut remap fix.
- `session-ses_*.md`: raw transcripts of earlier agent sessions. They are very large (6k–15k lines each), so grep them rather than reading them whole.
  - `1ba0`: `--outdir`, `--batch_name`, `--replace`
  - `1b5f`: batch summary columns, `-o`, per-run workspaces
  - `1afa`: JSON tree migration fixes
  - `1920`: roll-ups, column order, GenIn panel
  - `0e76`: duplicate-header suffixing

## Conventions

- Keep the README in sync when changing CLI flags, output layout or mutations of interest.
- The user writes dated `CHANGELOG_YYYY-MM-DD.md` summaries after work sessions.
- Validate changes with the `unittest` suites and with a real single-sample run, plus a batch run when batch logic changes.
