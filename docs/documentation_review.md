# Documentation review of `src/`

Companion to the docstring reformatting pass. Section 1 records what was changed in the
source. Section 2 is a **flag list only**: nothing in it has been written into the code,
the suggested text is there for you to accept, edit or reject.

Line numbers refer to the files as they stand after the pass.

> **Note**
> This review was written before the package refactor. The modules have since moved
> from `src/` to `src/coolpackage/`, so a path written here as `analysis.py` is now
> `src/coolpackage/analysis.py`. The line numbers are unchanged.

---

## 1. What was changed

### 1.1 Docstring format (rule 1)

Every docstring in `src/` and `src/features/` used a numpydoc-style layout
(`Parameters` / `Returns` / `Method` / `Example` headings with `---` underlines, and
either `name : type` or `name -> type` entries). `sphinx.ext.napoleon` is **not** loaded
in `conf.py`, so none of that was being interpreted: the sections were rendered as
undifferentiated prose and the parameter names lost.

All docstrings now use plain reStructuredText:

| Was | Is now |
| --- | --- |
| `Parameters` / `Args` section | `:param <name>:` plus `:type <name>:` fields |
| `Returns` section | `:returns:` plus `:rtype:` fields |
| `Method` section | `.. rubric:: Method` followed by a paragraph or a bullet list |
| `Example` section | `.. rubric:: Example` followed by a `::` literal block |
| `>> call()` prompts | `>>> call()` |
| ASCII "Outline" tables inside `Returns` | literal blocks, so the column alignment survives |

Types were normalised (`string` → `str`, `integer` → `int`, `dataframe` /
`Pandas DataFrame` → `pandas.DataFrame`), and text was rewrapped to 100 columns.

`docs` now builds with **zero Sphinx warnings** (`python -m sphinx -b html source build`).

### 1.2 `.. todo::` annotations (rule 3)

Every `TODO` comment was rewritten as a Sphinx `.. todo::` directive in the docstring of
the scope it belonged to, and the original comment removed. `todo_include_todos = True`
is already set in `conf.py`, so they are rendered.

| File | Scope the todo now lives in |
| --- | --- |
| `analysis.py` | new module docstring (three todos, previously a `#### TODO Section ####` comment block above the imports) |
| `aggregation.py` | `Aggregation._reduce_aev_dimensions` |
| `measure.py` | `Measure._match_resid_codes` (two, previously inline in the docstring prose), `Measure.restart_measure_pdb_only` (two) |
| `protein.py` | `PDB.load_state`, `PDB.clean_and_split_pdb` |
| `patcher.py` | `_fasta_to_pir`, `_full_align`, `_patch_model`, `curate` |
| `features/charge.py` | `Charge.calculate` |
| `features/frustration.py` | `Frustration.calculate_frustration` |
| `features/sasa.py` | `SASA.calculate` |

### 1.3 Spelling and grammar (rule 4)

No information was added. Corrections applied: *enure* → *ensure*, *feauture* → *feature*,
*Defualt* → *Default*, *caclulate* → *calculate*, *moleucle* → *molecule*, *withing* →
*within*, *previoulsy* → *previously*, *seqeuence* → *sequence*, *encompasing* →
*encompassing*, *preprocesing* → *preprocessing*, *Molcule* → *Molecule*, *analyis* →
*analysis*, *isnt* → *is not*, *wont* → *will not*, *an parameter* → *a parameter*,
*the the protein* → *the protein*, *duplicates the the feature space* → *duplicates the
feature space*, *you you will* → *you will*, *outdirectory to used to* → *output directory
used to*, *the number individual functions* → *the number of individual functions*, *in
series of parallel* → *in series or in parallel*, *kept to updated* → *kept to update*,
*to workout* → *to work out*, *structure as obtained* → *structure was obtained*.

Four copy-and-paste errors were also corrected, because they named the wrong quantity:

| File | Was | Is now |
| --- | --- | --- |
| `features/evolution.py` | "Initialise the **Charge** class" | "Initialise the **Evolution** class" |
| `features/rmsf.py` | "Initialise the **Charge** class"; `:param df_prot:`; ":returns: … **seqcharge** output" | "Initialise the **RMSF** class"; `:param df_proteins:` (the actual argument name); ":returns: … **rmsf** output" |
| `features/aev.py`, `charge.py`, `flexibility.py`, `frustration.py`, `structure.py`, `feature.py` | ":param path: The path of the pdb file that **DAS**/**SASA** is being calculated for" | the feature the class actually measures |
| `measure.py` | `:param features:` in `_setup_measures` | `:param features_list:` (the actual argument name) |

### 1.4 Feature class preambles (rule 5)

The "Class to house the different methods for calculating …" preamble was removed from all
thirteen feature classes, leaving the substantive description as the summary sentence.
`PROPKA` and `PKAANI` would otherwise have carried identical descriptions, so each now
names the program it calls (PROPKA3 and pKaANI respectively), which is stated in the
`calculate()` docstring of each.

The same preamble pattern (*"Class to handle …"*, *"Class encompassing …"*, *"Class for
taking …"*) is still present on `Measure`, `Aggregation`, `Preprocessing`, `PDB` and
`Uniprot`. Those are not feature classes, so they were left alone; say the word if you
want the same treatment there.

### 1.5 One code change

`src/coolpackage/features/__init__.py` resolved `feature_folder` as the literal relative path
`'src/features'`, i.e. against the working directory. Sphinx builds from `docs/`, so
`os.listdir` raised `FileNotFoundError` and `import measure` failed, leaving the Measure
page empty. It is now resolved from `__file__`, which also makes the package importable
from any working directory at runtime.

---

## 2. Flagged, not applied

### 2.1 Parameters with an empty description

Three parameters have a `:param:` entry with nothing in it.

| Location | Suggested text |
| --- | --- |
| `analysis.py:221` `Analysis.get_outliers`, `whis` | `:param whis: Multiplier applied to the interquartile range when setting the outlier bounds, which are the 25th percentile minus ``whis`` times the IQR and the 75th percentile plus ``whis`` times the IQR. The default is 1.5.` |
| `preprocessing.py:66` `Preprocessing.clean`, `df` | `:param df: The dataframe of measurements to clean.` |
| `protein.py:208` `PDB._curate_row`, `row_details` | `:param row_details: One row of the input dataframe. When ``PDB_only`` is False it unpacks to ``[uniprot_code, pdb_code, method_obtained, resolution, chains]``, when it is True it holds the PDB code alone.` |

### 2.2 Parameters that are documented but do not exist

| Location | Note |
| --- | --- |
| `features/evolution.py:23` `Evolution.__init__` | Documents `df_prot`, which is not in the signature (`include_modified`, `aa_properties`, `error_filename`). Either the argument is missing from the signature or the field should be deleted. I could not tell which from the code, so nothing was changed. |
| `patcher.py:181` `_trim_align` | Documents `fbasename`, but the signature is `(tmp_folder, align_file)`. Note also that `align_file` is overwritten on the first line of the body with `f'{tmp_folder}{os.sep}alignment.seg.ali'`, so the argument passed in is ignored. Suggested, once the code is settled: `:param align_file: Path to the alignment file to trim. Currently ignored: the function always reads ``alignment.seg.ali`` from ``tmp_folder``.` |

### 2.3 Undocumented parameters

| Location | Suggested text |
| --- | --- |
| `aggregation.py:406` `_aggregate_choose`, `df_stats` | `:param df_stats: The dataframe of per-lysine statistics produced by ``_calculate_statistics()``.`<br>`:type df_stats: pandas.DataFrame` |
| `analysis.py:87` `GO_search_term`, `df` / `code` / `name` | `:param df: Dataframe to take the subset of rows from.`<br>`:type df: pandas.DataFrame`<br>`:param code: GO term code to search for. If empty, it is looked up from ``name``.`<br>`:type code: str`<br>`:param name: GO term name to search for. At least one of ``code`` and ``name`` must be given.`<br>`:type name: str` |
| `analysis.py:108` `GO_search_protein`, `uniprot_entry` | `:param uniprot_entry: Uniprot code to list the GO terms of.`<br>`:type uniprot_entry: str` |
| `analysis.py:119` `_GO_get_data`, `uniprot_code` / `lock` / `index` / `total` | `:param uniprot_code: Uniprot code to download the GO terms for.`<br>`:type uniprot_code: str`<br>`:param lock: Lock used to serialise the writes to ``self.GO_dict`` and ``self.name_to_code``.`<br>`:type lock: threading.Lock`<br>`:param index: Position of this code in the list being processed, used for the progress message.`<br>`:type index: int`<br>`:param total: Total number of codes being processed, used for the progress message.`<br>`:type total: int` |
| `analysis.py:164` `plot_graph`, `plot_type` / `feature` / `uniprot_entry` / `resid` | `:param plot_type: Type of plot to draw, either 'histogram' or 'boxplot'.`<br>`:type plot_type: str`<br>`:param feature: Name of the column to plot.`<br>`:type feature: str`<br>`:param uniprot_entry: Uniprot code to restrict the plot to. False plots the whole dataframe.`<br>`:type uniprot_entry: str, optional`<br>`:param resid: Residue number to restrict the plot to, used together with ``uniprot_entry``.`<br>`:type resid: int, optional` |
| `analysis.py:252` `remove_df`, `df_to_remove` | `:param df_to_remove: Dataframe whose rows, matched by index, are dropped from ``self.df``.`<br>`:type df_to_remove: pandas.DataFrame` |
| `analysis.py:549` `get_contingency_table`, `GO_code` / `my_list` / `reference` | `:param GO_code: GO term the table is built for.`<br>`:type GO_code: str`<br>`:param my_list: Uniprot codes of the selection of interest.`<br>`:type my_list: list`<br>`:param reference: Uniprot codes of the background set the selection is compared against.`<br>`:type reference: list` |
| `analysis.py:572` `enrichment_analysis`, `feature_one` / `feature_two` / `uniprot_cnt_cutoff` | `:param feature_one: Three elements, the feature name and the lower and upper bounds of the region to select on. The default is ``['propka', 7, 11]``.`<br>`:type feature_one: list`<br>`:param feature_two: As ``feature_one``, for the second axis of the region. The default is ``['sasa', 0, 10]``.`<br>`:type feature_two: list`<br>`:param uniprot_cnt_cutoff: Minimum number of Uniprot codes a GO term must be associated with to be tested. The default is 1.`<br>`:type uniprot_cnt_cutoff: int` |
| `measure.py:1364` `_mv_files` (nested in `_cleanup_calculation_files`) | `:param files: Files to move.`<br>`:type files: list`<br>`:param dest: Directory, relative to ``self.outdir``, the files are moved into. It is created if missing.`<br>`:type dest: str` |
| `patcher.py:181` `_trim_align`, `align_file` | see 2.2 |
| `preprocessing.py:90` `calculate_vif`, `multi_vif` | `:param multi_vif: Set to True when the method is called repeatedly from ``calculate_diff_features()``, in which case the AEV columns have already been expanded and the constant column already added, so that preparation step is skipped. The default is False.`<br>`:type multi_vif: bool` |
| `features/pkaani.py:13` `PKAANI.__init__`, `calc_method` | The argument is accepted but never used anywhere in the class. Either remove it from the signature or document what it is meant to select. Nothing was changed, as the intent is not recoverable from the code. |

### 2.4 Undocumented return values

`:returns:` and `:rtype:` are missing from the following. Suggested text, one line each:

| Location | Suggested text |
| --- | --- |
| `aggregation.py:114` `aggregate_data` | `:returns: The aggregated dataframe, one row per lysine residue.`<br>`:rtype: pandas.DataFrame` |
| `aggregation.py:343` `_calculate_statistics` | `:returns: Dataframe carrying every statistic that can be used for aggregation.`<br>`:rtype: pandas.DataFrame` |
| `aggregation.py:406` `_aggregate_choose` | `:returns: ``df_stats`` with the feature columns the user did not select dropped.`<br>`:rtype: pandas.DataFrame` |
| `aggregation.py:460` `_aggregate_avg_less_avgaev` | `:returns: The aggregated dataframe, with the average AEV of a general lysine subtracted from the AEV columns.`<br>`:rtype: pandas.DataFrame` |
| `analysis.py:87` `GO_search_term` | `:returns: The rows of ``df`` whose Uniprot code is associated with the GO term, or a message string if the input is insufficient or the code and name do not match.`<br>`:rtype: pandas.DataFrame or str` |
| `analysis.py:108` `GO_search_protein` | `:returns: Names of the GO terms associated with the Uniprot code.`<br>`:rtype: list` |
| `analysis.py:221` `get_outliers` | `:returns: The rows of the selection that fall outside the interquartile bounds.`<br>`:rtype: pandas.DataFrame` |
| `analysis.py:356` `remove_not_important_residues` | `:returns: The dataframe reduced to the residues listed in ``req_resid_table``.`<br>`:rtype: pandas.DataFrame` |
| `analysis.py:549` `get_contingency_table` | `:returns: The 2x2 contingency table ``[[in term and in list, in term not in list], [not in term but in list, in neither]]``.`<br>`:rtype: list` |
| `analysis.py:572` `enrichment_analysis` | `:returns: One row per GO term, with columns 'GO ID', 'GO Term', 'raw p value', 'FDR', 'num in the region' and 'num in the bkgd', sorted by p-value.`<br>`:rtype: pandas.DataFrame` |
| `measure.py:349` `_match_resid_codes` | `:returns: The ``aa_properties`` dictionary matching the 3 letter code given.`<br>`:rtype: dict` |
| `measure.py:430` `_setup_report_errors_file` | `:returns: Name of the error file that was created.`<br>`:rtype: str` |
| `measure.py:554` `measure_dataframe` | `:returns: The measurements for every curated structure, or the string 'Call PDB_only method instead' when the class was set up with PDB_only.`<br>`:rtype: pandas.DataFrame or str` |
| `measure.py:645` `_measure_file` | `:returns: The measurements for the single file given, empty if none could be calculated.`<br>`:rtype: pandas.DataFrame` |
| `measure.py:880` `restart_measure` | `:returns: The string 'restart_measure() function not callable when using PDB_only' when the class was set up with PDB_only, otherwise nothing.`<br>`:rtype: str or None` |
| `measure.py:933` `_combine_dataframes` | `:returns: ``target`` with the values of ``to_merge`` inserted in column ``col_name``.`<br>`:rtype: pandas.DataFrame` |
| `measure.py:991` `measure_PDB_only` | `:returns: The string 'You have called the wrong method for measuring data, call the general measures function instead' when the class was not set up with PDB_only, otherwise nothing.`<br>`:rtype: str or None` |
| `measure.py:1123` `restart_measure_pdb_only` | `:returns: An explanatory string when the method cannot run, either because PDB_only is not set or because no previous measurement was found in the log file. Otherwise nothing.`<br>`:rtype: str or None` |
| `patcher.py:121` `_full_align` | `:returns: The PDB name of the target sequence, read from the alignment file.`<br>`:rtype: str` |
| `patcher.py:366` `fragment` | `:returns: One row per chain, each holding the number of gaps, the number of missing residues and the largest sequence gap.`<br>`:rtype: numpy.ndarray` |
| `protein.py:145` `gather_proteins` | `:returns: The exception raised while converting a list of PDB codes into a dataframe, if one was raised. Otherwise nothing.`<br>`:rtype: Exception or None` |
| `protein.py:208` `_curate_row` | `:returns: A single-row dataframe describing the curated structure, or None if the curation failed.`<br>`:rtype: pandas.DataFrame or None` |
| `protein.py:388` `clean_and_split_pdb` | `:returns: The largest sequence gap found across the conformers of the structure.`<br>`:rtype: int` |
| `features/evolution.py:91` `Evolution.calculate` | `:returns: Dataframe with the chain, residue number and ESM vector of every residue of interest, empty if the model could not be loaded.`<br>`:rtype: pandas.DataFrame` |
| `features/nmr.py:70` `NMR.calculate_legolas` | `:returns: Dataframe with the chain, residue number and legolas 15N shift of every residue of interest, empty if the calculation failed.`<br>`:rtype: pandas.DataFrame` |

### 2.5 Missing docstrings

**Modules.** None of the following carry a module docstring, so the module pages have no
introductory text: `aggregation.py`, `alphafold.py`, `helper.py`, `measure.py`,
`preprocessing.py`, `protein.py`, `uniprot.py`, `viewer.py`, and every file in
`src/features/` except `__init__.py` and `error_reporting.py`. For the feature files the
class docstring already carries the description, so a module docstring would only repeat
it. The ones that would earn their place:

- `alphafold.py`: "Download of AlphaFold structures and extraction of their PLDDT scores."
- `helper.py`: "Small utilities shared across the pipeline."
- `patcher.py` already has one and can serve as the model.

**Classes.**

| Location | Note |
| --- | --- |
| `analysis.py:44` `Analysis` | No class docstring. The `__init__` docstring describes it ("allows you to create graph and go over other metrics such as GO terms"); a class-level summary would let the module todos and the class description sit apart. |
| `viewer.py:9` `Viewer` | No class docstring, and no docstring anywhere else in the file. |
| `helper.py:23` `ShutUp` | No class docstring. Suggested: "Context manager redirecting stdout to os.devnull, used to silence the output of external programs." |

**Functions.**

| Location | Note |
| --- | --- |
| `analysis.py:77` `get_data`, `analysis.py:82` `get_data_alphafold` | Both carry a `# GW 05.12.24 function potentially unused - remove?` comment. Worth resolving before documenting them. |
| `analysis.py:248` `get_extreme_values` | Returns the rows where ``feature`` is below ``lower`` or above ``upper``. Note the defaults (1 and 14) are pKa-specific. |
| `helper.py:28` `ShutUp.__exit__` | Dunder, no docstring needed. |
| `uniprot.py:80` `_get_next_link`, `uniprot.py:86` `_get_batch` | Nested helpers of `get_organism_proteins`, paging through the Uniprot REST API. Not rendered by Sphinx. |
| `viewer.py` (10 nested callbacks) | Widget callbacks inside `Viewer.__init__` and `Viewer.advanced_plot`. Not rendered by Sphinx. |

### 2.6 Sentences that stop mid-way

Three descriptions are truncated in the source. The missing text cannot be recovered from
the code, so nothing was changed.

| Location | Text |
| --- | --- |
| `measure.py` `measure_dataframe`, Method | "… Iterate over the list of the files, check if structure file is" |
| `features/frustration.py` `calculate_frustration`, Method | "… Use this model to calculate the " |
| `protein.py` `PDB.clean`, last bullet | "check for modified residues in structure, convert back to" |

### 2.7 One inconsistency found while writing the feature table

`Measure._setup_measures()` expands `features=['all']` to

```
['propka', 'pkaANI', 'sasa', 'depth', 'aev', 'seqcharge', 'legolas', 'melodia',
 'aev_legolas', 'frustration', 'density', 'das', 'flexibility', 'esm', 'rmsf']
```

`'esm'` is not handled by any branch of the loop, and it is not the name of a class exported
by `src.features` (the ESM class is called `Evolution`, and its keyword is `'evolution'`).
It therefore falls through to the final `else` and raises `Exception("measure esm unknown")`,
so `features=['all']` cannot currently run to completion. Changing `'esm'` to `'evolution'`
in that list would fix it. Nothing was changed, as this is a code fix rather than a
documentation one.

The `:param features:` list in `Measure.__init__` was corrected as part of the pass: it
previously omitted `'flexibility'`, `'evolution'`, `'rmsf'` and `'aev_legolas'`, and it did
not mention that any class dropped into the features folder can be requested by its class
name. The feature table on the "Feature measurements" page is built from the same mapping.
