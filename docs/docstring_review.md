# Documentation review — `src/`

Companion to the Sphinx/docstring pass. Sections 1 and 2 record what was changed;
sections 3 to 6 are **flags only**, nothing in them has been applied to the code.

Line numbers refer to the files as they stand after the reformatting pass.

---

## 1. Sphinx configuration (applied)

The move of the sources from the repository root and `features/` into `src/` and
`src/features/` broke module resolution in three separate places. All are fixed.

| File | Change |
| --- | --- |
| `docs/source/conf.py` | `sys.path` now carries both the repository root (the modules import each other as `src.<module>`) and `src/` (the `autodoc` directives in the `.rst` files address the modules by bare name, e.g. `protein.PDB`). |
| `docs/source/conf.py` | `autoapi_dirs` pointed at the deleted `../../features`; it now points at `src/features`. |
| `docs/source/conf.py` | `autodoc_mock_imports` extended from `["modeller"]` to also cover `nglview`, `torchani`, `esm`, `melodia_py`, `frustratometer` and `pkaani`. Without these, `viewer.py` and the GPU/optional feature modules cannot be imported on a machine that only has the core dependencies, and the corresponding pages come out empty. |
| `docs/source/conf.py` | `todo_include_todos = True` added, so the `.. todo::` notes described in section 4 stay visible in the rendered pages instead of being silently dropped. |
| `src/features/__init__.py` | `feature_folder` was the literal relative path `'src/features'`, resolved against the current working directory. Sphinx builds from `docs/`, so `os.listdir` raised `FileNotFoundError` and `import src.measure` failed. It is now resolved from `__file__`. This also makes the package importable from any working directory at runtime, not only from the repository root. |

New pages were added for the modules that the move left with no coverage at all:
`aggregation.rst`, `preprocessing.rst`, `patcher.rst`, `viewer.rst`, `helper.rst`.
An empty `docs/source/_static/` was created to clear the one remaining build warning.

`index.rst` was then reorganised into one captioned toctree per pipeline step
(1. Uniprot, 2. Protein structure processing, 3. Measure, 4. Analysis,
5. Viewer, plus an Appendix for `helper` and the FAQ). Steps 4 and 5 carry a
"work in progress" marker in their caption and a `.. note::` at the top of the
page.

### The "Feature measurements" table

The same commit that moved the sources (`81a4711`) also added
`src/features/__init__.py`, which `features/` did not previously have. That one
file changed what AutoAPI reports: with no `__init__.py` each feature file was a
*top-level module*, so the landing-page template's
`pages|selectattr("is_top_level_object")` yielded one row per feature. With the
`__init__.py` in place `features` is a package and is itself the only top-level
object, so the table collapsed to a single row pointing at an intermediate
package page that lists the features without a table.

Two changes restore it:

* `_templates/autoapi/index.rst` now iterates
  `pages|selectattr("type", "equalto", "module")`, i.e. the submodules of the
  package, rather than top-level objects. Class names are rendered with a
  `~` prefix so the table reads `AEV` rather than `features.aev.AEV`.
* `exclude_patterns = ['autoapi/features/index.rst']` in `conf.py` drops the
  intermediate package page, so the table links straight to each feature and no
  document ends up in two toctrees.

The table is back to 15 rows, one per feature class, each linking to its own
class anchor.

Verified: `sphinx-build -b html source build` completes with **0 warnings**, and every
module under `src/` imports successfully under the mocks.

**Not changed, but worth a decision:** `project`, `htmlhelp_basename` and the
`index.rst` title are all still the `sphinx-quickstart` placeholder `coolpackage`.

---

## 2. Docstring format (applied)

All **144** docstrings across the 26 modules were converted from NumPy-style section
headers to reST field lists, which is what Sphinx parses natively (`napoleon` is not
among the enabled extensions, so the previous format was being rendered as
undifferentiated prose).

The mapping used throughout:

| Was | Is now |
| --- | --- |
| `Parameters` / `----------` block | `:param name:` + `:type name:` |
| `Returns` / `-------` block | `:returns:` + `:rtype:` |
| `Method` / `------` block | `**Method**` paragraph |
| `Example` / `-------` block | `**Example**` followed by a `::` literal block |
| `name -> type` (used in ~half the files) | `:type name: type` |

Two mechanical problems the old format was hiding, now fixed by the conversion:

* Roughly half of the parameter entries used `name -> type` rather than the
  NumPy `name : type`. Even under `napoleon` those would not have been recognised.
* Option lists (`- 'avg': ...` in `Aggregation.__init__`, `- AtomDepth ...` in
  `Depth.__init__`, the `Method` bullet lists in `Preprocessing.undersampling`)
  had no blank line before them and would have been folded into the preceding
  paragraph. They now render as lists.

Example output blocks are wrapped in literal blocks, so the tabular DataFrame
dumps keep their alignment instead of being reflowed.

---

## 3. Undocumented parameters and return values

Nothing below has been added. Suggested text is offered for each.

### 3.1 `src/viewer.py` — no documentation at all

The module has 12 functions and one class, none with a docstring. This is the
largest single gap; the new `viewer.rst` page currently renders an empty class.

* **`Viewer` (line 9)** — no class docstring.
  > Class providing an interactive Jupyter interface over an :class:`Analysis`
  > instance. Couples a plotly scatter of two features against `ipywidgets`
  > range sliders and an `nglview` structure viewer, so that a region of the
  > feature space can be selected, inspected residue by residue and exported.

* **`Viewer.__init__(self, analysis, outdir='result')` (line 11)**
  > Initialise the Viewer class.
  >
  > `:param analysis:` The Analysis instance whose dataframe (`analysis.df`) is
  > to be displayed.
  > `:type analysis: analysis.Analysis`
  > `:param outdir:` The directory the regional data export is written into.
  > Default is `result`.
  > `:type outdir: str`

* **`Viewer.advanced_plot(self, export_path='Regional_Data.csv', cutoff=0)` (line 170)**
  > Build the interactive widget assembly: pKa/SASA range sliders, a GO term
  > dropdown, an export button and a structure viewer, and return it for display
  > in a notebook. The assembled widget is also stored on `self.plot`.
  >
  > `:param export_path:` Name of the csv file, written inside `self.outdir`,
  > that the selected region is exported to. Default `Regional_Data.csv`.
  > `:type export_path: str`
  > `:param cutoff:` Minimum number of UniProt entries a GO term must be
  > associated with to appear in the dropdown.
  > `:type cutoff: int`
  > `:returns:` The assembled widget container.
  > `:rtype: ipywidgets.VBox`

* Nested callbacks — `call_back_buttom_export` (72), `call_back_buttom` (99),
  `call_back_buttom_open_url` (123), `interact_slides` (174),
  `interact_dropdown` (250), `interact_3D` (311), `update_point_2` (329),
  `update_point_3` (367), `call_back_buttom_2` (385). These are widget handlers
  and are not part of the public API; a one-line docstring each is enough, e.g.
  > Callback for the export button: write the currently selected region to
  > `self.export_path`.
  >
  > `:param b_export:` The button that triggered the callback (unused).
  > `:type b_export: ipywidgets.Button`

  Note the recurring misspelling **`buttom`** in these identifiers. It is a name,
  not prose, so it was left alone, but it is worth a rename.

### 3.2 `src/analysis.py`

* **`Analysis` (line 31)** — no class docstring. This class *is* on a
  hand-written documentation page (`analysis.rst`), so the page currently opens
  with nothing.
  > Class collecting the post-measurement analysis tools: plotting of individual
  > features, retrieval and enrichment analysis of GO terms from UniProt,
  > outlier extraction, and merging or trimming of measurement dataframes.

* **`get_data(self, uniprot_entry, resid)` (line 64)** — no docstring, return
  undocumented. Marked `# GW 05.12.24 function potentially unused - remove?`
  at line 63; see section 4.
  > Return the subset of the measurements dataframe matching one UniProt entry
  > and one residue number.
  >
  > `:param uniprot_entry:` Uniprot code of interest.
  > `:type uniprot_entry: str`
  > `:param resid:` Residue number of interest.
  > `:type resid: int`
  > `:returns:` The matching rows of `self.df`.
  > `:rtype: pandas.DataFrame`

* **`get_data_alphafold(self)` (line 69)** — no docstring, return undocumented.
  Also marked potentially unused.
  > Return the subset of the measurements dataframe obtained from predicted
  > (AlphaFold) structures, i.e. the rows whose `Method` column is `Predicted`.
  >
  > `:returns:` The matching rows of `self.df`.
  > `:rtype: pandas.DataFrame`

* **`GO_search_term(self, df, code='', name='')` (line 74)** — all three
  parameters and the return undocumented.
  > `:param df:` Dataframe to take the subset from.
  > `:type df: pandas.DataFrame`
  > `:param code:` GO term code to search for. Either this or `name` must be given.
  > `:type code: str`
  > `:param name:` GO term name to search for, resolved to a code through
  > `self.name_to_code`. Either this or `code` must be given.
  > `:type name: str`
  > `:returns:` The rows of `df` whose Uniprot entry is associated with the GO
  > term. If neither `code` nor `name` is given, or if the two disagree, an
  > explanatory string is returned instead of a dataframe.
  > `:rtype: pandas.DataFrame or str`

  The mixed return type (dataframe on success, message string on bad input) is
  worth reconsidering; raising `ValueError` would be easier for a caller to handle.

* **`GO_search_protein(self, uniprot_entry)` (line 95)**
  > `:param uniprot_entry:` Uniprot code to list the GO terms of.
  > `:type uniprot_entry: str`
  > `:returns:` The names of the GO terms associated with that entry.
  > `:rtype: list`

* **`_GO_get_data(self, uniprot_code, lock, index, total)` (line 106)**
  > `:param uniprot_code:` Uniprot code to fetch GO annotations for.
  > `:type uniprot_code: str`
  > `:param lock:` Lock guarding the shared `self.GO_dict` and
  > `self.name_to_code` dictionaries against concurrent writes.
  > `:type lock: threading.Lock`
  > `:param index:` Position of this code in the overall list, used for the
  > progress message only.
  > `:type index: int`
  > `:param total:` Total number of codes being processed, used for the
  > progress message only.
  > `:type total: int`

* **`GO_get_data(self)` (line 141)** — no docstring.
  > Retrieve the GO annotations for every unique Uniprot entry in the
  > measurements dataframe, in parallel over a thread pool, and populate
  > `self.GO_dict`, `self.name_to_code` and `self.code_to_name`.

* **`plot_graph(self, plot_type, feature, uniprot_entry=False, resid=False)` (line 151)**
  > `:param plot_type:` Type of plot to draw, either `histogram` or `boxplot`.
  > `:type plot_type: str`
  > `:param feature:` Name of the measurements column to plot.
  > `:type feature: str`
  > `:param uniprot_entry:` If given together with `resid`, restrict the plot to
  > that single residue. Default False plots the whole dataframe.
  > `:type uniprot_entry: str, optional`
  > `:param resid:` Residue number, used together with `uniprot_entry`.
  > `:type resid: int, optional`

* **`get_outliers(...)` (line 208)** — `whis` is named but has no description,
  and the return is undocumented.
  > `:param whis:` Multiple of the interquartile range beyond the first and
  > third quartiles past which a point is treated as an outlier. Default 1.5.
  > `:type whis: float`
  > `:returns:` The rows lying outside those bounds.
  > `:rtype: pandas.DataFrame`

* **`get_extreme_values(self, feature, lower=1, upper=14)` (line 235)** — no docstring.
  > Extract the rows whose value for a given feature falls outside a fixed range.
  > Unlike :meth:`get_outliers` the bounds are absolute, not derived from the
  > distribution.
  >
  > `:param feature:` Name of the measurements column to threshold.
  > `:type feature: str`
  > `:param lower:` Values strictly below this are returned. Default 1.
  > `:type lower: float`
  > `:param upper:` Values strictly above this are returned. Default 14.
  > `:type upper: float`
  > `:returns:` The rows outside the range.
  > `:rtype: pandas.DataFrame`

  The defaults `1` and `14` are pKa-specific; the parameter name is generic.

* **`remove_df(self, df_to_remove)` (line 239)** — the docstring is present but
  **empty** (whitespace only). See section 4.
  > Drop from the measurements dataframe every row present in the given
  > dataframe, matched by index, and report the row counts before and after.
  >
  > `:param df_to_remove:` Dataframe whose index labels are to be dropped from
  > `self.df`.
  > `:type df_to_remove: pandas.DataFrame`

* **`remove_not_important_residues(...)` (line 344)** — return undocumented.
  > `:returns:` The trimmed measurements dataframe, also stored on `self.df`
  > and written to `outname` in the output directory.
  > `:rtype: pandas.DataFrame`

* **`get_contingency_table(self, GO_code, my_list, reference)` (line 539)**
  > `:param GO_code:` The GO term the table is built for.
  > `:type GO_code: str`
  > `:param my_list:` The Uniprot codes of the selected region.
  > `:type my_list: list`
  > `:param reference:` The Uniprot codes of the background set.
  > `:type reference: list`
  > `:returns:` A 2x2 contingency table, `[[in term and in list, in term and not
  > in list], [not in term and in list, not in term and not in list]]`.
  > `:rtype: list`

* **`enrichment_analysis(...)` (line 562)** — all parameters and the return
  undocumented.
  > `:param feature_one:` Three-element specification of the first axis of the
  > region, `[feature name, lower bound, upper bound]`. Default
  > `['propka', 7, 11]`.
  > `:type feature_one: list`
  > `:param feature_two:` Same for the second axis. Default `['sasa', 0, 10]`.
  > `:type feature_two: list`
  > `:param uniprot_cnt_cutoff:` Minimum number of Uniprot entries a GO term
  > must be associated with to be tested. Default 1.
  > `:type uniprot_cnt_cutoff: int`
  > `:returns:` One row per GO term with its raw Fisher p value, the
  > Benjamini-Hochberg adjusted p value, and the counts in the region and in
  > the background. Columns: `GO ID`, `GO Term`, `raw p value`, `FDR`,
  > `num in the region`, `num in the bkgd`.
  > `:rtype: pandas.DataFrame`

### 3.3 `src/measure.py`

Every one of the following returns a **sentinel string** on the wrong-mode path
rather than raising, which is exactly why documenting the return matters here: a
caller assigning the result gets a `str` instead of a `DataFrame` with no error.

* **`_match_resid_codes` (line 352)**
  > `:returns:` `aa_properties` -- the 3 letter codes and atom names for the
  > requested residue, in the form consumed by the feature classes.
  > `:rtype: dict`
* **`_setup_report_errors_file` (line 433)**
  > `:returns:` `new_file_name` -- the name of the error file created for this
  > run, of the form `measure_errors_{date}_{x}.txt`.
  > `:rtype: str`
* **`measure_dataframe_parallel` (line 561)** and **`measure_dataframe` (line 655)**
  > `:returns:` Nothing on the normal path; the results are stored on `self.df`.
  > If the class was set up for a different mode, an explanatory string is
  > returned instead (`'Call PDB_only method instead'`, `'Series measurement
  > method called, call measure_dataframe_parallel instead'`).
  > `:rtype: str or None`
* **`_measure_file` (line 807)**
  > `:returns:` The measurements for that one file, one row per residue of
  > interest. Empty, with only the identifying columns, if the file yielded
  > nothing.
  > `:rtype: pandas.DataFrame`
* **`restart_measure` (line 1046)**, **`measure_PDB_only` (line 1156)**,
  **`restart_measure_pdb_only` (line 1288)** — same sentinel-string pattern as
  above; document as `:rtype: str or None`.
* **`_combine_dataframes` (line 1098)**
  > `:returns:` `target` -- the target dataframe with the merged column added.
  > `:rtype: pandas.DataFrame`
* **`_setup_measures` (line 180)** — the parameter is called `features_list` in
  the signature but the docstring documents **`features`**. One of the two names
  needs to change; left as the author wrote it (see section 6).
* **`_cleanup_calculation_files._mv_files(files, dest)` (line 1531)**
  > `:param files:` Paths of the leftover files to move.
  > `:type files: list`
  > `:param dest:` Directory, relative to `self.outdir`, to move them into.
  > `:type dest: str`

### 3.4 `src/aggregation.py`

* **`aggregate_data` (line 117)**
  > `:returns:` `self.df_agg` -- the aggregated dataframe, one row per residue.
  > `:rtype: pandas.DataFrame`
* **`_calculate_statistics` (line 346)**
  > `:returns:` `df_stats` -- one row per residue carrying every statistic
  > (min, max, median, average, standard deviation, range, random draw) for
  > every feature.
  > `:rtype: pandas.DataFrame`
* **`_aggregate_choose(df_stats)` (line 409)** — parameter and return undocumented.
  > `:param df_stats:` The full statistics dataframe produced by
  > :meth:`_calculate_statistics`, from which the user's choices are kept.
  > `:type df_stats: pandas.DataFrame`
  > `:returns:` `df_stats` reduced to the statistics the user selected.
  > `:rtype: pandas.DataFrame`
* **`_aggregate_avg_less_avgaev` (line 463)**
  > `:returns:` `df_data_agg` -- the aggregated dataframe with the average
  > lysine AEV subtracted from the AEV columns.
  > `:rtype: pandas.DataFrame`

### 3.5 `src/patcher.py`

* **`_full_align` (line 123)** — return undocumented.
  > `:returns:` `seq_name` -- the PDB name of the target sequence read out of
  > the alignment.
  > `:rtype: str`
* **`_trim_align(tmp_folder, align_file)` (line 180)** — the docstring documents
  `fbasename`, which is **not** a parameter of this function; `align_file` is
  undocumented. Left as written (section 6). Suggested replacement:
  > `:param align_file:` Path to the alignment file to trim in place.
  > `:type align_file: str`
* **`fragment` (line 368)** — return undocumented.
  > `:returns:` One `[number of gaps, number of missing residues, largest
  > sequence gap]` triple per chain.
  > `:rtype: numpy.ndarray`

### 3.6 `src/protein.py`

* **`gather_proteins` (line 143)** — return undocumented, and the only return is
  an oddity: on failure to build the dataframe the function `return e`, i.e. it
  **returns the exception object** rather than raising it (line 165).
  > `:returns:` Nothing on success; the results are stored on `self.df`. If the
  > list of PDB codes could not be converted to a dataframe, the exception is
  > returned rather than raised.
  > `:rtype: Exception or None`

  Worth fixing rather than documenting: `raise` would be the expected behaviour.

* **`_curate_row` (line 206)** — `row_details` is named but has **no
  description**, and the return is undocumented. From the unpacking at lines
  229-232:
  > `:param row_details:` One row of the input dataframe as a list. When
  > `PDB_only` is False, `[uniprot_code, pdb_code, method_obtained, resolution,
  > chains]`; when True, a single-element list holding the PDB code.
  > `:type row_details: list`
  > `:returns:` A single-row dataframe describing the curated structure, or
  > None if the structure was skipped or could not be curated.
  > `:rtype: pandas.DataFrame or None`

### 3.7 `src/uniprot.py`

* **`get_organism_proteins._get_next_link(headers)` (line 85)** — no docstring.
  > Extract the URL of the next page of results from the response headers of a
  > paginated UniProt API call.
  >
  > `:param headers:` Response headers of the previous request.
  > `:type headers: requests.structures.CaseInsensitiveDict`
  > `:returns:` The next page URL, or None when the last page has been reached.
  > `:rtype: str or None`
* **`get_organism_proteins._get_batch(batch_url)` (line 91)** — no docstring.
  > Generator walking the paginated UniProt API, yielding one response and the
  > total result count per page until the pages are exhausted.
  >
  > `:param batch_url:` URL of the first page of results.
  > `:type batch_url: str`
  > `:yields:` `(response, total)` for each page.

### 3.8 `src/helper.py`

* **`ShutUp` (line 23)** — no class docstring.
  > Context manager redirecting standard output to the null device, used to
  > silence the third-party packages that print unconditionally.
* **`__enter__` (line 24)** and **`__exit__` (line 28)** — no docstrings. The
  class docstring above covers the behaviour; if entries are wanted:
  > Redirect `sys.stdout` to the null device, keeping a reference to the
  > original stream.
  >
  > Close the null device and restore the original `sys.stdout`.
  > `:param args:` Exception type, value and traceback, ignored.

### 3.9 `src/preprocessing.py`

* **`calculate_vif(data, features, multi_vif=False)` (line 92)** — `multi_vif`
  undocumented.
  > `:param multi_vif:` When False (default), the dataframe is first reduced to
  > `features`, a constant column is added and any `aev` column is expanded into
  > one column per AEV component with the all-zero components dropped. When
  > True, `data` is used as given, having already been prepared by a previous
  > call. Default False.
  > `:type multi_vif: bool`

### 3.10 `src/features/`

* **`pkaani.PKAANI.__init__` (line 13)** — `calc_method='propka'` undocumented.
  It is also **never read anywhere in the class**; it looks like a leftover from
  when PROPKA and pKaANI shared one module. Either remove it or document it as
  accepted and ignored.
* **`rmsf.RMSF.__init__` (line 14)** — the signature takes `df_proteins`, the
  docstring documents `df_prot`. See section 6.
* **`evolution.Evolution.__init__` (line 23)** — documents `df_prot`, which is
  **not a parameter of this function at all**. See section 6.
* **`evolution.Evolution.calculate` (line 91)** — return undocumented.
  > `:returns:` `df_esm` -- dataframe with information on chain, residue number
  > and the ESM embedding of the local sequence. Outline: `Chain  Resid
  > evolution`. Empty with those columns if the structure yielded no residues
  > of interest.
  > `:rtype: pandas.DataFrame`
* **`nmr.NMR.calculate_legolas` (line 70)** — return undocumented, although the
  Example block already shows the shape.
  > `:returns:` `df_legolas` -- dataframe with information on chain, residue
  > number and the predicted 15N chemical shift. Empty with those columns if
  > legolas produced no output.
  > `:rtype: pandas.DataFrame`

---

## 4. Incomplete comments and docstrings

### 4.1 Truncated sentences left mid-phrase

These are docstring text, not code comments, and each stops in the middle of a
clause. They were carried over verbatim.

| Location | Text |
| --- | --- |
| `src/features/frustration.py:62` | "Use this model to calculate the" |
| `src/protein.py:532` (`PDB.clean`) | "check for modified residues in structure, convert back to" |
| `src/measure.py:574` and `src/measure.py:668` | "Iterate over the list of the files, check if structure file is" |
| `src/aggregation.py:33-34` | "`average subtract aev`: Take the average of all aev features except the aevs **WHAT IS USED HERE**" |

### 4.2 Explicit "finish this" markers inside docstrings

`Measure.measure_dataframe_parallel` (line 561) and `Measure.measure_dataframe`
(line 655) both opened with a bare `TODO FINISH THIS` line. These are the two
main entry points of the class. They are now `.. todo::` directives so that
Sphinx renders them as flagged admonitions rather than as body prose, and
`todo_include_todos = True` keeps them visible.

`Measure._match_resid_codes` (line 352) carried two TODO lines in the middle of
its description; same treatment. Their content:

* default to carbamylation data for an unrecognised residue code should be
  replaced by stopping the run;
* checking which of PROPKA (ASP, GLU, HIS, CYS, TYR, LYS, ARG) and pKaANI (ASP,
  GLU, HIS, TYR, LYS) can actually run for the requested residue.

### 4.3 Empty docstring

`Analysis.remove_df` (`src/analysis.py:239`) has a triple-quoted docstring
containing only whitespace. This is worse than having none: Sphinx renders the
member with a blank description and `ast.get_docstring` returns an empty string,
so it passes "has a docstring" checks. Suggested text is in section 3.2.

### 4.4 TODO comments in code (not docstrings, left untouched)

| Location | Note |
| --- | --- |
| `src/analysis.py:14-20` | A four-item "TODO Section" header block: GO term analysis only works for propka not pkaani; the GO term functions may no longer work with later additions; the `dropna` on init was removed and needs replacing with a targeted cleaner. |
| `src/analysis.py:63, 68` | `get_data` and `get_data_alphafold` marked "function potentially unused - remove?" since 05.12.24. Neither is called anywhere in `src/`. |
| `src/aggregation.py:218` | "implement the autoencoder method here" — `autoencoder` is an advertised option of `aev_red_method` in the `Aggregation.__init__` docstring but is not implemented. |
| `src/features/charge.py:109` | ability to use other letter codes and charges. |
| `src/features/sasa.py:129` | element codes are LYS-specific, needs generalising. |
| `src/features/frustration.py:89` | `openmm` and `pdbfixer` still to be added to the README. |
| `src/measure.py:1346, 1376` | removing a measurement from `measures_log.txt`; switching to a "completed" column. |
| `src/patcher.py:118, 167, 309` | `env.io.two_char_chain = True` — "check locations of these to see if they do anything" (three identical copies). |
| `src/patcher.py:646` | "do we actually need to launch modeller if gaps are 0?" |
| `src/protein.py:121` | could `PDB_only` be inferred from the measures dataframe format. |
| `src/protein.py:427` | realignment only possible for structures with a Uniprot code. |
| `src/features/nmr.py:59` | Hard-coded developer path `/home/gweston/Documents/extra_packages/legolas-main/test/legolas.py`, with the comment "if you are not GW and running this, you will need to change this path". It is overridable via `LEGOLAS_PATH`, but the fallback should probably not be a personal path. |

---

## 5. Spelling corrections applied to docstrings

Only these were changed; no wording was otherwise altered and no information was
added or removed.

| File | Was | Now |
| --- | --- | --- |
| `src/features/feature.py` | `feauture` | feature |
| `src/features/feature.py` | `please enure that the class name if the same` | please ensure that the class name is the same |
| `src/features/sasa.py` | `moleucle` | molecule |
| `src/features/aev.py` | `For each NZ atom withing the lysines` | within |
| `src/features/nmr.py` | `within the the protein structure` | within the protein structure |
| `src/features/PROPKA.py` | `did not caclulate correctly` | calculate |
| `src/features/depth.py` | `PDBparser` | PDBParser |
| `src/aggregation.py` | `Defualt is set to 100` | Default |
| `src/aggregation.py` | `to workout which features show variation` | to work out |
| `src/preprocessing.py` | `Class encompasing methods` | encompassing |
| `src/preprocessing.py` | `used for preprocesing the data` | preprocessing |
| `src/analysis.py` | `remove anything from the dataframe that isnt in this list` | isn't |
| `src/analysis.py` | `List of features that the analyis should extract` | analysis |
| `src/protein.py` | `Finds all pdb files created previoulsy` | previously |
| `src/protein.py` | `alignment of the structure to the uniprot seqeuence` | sequence |
| `src/protein.py` | `Using the replacement dict given as an parameter` | as a parameter |
| `src/patcher.py` | `bb.Molcule instance` | `bb.Molecule` |
| `src/patcher.py` | `a gap bigger than cutoff which wont be patched` | won't |
| `src/measure.py` | `reducing the number individual functions` (x2) | the number of individual functions |
| `src/alphafold.py` | `The outdirectory to used to ...` (x2) | The output directory used to |
| `docs/source/index.rst` | `featurisation of aminoacids` | amino acids |

---

## 6. Wrong content, flagged but deliberately not changed

These are factual errors rather than typos, so correcting them would change what
the documentation claims. They need your decision.

### 6.1 Copy-paste class names in `__init__` summaries

Both of these say they initialise the **Charge** class:

* `src/features/evolution.py:30` — `Evolution.__init__`: "Initialise the Charge class"
* `src/features/rmsf.py:21` — `RMSF.__init__`: "Initialise the Charge class"

### 6.2 Copy-paste feature names in `calculate` parameter descriptions

Six `calculate` methods describe `path` as "The path of the pdb file that **DAS**
is being calculated for", inherited from `das.py`:

| File | Line | Should presumably read |
| --- | --- | --- |
| `src/features/charge.py` | 55 | seqcharge |
| `src/features/flexibility.py` | 52 | flexibility |
| `src/features/frustration.py` | 64 | frustration |
| `src/features/rmsf.py` | 61 | RMSF |
| `src/features/structure.py` | 69 | melodia |
| `src/features/feature.py` | 46 | the feature (this is the template class) |

`src/features/aev.py:90` has the same problem the other way round: it says
"The path of the pdb file that **SASA** is being calculated for" in `AEV.calculate`.

Two return descriptions also carry the wrong feature name:

* `src/features/rmsf.py:71` — `df_rmsf` described as "dataframe with information
  on chain, residue number and **seqcharge** output".
* `src/features/feature.py:52` — the template's return outline shows a column
  named `das` rather than a generic placeholder.

`src/features/feature.py` is the template new contributors copy, so its two
errors propagate.

### 6.3 Parameter names that do not match the signature

| Location | Signature | Docstring |
| --- | --- | --- |
| `src/measure.py:180` | `_setup_measures(self, features_list)` | documents `features` |
| `src/features/rmsf.py:14` | `__init__(self, df_proteins, ...)` | documents `df_prot` |
| `src/patcher.py:180` | `_trim_align(tmp_folder, align_file)` | documents `fbasename`, not `align_file` |
| `src/features/evolution.py:23` | `__init__(self, include_modified, aa_properties, error_filename)` | documents `df_prot`, which the function does not take at all |

In reST field-list form these now render as documented parameters that do not
exist in the rendered signature, which is more visible than it was under the old
format. Fixing the names is a one-line change in each case; I have left the
author's text so the change is yours to make.

### 6.4 Documented but unimplemented option

`Aggregation.__init__` lists `autoencoder` as a value of `aev_red_method`
("Uses an autoencoder to reduce the dimensions, better for non-linear data"),
but `src/aggregation.py:218` is `# TODO: implement the autoencoder method here`.
Selecting it currently does nothing.
