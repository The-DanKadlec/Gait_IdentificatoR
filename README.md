# Gait_IdentificatoR
Automated pipeline for continuous gait-event identification from trunk-mounted GPS resultant acceleration during overground running.

# Trunk-accelerometer step and stride interval detection

Detects gait impact landmarks from a single trunk-mounted accelerometer
(resultant signal) during steady-state running and exports step and/or
stride interval time-series. This is the gait-event identification pipeline described in [DOI TBC]

## Input format
A folder of CSV files. Each file needs two columns:
- a time column in seconds (default name `time_s`)
- a resultant-acceleration column (default name `resultant`)

Column names and sampling frequency are set in the USER SETTINGS block.

## Usage
1. Open the script and edit the USER SETTINGS block: set `input_dir`,
   your column names, `fs`, and `export_type` ("step", "stride", or "both").
2. Run the script. It creates an `intervals_output` sub-folder.
3. Each input file produces `<name>_step.csv` and/or `<name>_stride.csv`,
   plus `_run_summary.csv` listing landmark counts and any gaps.

## Output
Each output CSV has one row per landmark: the landmark index, its time
(s), and the step or stride interval (s). Missing intervals are kept as
explicit `NA` in place and never interpolated, so the series keeps its
temporal position. The run summary and console warnings report the number
and location of any gaps for manual inspection.

## Parameters
All processing parameters (filter band, physiological interval limits,
Viterbi prior strength, sinc refinement, gap recovery) are documented in
the USER SETTINGS block and can be adjusted. Defaults match the published
pipeline. Set `max_insert <- 0` to disable gap recovery.

## Example
An example recording (`example_GPS.csv`) is included. With `input_dir`
set to its folder and default settings, running the script reproduces a
complete step and stride series.

## Citation
[Manuscript in preparation]

## Licence
MIT License
