# LC-MS/MS Viewer V2
Analysis of targeted LC-MS/MS data recorded on a high-resolution QTOF (one 
MS1 scan should be followed by three MS2 scans with different collision energies,
e.g. CE=10, CE=20, CE=40):
import of .mzXML files , chromatogram extraction (EIC) for all target metabolites, 
examination of precursor purity, denoising of MS2 spectra, comparison of MS2 
fragments with predicted fragments from CFM-ID 4.0, and export of the result.

## What it is for

Confirmation of hits (metabolite accumulations in CRISPRi strains) from FI-MS 
screen via targeted LC-MS/MS on a high-resolution QTOF (see STAR protocol for 
further instructions). Each raw file targets one metabolite, in one CRISPRi 
knockdown strain, in one polarity.

- Was the correct precursor mass selected in the method?
- Are there peaks with another mass than the target within the isolation 
width of the quadrupole?
- How strong is the accumulation of a metabolite compared to all other samples
(calculation of fold-changes)?
- Do the MS2 fragments matches fragments form predicted spectra?

## Requirements

- MATLAB R2019b or newer
- Signal Processing Toolbox (`findpeaks`)

## Installation

1. Convert the raw data to `.mzXML` with MSConvert (see MSConvert_Settings.png)
2. Select folder with downloaded software.
3. Right-click this folder in MATLAB, choose
   *Add to Path > Selected Folders and Subfolders*.
4. Run the app: Open LC_MSMS_Viewer_V2.m with double-click. Press Run. GUI will open.
5. "Load Raw Data" -> Select .mzXML files -> Select excel sheet with hitList
   -> Select folder with predicted spectra from CFM-ID 4.0 -> Enter filename 
   for saving impirted data.


### Needed

| Input | Description |
|---|---|
| `.mzXML` files | named `<abbreviation>_<gene>_<positionBatch>_<polarity>.mzXML`, e.g. `3psme_aroC_P4C6msAV958_neg.mzXML` |
| hit list `.xlsx` | one row per targeted metabolite; columns `Gene`, `Metabolite`, `MetaboliteAbbreviation`, `Polarity`, `Mode`, `Mass` |
| predicted spectra | CFM-ID `.txt` files named `<abbreviation>_<polarity>.txt` |

See Example Data on Zenodo: https://zenodo.org/records/21990180


## Using the app

| Menu | What it does |
|---|---|
| **Load > Raw Data** | Run the whole pipeline and show the result |
| **Load > Analysed Data** | Reopen a data set saved earlier |
| **Export > To Workspace** | Copy the open data set to the base workspace as `LC_MSMS_data` |

The table lists every metabolite with its fold change (compared to the median 
of all samples and the two quality check-ups. Selecting a row draws:
- its extracted ion chromatogram against the median intensity of all samples 
(sample is shown in red; red dot shows timepoint where MS2 spectra are extracted;
 dashed black line shows the median intensity of all other samples)
- the MS1 spectrum of the scan with the highest intensity for the target ion 
within the isolation width of the quadrupole (+/- 0.65 m/z)
- MS2 spectra at the three different collision energies. Fragments matching 
predicted fragments are labelled with a stra ("*"). 


## Using it as a script

The whole analysis can run without the GUI:

```matlab
lc = run_lc_msms_pipeline( ...
    'C:\data\raw', ...
    'C:\data\hitList.xlsx', ...
    'C:\data\PredictedSpectra', ...
    'C:\results\myrun');
```

To change the method, edit a copy of the parameters rather than the source:

```matlab
params = lc_default_params();
params.ms2MinPeakHeight = 200;    % keep smaller fragments
params.eicToleranceDa   = 0.005;  % wider extraction window
lc = run_lc_msms_pipeline(raw, hits, predictions, out, params);
```

## Layout

```
LC_MSMS_Viewer_V2.m              the app
lc_default_params.m              all parameters

Functions/
  run_lc_msms_pipeline.m         runs steps 0 to 9
  read_mzxml.m                   read .mzXML files
  step0_read_hit_list.m          read the hit list, match the raw files
  step1_extract_spectra.m        EIC chromatograms, apex, MS1 and MS2 spectra
  step2_fold_changes.m           fold change against the median of all samples
  step3_precursor_purity.m       purity of precursor ion inside the isolation window of the quadrupole
  step4_mass_accuracy.m          mass deviation of the isolated precursor
  step5_clean_ms2_fragments.m    removing non-fragmented precursor ion and low intensity MS2 signals
  step6_normalise_fragments.m    normalise intensity of MS2 fragments to highest fragment
  step7_match_predicted_spectra.m compare with predicted spectra from CFM-ID 4.0
  step8_export_qc.m              QC summary excel sheet
  step9_export_fragments.m       export excel sheet with detailed results (e.g. fragments, predicted fragments)
  read_predicted_spectrum.m      open fragmentation tree (downloaded as .txt files) from CFM-ID4.0
  isolation_window.m             isolation window of the quadrupole window
  plot_sample.m                  the five panels the app draws

```

## Output

Two excel sheets are written next to the saved `.mat`:

- **`<name>_QC.xlsx`** — one row per sample: the precursor check-ups, retention
  time, purity check-up, fold change, fragment counts and the fraction of
  fragments matched to a prediction.
- **`<name>.xlsx`** — the same, plus the MS2 spectra themselves: measured
  masses, intensities, normalised intensities, matched predicted masses,
  SMILES and carbon counts, one semicolon-separated list per collision energy,
  nCandidates (how many different SMILES are available for this m/z from prediction).

## References

STAR protocol: Add Reference for STAR protocol here.

MSConvert: Chambers, M.C., Maclean, B., Burke, R., Amodei, D., Ruderman, D.L., Neumann, S., Gatto, L., Fischer, B., Pratt, B., Egertson, J., et al. (2012). A cross-platform toolkit for mass spectrometry and proteomics. Nat Biotechnol 30, 918–920. https://doi.org/10.1038/nbt.2377.

CFM-ID 4.0: Wang, F., Liigand, J., Tian, S., Arndt, D., Greiner, R., Wishart, D.S., 2021. CFM-ID 4.0: More Accurate ESI-MS/MS Spectral Prediction and Compound Identification. Analytical Chemistry 93, 11692–11700. https://doi.org/10.1021/acs.analchem.1c01465

Code has been optimized with Claude Code Opus 5.
