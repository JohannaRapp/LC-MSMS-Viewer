function params = lc_default_params()
%LC_DEFAULT_PARAMS  Central configuration for the LC-MS/MS pipeline.
%
%   PARAMS = LC_DEFAULT_PARAMS() returns a struct holding every constant 
%   the pipeline uses. 
%
%   See also STEP0_READ_HIT_LIST, STEP1_EXTRACT_SPECTRA, READ_MZXML.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp.

% ------------------------------------------------------------------------
% Raw file naming
% ------------------------------------------------------------------------

% Raw files are named <metaboliteAbbreviation>_<gene>_<positionBatch>_<polarity>.mzXML,
% for example "3psme_aroC_P4C6msAV958_neg.mzXML". The first three fields
% identify the sample; the abbreviation and gene together with the polarity
% are what link a file to its row in the hit list.
params.fileNamePattern = ['^(?<abbreviation>.+?)_(?<gene>[^_]+)_' ...
                          '(?<positionBatch>[^_]+)_(?<polarity>pos|neg)\.mzXML$'];

% ------------------------------------------------------------------------
% Extracted ion chromatograms (see STEP1_EXTRACT_SPECTRA)
% ------------------------------------------------------------------------

% Half-width of the window used to pull a target mass out of an MS1 scan.
% A point counts towards the EIC if it lies within this of the expected m/z.
params.eicToleranceDa = 0.003;

% Tolerance for the check that the precursor the instrument actually
% isolated matches the mass requested in the hit list. Deliberately looser
% than eicToleranceDa: this is catching a mis-entered method, not measuring
% mass accuracy, and it should not fire on ordinary calibration drift.
params.precursorSetupToleranceDa = 0.01;

% Collision energies the method is expected to acquire, in eV. Each targeted
% precursor is fragmented once at each of these.
params.expectedCollisionEnergies = [10 20 40];

% How many MS2 scans to take from around the chromatographic peak. The
% acquisition cycles one MS1 scan followed by one MS2 scan per collision
% energy, so the scans bracketed by the two most intense MS1 scans are the
% ones recorded at the top of the peak.
params.nMs2ScansPerPeak = 3;

% ------------------------------------------------------------------------
% MS2 fragment peak picking (see STEP1_EXTRACT_SPECTRA)
% ------------------------------------------------------------------------
% parameters for peak-picking via findpeaks.m function
params.ms2MinPeakHeight = 400;
params.ms2MinPeakProminence = 400;

% ------------------------------------------------------------------------
% Precursor purity and mass accuracy (see STEP3, STEP4)
% ------------------------------------------------------------------------

% Half-width of the quadrupole isolation window, in Da. The instrument
% isolates 1.3 Da around the target mass, so anything within +/- 0.65 Da of the
% precursor is co-isolated and fragmented along with it.
params.isolationHalfWidthDa = 0.65;

% Minimum height for a peak (using findpeaks.m function) to be considered 
% when judging purity inside the isolation window.
params.purityMinPeakHeight = 500;

% A co-isolated peak (within the isolation width of the quadrupole) 
% taller than this fraction of the precursor peak fails the purity check: 
% the MS2 spectrum is then a mixture.
params.purityRelativeThreshold = 0.2;

% Below this absolute intensity the precursor is too weak to trust the MS2
% spectrum, whatever its purity.
params.purityMinPrecursorIntensity = 8000;

% ------------------------------------------------------------------------
% MS2 spectrum cleaning (see STEP5, STEP6)
% ------------------------------------------------------------------------

% Fragments within this of the precursor mass are unfragmented precursor,
% not fragments, and are removed along with all larger fragments.
params.precursorMatchToleranceDa = 0.003;

% Fragments with normalized intensities below this fraction of the normalized 
% peak intensity of the main peak are dropped as noise.
params.fragmentMinRelativeIntensity = 0.05;

% ------------------------------------------------------------------------
% Comparison with predicted spectra (see STEP7)
% ------------------------------------------------------------------------

% Tolerance for calling a measured fragment a match to a CFM-ID prediction.
params.predictionToleranceDa = 0.003;

% Isobaric metabolites that flow injection cannot separate share one hit
% list row, with their abbreviations joined by this character: the row
% "26dapLL-26dapM" needs both 26dapLL_pos.txt and 26dapM_pos.txt.
%
% Note that this means an abbreviation may not itself contain the separator.
params.isobarSeparator = '-';

% Folder holding the CFM-ID .txt predictions, named <abbreviation>_<polarity>.txt.
% Empty means "ask", which is what the app does. Setting it makes a run
% reproducible without any dialog.
params.predictedSpectraFolder = '';

% ------------------------------------------------------------------------
% Provenance
% ------------------------------------------------------------------------

% Written into every saved data set and every exported workbook.
params.softwareVersion = "LC_MSMS_Viewer_V2.0";

% Version of the LC_MSMS_data struct layout itself.
params.formatVersion = 2;

end
