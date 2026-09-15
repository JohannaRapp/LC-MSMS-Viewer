function [lcData, savedFile] = run_lc_msms_pipeline(rawFolder, hitListFile, predictedFolder, outputBase, params)
%RUN_LC_MSMS_PIPELINE  Import and analyse a set of targeted LC-MS/MS runs.
%
%   [LCDATA, SAVEDFILE] = RUN_LC_MSMS_PIPELINE() asks for the raw files (.mzXML), the
%   hit list, the predicted spectra folder and where to save the imported data,
%   then runs %   steps 0 to 9 in order. LCDATA is empty if any dialog is cancelled.
%
%   [...] = RUN_LC_MSMS_PIPELINE(RAWFOLDER, HITLISTFILE, PREDICTEDFOLDER,
%   OUTPUTBASE, PARAMS) runs without dialogs. OUTPUTBASE is the path and
%   name the results are written under, without an extension: the data set
%   is saved as <base>.mat and the two workbooks as <base>_QC.xlsx and
%   <base>.xlsx.
%
%   The steps are
%
  % run_lc_msms_pipeline.m         runs steps 0 to 9
  % read_mzxml.m                   read .mzXML files
  % step0_read_hit_list.m          read the hit list, match the raw files
  % step1_extract_spectra.m        EIC chromatograms, apex, MS1 and MS2 spectra
  % step2_fold_changes.m           fold change against the median of all samples
  % step3_precursor_purity.m       purity of precursor ion inside the isolation window of the quadrupole
  % step4_mass_accuracy.m          mass deviation of the isolated precursor
  % step5_clean_ms2_fragments.m    removing non-fragmented precursor ion and low intensity MS2 signals
  % step6_normalise_fragments.m    normalise intensity of MS2 fragments to highest fragment
  % step7_match_predicted_spectra.m compare with predicted spectra from CFM-ID 4.0
  % step8_export_qc.m              QC summary excel sheet
  % step9_export_fragments.m       export excel sheet with detailed results (e.g. fragments, predicted fragments)
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    rawFolder (1,:) char = ''
    hitListFile (1,:) char = ''
    predictedFolder (1,:) char = ''
    outputBase (1,:) char = ''
    params (1,1) struct = lc_default_params()
end

lcData = [];
savedFile = '';

%% ------------------------------------------------------------------
%  Where everything is
%  ------------------------------------------------------------------

if isempty(rawFolder)
    [chosenFiles, chosenFolder] = uigetfile('*.mzXML', ...
        'Select .mzXML LC-MS/MS data files', 'MultiSelect', 'on');
    if isequal(chosenFiles, 0)
        return
    end
    if ischar(chosenFiles)
        chosenFiles = {chosenFiles};
    end
    rawFolder = chosenFolder;
    rawFileNames = sort(chosenFiles(:)).';
else
    found = dir(fullfile(rawFolder, '*.mzXML'));
    if isempty(found)
        error('runPipeline:noRawFiles', 'No .mzXML files in %s.', rawFolder);
    end
    rawFileNames = {found.name};
end

if isempty(hitListFile)
    [chosenName, chosenFolder] = uigetfile('*.xlsx', 'Select the hit list');
    if isequal(chosenName, 0)
        return
    end
    hitListFile = fullfile(chosenFolder, chosenName);
end

if isempty(predictedFolder)
    predictedFolder = uigetdir(rawFolder, 'Select the folder with predicted spectra from CFM-ID');
    if isequal(predictedFolder, 0)
        return
    end
end

if isempty(outputBase)
    [chosenName, chosenFolder] = uiputfile('*.mat', 'Enter filename for the analysed data set');
    if isequal(chosenName, 0)
        return
    end
    stamp = char(datetime('today', 'Format', 'yyyyMMdd'));
    [~, stem] = fileparts(chosenName);
    outputBase = fullfile(chosenFolder, [stamp '_LC_MSMS_data_' stem]);
end

%% ------------------------------------------------------------------
%  Run the steps
%  ------------------------------------------------------------------

fprintf('step 0: reading the hit list\n');
[~, samples, eicTargets] = step0_read_hit_list(hitListFile, rawFileNames, params);

fprintf('step 1: reading %d raw files\n', numel(rawFileNames));
lcData = step1_extract_spectra(samples, eicTargets, rawFolder, params);

fprintf('step 2: fold changes\n');
lcData = step2_fold_changes(lcData);

fprintf('step 3: precursor purity\n');
lcData = step3_precursor_purity(lcData);

fprintf('step 4: mass accuracy\n');
lcData = step4_mass_accuracy(lcData);

fprintf('step 5: cleaning MS2 spectra\n');
lcData = step5_clean_ms2_fragments(lcData);

fprintf('step 6: normalising fragments\n');
lcData = step6_normalise_fragments(lcData);

fprintf('step 7: matching predicted spectra\n');
lcData = step7_match_predicted_spectra(lcData, predictedFolder);

fprintf('step 8: writing the QC summary\n');
step8_export_qc(lcData, [outputBase '_QC.xlsx']);

fprintf('step 9: writing the fragment detail\n');
step9_export_fragments(lcData, [outputBase '.xlsx']);

savedFile = [outputBase '.mat'];
LC_MSMS_data = lcData;
save(savedFile, 'LC_MSMS_data', '-v7.3');
fprintf('Saved %s\n', savedFile);

end
