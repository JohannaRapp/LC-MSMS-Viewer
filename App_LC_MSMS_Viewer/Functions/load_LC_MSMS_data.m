%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: load_LC_MSMS_data.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Main script for "App_LC_MSMS_Viewer.mlapp".
% Combines different steps in loading and analysing targeted
% LC-MS/MS data. Select raw data for import and select a 
% filename for saving imported and analysed data.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function load_LC_MSMS_data
% select converted targeted LC-MS/MS datafiles
[files, path] = uigetfile('.mzXML', 'Select .mzXML LC-MS/MS datafiles','MultiSelect','on');

% select hitList
[file_hitList, path_hitList] = uigetfile('.xlsx', 'Select hitList .xlsx sheet');

% select filename for saving
[save_file, save_path] = uiputfile('.mat', 'Select filename for saving imported LC-MS/MS data');
Date = string(datetime('today', 'Format', 'yyyyMMdd'));
filename = append(save_path, Date, "_LC_MSMS_data_", save_file);
filename1 = replace(filename, '.mat','');

% step 0 - Add HitList
[Sample_inf_table,EIC_pos, EIC_neg, hitList] = step0_Add_hitList(path_hitList, file_hitList, files);

% step 1: Extract data from mzXML files and create EICs
[LC_MSMS_data] = step1_Extract_mzXML(Sample_inf_table,EIC_pos, EIC_neg, path, files, hitList);

% step 2 - Calculate fold-changes
[LC_MSMS_data] = step2_foldchanges_EIC(LC_MSMS_data);

% step 3 - Quality Control MS1
[LC_MSMS_data] = step3_QualityControl_MS1(LC_MSMS_data);

% step 4 - Check main MS1 peak in isolation window + compare with precMz
[LC_MSMS_data] = step4_CheckMainMS1peak(LC_MSMS_data);

% step 5 - Remove MS2 fragments which have the same mass as/are larger than precursor ion
[LC_MSMS_data] = step5_RemoveMS2fragments_PrecursorIons(LC_MSMS_data);

% step 6 - Add nomrlaized Intensities, delete normInt <0.05
[LC_MSMS_data] = step6_Add_normalizedInt(LC_MSMS_data);

% step 7 - Add prediction 
[LC_MSMS_data] = step7_Add_PredSpectra(LC_MSMS_data);

% step8 - export QC
[LC_MSMS_data] = step8_export_QC_Pred(LC_MSMS_data, filename1);

% step 9 - exportFragments_ScanNo
[LC_MSMS_data] = step9_exportFragments_ScanNo(LC_MSMS_data, filename1);

save(filename, 'LC_MSMS_data', '-v7.3')