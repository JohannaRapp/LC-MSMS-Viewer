%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step0_Add_hitList.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Select hit List. This hit List contains infromation on which
% Precursor ions are used for targeted LC-MS/MS analysis. Example data
% sheet is found in the folder ("hitList_exampleData.xlsx".
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Sample_inf_table,EIC_pos, EIC_neg, hitList] = step0_Add_hitList(path_hitList, file_hitList, files)

hitList = readtable([path_hitList, file_hitList]); 
hitList.MetStrainPol = append(hitList.MetaboliteAbbreviation, "_", hitList.Gene,"_", hitList.Polarity);
hitList.MetMode = append(hitList.MetaboliteAbbreviation, "_", hitList.Polarity);

[EIC_pos,EIC_neg] = extract_EIC(hitList, files);

% Specify how th raw files are named (e.g.
% metabolie_gene_position_polarity.mzXML -->
% 3psme_aroC_P4C6msAV958_neg.mzXML)
Header_SampleInf=["Metabolite" "Gene" "PositionBatch" "Polarity"];

% Extract Sample Information
for c = 1:size(files,2)
    dummy1 = files{c}; % whole File Name
    dummy2 = {dummy1(1:end-6)}; % File Name without .mzXML
    filesInfo = cellstr(strsplit(string(dummy2),'_'));
    Sample_inf(c,:) = filesInfo;
end

Sample_inf_table = cell2table(Sample_inf, 'VariableNames', Header_SampleInf);
[Sample_inf_table] = AddInfo_allEICs(hitList, Sample_inf_table);
end