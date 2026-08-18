%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: AddInfo_allEICs.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Sample_inf_table] = AddInfo_allEICs(hitList, Sample_inf_table)
Sample_inf_table.MetStrainPol = append(Sample_inf_table.Metabolite,'_', Sample_inf_table.Gene, '_', Sample_inf_table.Polarity);
for i = 1:size(Sample_inf_table,1)
    dummy = find(string(hitList.MetStrainPol)==string(Sample_inf_table.MetStrainPol(i)));
    Sample_inf_table.mode(i) = hitList.Mode(dummy);

end
Sample_inf_table.metMode = append(Sample_inf_table.Metabolite,'_', Sample_inf_table.mode);
end