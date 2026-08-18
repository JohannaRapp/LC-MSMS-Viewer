%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: extract_EIC.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Extracts masses and polarity for EICs.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [EIC_pos,EIC_neg] = extract_EIC(hitList, files)

for i = 1:size(files,2)
    files_split(i,:) = strsplit(string(files(i)), {'_', '.'});
end
files_inf = array2table(files_split, 'VariableNames', {'abb', 'Strain', 'PositionBatch', 'Polarity', 'format'});
files_inf.MetStrainPol = append(files_inf.abb,'_', files_inf.Strain,'_', files_inf.Polarity);

% extract mass from hitList
for i = 1:size(files_inf,1)
    idx = find(files_inf.MetStrainPol(i)==hitList.MetStrainPol);
    files_inf.mode(i) = hitList.Mode(idx);
    files_inf.mass(i) = hitList.Mass(idx);
    files_inf.metName(i) = hitList.MetaboliteAbbreviation(idx);
end
files_inf.metMode = append(files_inf.abb, '_', files_inf.mode);
files_inf.metStrainMode = append(files_inf.abb, '_', files_inf.Strain, '_',files_inf.mode);

idx_neg = files_inf.Polarity=="neg";
idx_pos = files_inf.Polarity=="pos";
EIC_neg = files_inf(idx_neg,:);
EIC_pos = files_inf(idx_pos,:);

end

