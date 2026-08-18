%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: import_pred_spectra.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pred_spectra, mz_pred, smiles_pred]=import_pred_spectra(Met_kegg_pol)

% import predicted spectrum (pred_spectra=readtable(Met_CE10_pol,
% 'Delimiter', ' ')) is also possible, but in a few cases spectra
% are imported wrong --> therefor code below
pred_spectra = {};
mz_pred = [];
smiles_pred = {}; 
for i = 1:size(Met_kegg_pol,1)
    clearvars -except Met_kegg_pol i pred_spectra mz_pred smiles_pred
    if exist(Met_kegg_pol(i))
        % open text file
        delimiter = {'\n'}; % split after line break
        fileID = fopen(Met_kegg_pol(i));
        data = textscan(fileID, '%q', 'Delimiter',delimiter);
        data_import = data{1,1};
        data = data{1,1};
        fclose(fileID);
        idx = find(cellfun(@isempty, data)==1); % first idx +1 is start of possible fragments
        data(1:idx(1)) = [];
        idx1 = find(cellfun(@isempty, data)==1);
        if ~isempty(idx1)
            data(idx1(1):end,:) = [];
        end
    
        for z = 1:size(data,1)
            data_split(z,:) =strsplit(string(data(z)),' ');
        end
        pred_spectra_temp = array2table(data_split, 'VariableNames',{'No', 'mz', 'SMILES'});
        mz_pred_temp = double(pred_spectra_temp.mz);
        smiles_pred_temp  = cellstr(pred_spectra_temp.SMILES);
        pred_spectra = [pred_spectra; pred_spectra_temp];
        mz_pred = [mz_pred; mz_pred_temp];
        smiles_pred = [smiles_pred; smiles_pred_temp];
    
    else % if there is no predicted spectrum
        pred_spectra_temp = [];
        mz_pred_temp = [];
        smiles_pred_temp  = [];
        pred_spectra = [pred_spectra; pred_spectra_temp];
        mz_pred = [mz_pred; mz_pred_temp];
        smiles_pred = [smiles_pred; smiles_pred_temp];
    end
end
if size(Met_kegg_pol,1)>1
    % find double entries 
    [~, idx_unique] = unique(smiles_pred,'stable');
    mz_pred = mz_pred(idx_unique);
    smiles_pred = smiles_pred(idx_unique);
    pred_spectra = pred_spectra(idx_unique,:);
end

