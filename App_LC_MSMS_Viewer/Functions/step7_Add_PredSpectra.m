%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step7_Add_PredSpectra.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Imports predicted spectra for all metabolites defined in the
% hitList. .txt files with predicted spectra are downloaded from CFM-ID
% (example files are shown in the folder "Predicted spectra"). Predicted 
% Spectra of isobaric metabolites are merged.Measured MS2 fragments are 
% compared with predicted spectra with a mass tolerance of 3 mDa.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step7_Add_PredSpectra(LC_MSMS_data)
for i = 1:2 % pos/neg loop
    for k = 1:size(LC_MSMS_data(i).SampleInf,1)
        clearvars -except LC_MSMS_data i k

        % check if there are isobaric compounds
        isobar_met = strsplit(string(LC_MSMS_data(i).SampleInf.Metabolite(k)),{'-'});
        if size(isobar_met,2)>1
           for s = 1:size(isobar_met,2)
                Met_pol(s,1) = append(isobar_met(s),'_',...
                    LC_MSMS_data(i).SampleInf.Polarity(k),'.txt');
           end
        else
            Met_pol = string(append(LC_MSMS_data(i).SampleInf.Metabolite(k),...
                '_', LC_MSMS_data(i).SampleInf.Polarity(k),'.txt'));
        end
        
        % import predicted spectrum
        [~, mz_pred, smiles_pred] = import_pred_spectra(Met_pol);

        % check if predicted spectra was found
        for b = 1:size(Met_pol,1)
            if exist(Met_pol(b))
                PredSpectraFound(b,1) ="Yes";
            else
                PredSpectraFound(b,1) ="No";
            end
        end
        LC_MSMS_data(i).SampleInf.Found(k) = join(PredSpectraFound,'; ');
        clearvars PredSpectraFound
        
        for CE = 1:3 % loop over all three CEs
            if CE == 1
                Var = "CE10";       
            elseif CE ==2
                Var = "CE20";   
            elseif CE == 3
                Var = "CE40";   
            end
            frag=LC_MSMS_data(i).fragMS2(k).(Var); % fragments of current sample
            if ~isempty(frag)
                mz = frag(:,1);
            else
                mz = [];
            end
        
            if ~isempty(mz)&& ~isempty(mz_pred) % Check if there are matches between measured mz and predicted mz
                for j=1:size(mz,1)
                    mz_temp=mz(j,1);
                    [match, id]=ismembertol(mz_temp, mz_pred, 0.003, 'DataScale', 1);
                                          
                    if match~=0 % for the case that there is a match between 12C frag and mz_pred
                       pred(j,1) = mz_pred(id);
                       smiles(j,1) = smiles_pred(id);
                       CarbonSmiles(j,1) = count(smiles_pred(id), 'C', 'IgnoreCase',true);
                       
                    else % for the case that there is no matching mz-value in the prediction
                        pred(j,1) = NaN;
                        smiles(j,1) = cellstr('NaN');
                        CarbonSmiles(j,1) = NaN;
                    end

                    clearvars mz_temp match id
                end
                    LC_MSMS_data(i).prediction(k).mz_pred.(Var) = pred;
                    LC_MSMS_data(i).prediction(k).smiles.(Var) = smiles;
                    LC_MSMS_data(i).prediction(k).CarbonSmiles.(Var) = CarbonSmiles;
                    clearvars pred smiles CarbonSmiles
            elseif isempty(mz) % when there is no measured fragment
                LC_MSMS_data(i).prediction(k).mz_pred.(Var) = [];
                LC_MSMS_data(i).prediction(k).smiles.(Var) = {};
                LC_MSMS_data(i).prediction(k).CarbonSmiles.(Var) = [];
            elseif isempty(mz_pred) % when there are no fragments
                LC_MSMS_data(i).prediction(k).mz_pred.(Var) = NaN(size(mz,1),1);
                LC_MSMS_data(i).prediction(k).smiles.(Var) = repmat(cellstr('NaN'),size(mz,1),1);
                LC_MSMS_data(i).prediction(k).CarbonSmiles.(Var) = NaN(size(mz,1),1);
            end
                clearvars mz frag
        end % end CE loop
        clearvars mz_pred smiles_pred mz
    end % end sample loop
end % end pos/neg loop
end
