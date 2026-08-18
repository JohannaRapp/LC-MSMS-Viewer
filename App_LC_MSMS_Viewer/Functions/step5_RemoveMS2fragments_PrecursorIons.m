%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step5_RemoveMS2fragments_PrecursorIons.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: MS2 spectra cleaning. MS2 fragments matching the precursor
% ion are removed. MS2 fragments larger than precursor ion are removed.
% Deleted fragemnts are saved in LC_MSMS_data.deletedMS2.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step5_RemoveMS2fragments_PrecursorIons(LC_MSMS_data)
for i = 1:2 % pos/neg loop
    for k = 1:size(LC_MSMS_data(i).SampleInf,1) % sample loop
        precMz_temp = LC_MSMS_data(i).PrecMz(k);
        for CE = 1:3 % CE loop
            if CE == 1
                Var = "CE10";       
            elseif CE ==2
                Var = "CE20";   
            elseif CE == 3
                Var = "CE40";   
            end
            frag_temp = LC_MSMS_data(i).fragMS2(k).(Var);
            if iscell(frag_temp) % in certain cases we do not have MS2 spectra --> then frag_temp is cell array not matrix
                frag_temp= zeros(0, 2);
            end
            mz_temp = frag_temp(:,1);
            delta = mz_temp-precMz_temp;
            idx = find(abs(delta)<0.003); % fragments which have similar mass than PrecMz
            idx1 = find(mz_temp>precMz_temp); % fragments which are larger than PrecMz
            idx_delete = unique([idx; idx1]); % unique because it is possible that fragment is similar to precMz and larger than precM
            if ~isempty(idx_delete)
                LC_MSMS_data(i).deletedMS2(k).(Var) = frag_temp(idx_delete,:);
                frag_temp(idx_delete,:) = [];
                LC_MSMS_data(i).fragMS2(k).(Var) = frag_temp;
            end
            clearvars -except i k LC_MSMS_data precMz_temp
        end
    end % end sample loop
end % end pos/neg loop
end