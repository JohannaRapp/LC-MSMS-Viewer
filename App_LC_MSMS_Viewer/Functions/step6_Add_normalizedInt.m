%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step6_Add_normalizedInt.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Normalizes peak intensities of all MS2 fragments to fragment
% with highest peak intensity. Delete MS2 fragemnts which showed normalized 
% peak intensity <0.05. Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step6_Add_normalizedInt(LC_MSMS_data)
for i = 1:2 % pos/neg loop
    clearvars -except i LC_MSMS_data
    for k = 1:size(LC_MSMS_data(i).SampleInf,1)
        for CE = 1:3 % loop over all three CEs
            if CE == 1
                Var = "CE10";       
            elseif CE ==2
                Var = "CE20";   
            elseif CE == 3
                Var = "CE40";   
            end
            if ~isempty(LC_MSMS_data(i).fragMS2(k).(Var))
                frag = LC_MSMS_data(i).fragMS2(k).(Var);
                mz = frag(:,1);
                int = frag(:,2);
                if ~isempty(int)
                    max_int = max(int);
                    norm_Int = int/max_int;

                    % find norm Int <0.05 and delete those entries
                    idx_delete = find(norm_Int<0.05);
                    norm_Int(idx_delete) = [];
                    mz(idx_delete) = [];
                    LC_MSMS_data(i).normFrag(k).(Var) = [mz,norm_Int];
                    
                    % delete also fragMS2 with non-normalized intensity
                    % where normInt is <0.05
                    int(idx_delete) = [];
                    LC_MSMS_data(i).fragMS2(k).(Var) = []; % delete old entry
                    LC_MSMS_data(i).fragMS2(k).(Var) = [mz,int];
                else
                    LC_MSMS_data(i).normFrag(k).(Var) = [];
                end           
            else
                LC_MSMS_data(i).normFrag(k).(Var) = [];
            end
            clearvars frag mz int max_int norm_Int idx_delete
        end % end CE loop
        clearvars -except i k LC_MSMS_data CE
    end % end sample loop
    clearvars -except i k LC_MSMS_data 
end % end pos/neg loop
end