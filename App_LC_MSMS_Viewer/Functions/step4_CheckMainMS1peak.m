%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step4_CheckMainMS1peak.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Calculates mass deviation of isolated precursor ion to
% expected mass. Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step4_CheckMainMS1peak(LC_MSMS_data)

for i = 1:2 % pos/neg loop
    for a = 1:size(LC_MSMS_data(i).SampleInf,1) % sample loop
        prec_mz = LC_MSMS_data(i).PrecMz(a);
        window_plus = prec_mz+0.65; % isolation width=1.3 --> 1.3/2=0.65
        window_minus = prec_mz-0.65; 
        [~, mz_idx_plus] = min(abs(LC_MSMS_data(i).MS1(a).mz - window_plus));
        [~, mz_idx_minus] = min(abs(LC_MSMS_data(i).MS1(a).mz - window_minus));
        [~, max_int_idx] = max(LC_MSMS_data(i).MS1(a).Int(mz_idx_minus:mz_idx_plus));
        mz_intervall = LC_MSMS_data(i).MS1(a).mz(mz_idx_minus:mz_idx_plus);
        mz_max_int = mz_intervall(max_int_idx);
        delta = abs(prec_mz-mz_max_int);
        LC_MSMS_data(i).QC(a).deviation = delta;
        clearvars -except i a LC_MSMS_data
    end % end sample loop
end % end pos/neg loop