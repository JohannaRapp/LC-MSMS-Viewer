%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step3_QualityControl_MS1.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: MS1 scan with highest peak intensity of precursor ion is
% used to analyze the precursor ion purity. MS1 spectrum is analyzed within
% the isolation width of the quadrupole (3 mDa --> 0.65 Da+/- exact mass of
% the precuror ion). Check if highest peak matches precursor ion (tolerance
% 3 mDa) and check if there are other peaks with more than 20% peak height 
% of the main peak. If there are no peaks bigger than 20% of the main peak 
% the QC is passes. QC is not passed if main peak height is lower than
% 8000.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step3_QualityControl_MS1(LC_MSMS_data)

for k = 1:2 % pos neg loop
    clearvars -except LC_MSMS_data k
    rel_threshold = 0.2; % for normalized peak intensities (on precMz peak)
    height_threshold = 8000; % for absolute intensities 
    
    for l = 1:size(LC_MSMS_data(k).SampleInf,1) % sample loop
        prec_mz = LC_MSMS_data(k).PrecMz(l);
        window_plus = prec_mz+0.65; % isolation width of quadrupole is 1.3 Da --> 1.3/2=0.65
        window_minus = prec_mz-0.65; 
        [~, mz_idx_plus] = min(abs(LC_MSMS_data(k).MS1(l).mz - window_plus));
        [~, mz_idx_minus] = min(abs(LC_MSMS_data(k).MS1(l).mz - window_minus));
        
        % findpeaks in isolation width % and normalization on highest peak
        [pks, locs] = findpeaks(LC_MSMS_data(k).MS1(l).Int(mz_idx_minus:mz_idx_plus), 'MinPeakHeight', 500);
        mz_isoWidth = LC_MSMS_data(k).MS1(l).mz(mz_idx_minus:mz_idx_plus);
        mzpeaks = mz_isoWidth(locs);
        [~, idx] = min(abs(mzpeaks-prec_mz)); % idx precMz peak
        normalized = pks./pks(idx);
        LC_MSMS_data(k).QC(l).PeaksInt = pks;
        LC_MSMS_data(k).QC(l).NormalPeaksInt = normalized;
        LC_MSMS_data(k).QC(l).PeaksMz = mzpeaks;
        LC_MSMS_data(k).QC(l).IntPrecMz = pks(idx);
        if pks(idx)<height_threshold
            LC_MSMS_data(k).QC(l).Passed = 'Too Low';
        elseif  isempty(idx) % when no peask are found
            LC_MSMS_data(k).QC(l).Passed = 'No peak found';
        elseif size(find(normalized>rel_threshold),1)>1 % if there is more than one peak which is gretaer than the threshold (see above) (20% of precMz peak)
            LC_MSMS_data(k).QC(l).Passed = 'SidePeaks';
        else
            LC_MSMS_data(k).QC(l).Passed = 'Yes';
        end
        
        if pks(idx)<height_threshold
            LC_MSMS_data(k).QC(l).IntPassed = 'Too Low';
        elseif  isempty(idx) % when no peask are found
            LC_MSMS_data(k).QC(l).IntPassed = 'Too Low';
        else
            LC_MSMS_data(k).QC(l).IntPassed = 'Yes';
        end

        if size(find(normalized>rel_threshold),1)>1
            LC_MSMS_data(k).QC(l).SidePeaksPassed = 'SidePeaks';
        else
            LC_MSMS_data(k).QC(l).SidePeaksPassed = 'Yes';
        end
        clearvars -except k l LC_MSMS_data rel_threshold height_threshold
    end % end sample loop
end % end pos/neg loop
