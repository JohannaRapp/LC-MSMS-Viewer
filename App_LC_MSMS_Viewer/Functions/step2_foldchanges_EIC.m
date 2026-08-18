%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step2_foldchanges_EIC.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Calculate fold-changes of peak heights from all EICs. 
% Start with accumulating metabolite in a CRISPRi strain, which was
% analyzed via LC-MS/MS. Check in EICs for this metabolite in all other
% datafiles. As retention time vectors might differ a bit between
% datafiles, interpolate rt-vector of all other datafiles to the one with
% the accumulating metabolite. Extract peak heights and calculate fold-changes 
% compared to median on a per batch basis.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step2_foldchanges_EIC(LC_MSMS_data)

for i = 1:2 % pos/neg loop
    clearvars -except LC_MSMS_data i
    
    % RT vectors can have different size --> calculate minimum size
    for x = 1:size(LC_MSMS_data(i).EIC,2)    
        length_RT_sample(x) = size(LC_MSMS_data(i).EIC(x).RT,1);
    end
    RT_length_min = min(length_RT_sample);

    for a = 1:size(LC_MSMS_data(i).SampleInf,1) % sample loop
        met_temp = LC_MSMS_data(i).SampleInf.Metabolite(a);
        idx_EIC = find(LC_MSMS_data(i).EIC(a).MetAbb==string(met_temp));
        idx_EIC = idx_EIC(1); % some EICs can occur more than once; because they have the same vector take first one
        RT_temp_sample = LC_MSMS_data(i).RT(a);
        idx_RT_temp = find(LC_MSMS_data(i).EIC(a).RT(:,idx_EIC)==RT_temp_sample);
        LC_MSMS_data(i).RTidx(a,1) = idx_RT_temp;

        RT_sample = LC_MSMS_data(i).EIC(a).RT(1:RT_length_min,idx_EIC); 
        Int_sample = LC_MSMS_data(i).EIC(a).Int(1:RT_length_min,idx_EIC);
        
        %% loop over all samples for this EIC
        for b = 1:size(LC_MSMS_data(i).EIC,2)
            RT_temp = LC_MSMS_data(i).EIC(b).RT(1:RT_length_min,idx_EIC); % RT vectors have different length --> use length of RT vector of current sample
            Int_temp = LC_MSMS_data(i).EIC(b).Int(1:RT_length_min,idx_EIC);
            Int_temp_interpolate = interp1(RT_temp, Int_temp, RT_sample, 'linear');
            Int_interpolate(:,b) = Int_temp_interpolate; 
        end

        median1 = median(Int_interpolate,2,'omitmissing');
        std1 = std(Int_interpolate,[],2,"omitmissing");
        RTidx = LC_MSMS_data(i).RTidx(a); % RTidx where peak is picked in current sample
        fc = Int_sample(RTidx)./median1(RTidx);
        LC_MSMS_data(i).foldChange(a,1) = fc;
        LC_MSMS_data(i).EIC(a).Int_interpol = Int_interpolate;
        LC_MSMS_data(i).EIC(a).RT_interpol = RT_sample;
        LC_MSMS_data(i).EIC(a).median1 = median1;
        LC_MSMS_data(i).EIC(a).std1 = std1;
        clearvars -except LC_MSMS_data i RT_length_min a
    end
end
end