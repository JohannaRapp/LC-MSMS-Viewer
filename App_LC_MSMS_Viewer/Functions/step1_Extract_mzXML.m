%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step1_Extract_mzXML.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: This script loads .mzXML data into a matlab structure using
% "mzXMLreadTS.m". Relevant information from the datafiles are extracted:
% - scan numbers of MS1 and MS2 scans and corresponding timepoints
% including collision energies
% - mass of precursor ion isolated in the quadrupole for targeted LC-MS/MS
% analysis
% - MS1 scans are extracted from the data. Extracted ion chromatograms (EICs)
% are created for all metabolites targeted within one polarity. 
% For this purpose, the m/z vector of each MS1 scan is compared with the 
% exact m/z-value of the precursor ion. m/z-values and corresponding intensities, 
% which differed less than Δ=0.003 Da from the theoretical m/z-value, are 
% used for the EIC. If several matching m/z-values are found, the one with 
% the highest intensity for the EIC is used.
% - The three MS2 scans that are in between the MS1 scans with the 
% highest and second highest intensity for the precursor ion are used to 
% extract MS2 spectra. Pick ion peaks in all three MS2 scans via findpeaks.m  
% - Make new matlab structure "LC_MSMS_data" where results are stored.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step1_Extract_mzXML(Sample_inf_table, EIC_pos, EIC_neg, path, files, hitList)

%% Load mzXML files
for i = 1:2 % pos neg loop
    if i == 1
        idx_pol = find(Sample_inf_table.Polarity=="pos");
        EIC = EIC_pos;
    elseif i == 2
        idx_pol = find(Sample_inf_table.Polarity=="neg");
        EIC = EIC_neg;
    end   
    %%
    for z = 1:size(idx_pol,1) % Sample loop (only selected polarity)
        File = char(append(path, files(1,idx_pol(z)))); %.mzXML datafile
    
        % load mzXML-File
        disp("load " + string(files(1,idx_pol(z))) + ": file "+ z + " from "+ size(idx_pol,1)+ " in polarity "+ i)
        mzXMLStruct = mzxmlreadTS(File); % this function imports .mzXML datafile
        
        % extract MS-level
        msLevel =  [mzXMLStruct.scan(:).msLevel]';
        idx_MS1 = find(msLevel==1);
        idx_MS2 = find(msLevel==2);
        
        % Extract RT
        RT = zeros(size(msLevel,1),1);
        for k = 1:size(msLevel,1)
            RT_raw = mzXMLStruct.scan(k).retentionTime;
            RT_extracted = string(RT_raw(1,3:end-1));
            RT(k,1) = str2double(RT_extracted); % retention times for all scans (MS1+MS2) in seconds
        end
               
        % Precursor mz from LC-MS2 datafile
        precMz = mzXMLStruct.scan(idx_MS2(1)).precursorMz.value; % take precursor mz of first MS2 scan (same precMz in each MS2 scan)
        
        % Check if Precursor mz was set correctly in LC-MS2 method
        met_temp = string(Sample_inf_table.metMode(idx_pol(z)));   
        idx_met_hitList = find(met_temp==hitList.MetMode); % use hitList to check mass
        mz_of_met_temp = hitList.Mass(idx_met_hitList);
        
        if abs(mz_of_met_temp-precMz)>0.01
           Sample_inf_table.PrecMzCheck(idx_pol(z)) = "Wrong"; 
        elseif abs(mz_of_met_temp-precMz)<0.01
            Sample_inf_table.PrecMzCheck(idx_pol(z)) = "Correct";
        end
    
        % check if all three Collision Energy (CE) 10,20,40 are measured
        CE = [mzXMLStruct.scan.collisionEnergy];
        CE_unique = unique(CE);
        if isequal(CE_unique, [10, 20, 40])
            Sample_inf_table.CECheck(idx_pol(z)) = "Correct";
        else
            Sample_inf_table.CECheck(idx_pol(z)) = "Wrong";
        end
        
        LC_MSMS_data(i).SampleInf(z,:) = Sample_inf_table(idx_pol(z),:);
        LC_MSMS_data(i).PrecMz(z,1) = precMz;

        % extract MS1 data an create EICs for all metabolites in the hitList (only metabolites measured in same polarity)
        for a = 1:size(EIC,1) % loop all EICs
            precMz_all_temp = EIC.mass(a);
            for b = 1:size(idx_MS1,1) % loop all MS1 scans from selected datafile
                mz = mzXMLStruct.scan(idx_MS1(b)).peaks.mz(1:2:end-1,:);
                int = mzXMLStruct.scan(idx_MS1(b)).peaks.mz(2:2:end,:);
                delta = mz-precMz_all_temp;
                idx_delta = find(abs(delta)<0.003); % find mz-values in MS1 data, which differ only 0.003 mz from precMz
                if ~isempty(idx_delta) % for the case that there are matches (precMZ and mz)
                    cmz_temp = mean(mz(idx_delta)); % Mean value of the actually measured mz values, which are only 0.003 away from the mz value of the reactant
                    int_temp = sort(int(idx_delta),1, 'descend');
                    int_final = int_temp(1); %take highest intensity
        
                    Max(b,1)=int_final; % intensity of highest peak, where delta<0.003
                    Max(b,2)=cmz_temp; % mz-Value (delta <0.003) with highest int
                    Max(b,3)=b; % shows MS1 scan number (NOT: Index from all scans!)
                    Max(b,4)=RT(idx_MS1(b)); % shows Retention time in seconds
            
                else % for the case that precMz does not match to any mz (MS1)
                    Max(b,1)=0;
                    Max(b,2)=0;
                    Max(b,3)=b;
                    Max(b,4)=RT(idx_MS1(b));
                end
            end % end loop MS1 Scans
            
            LC_MSMS_data(i).EIC(z).RT(:,a) = Max(:,4);
            LC_MSMS_data(i).EIC(z).Int(:,a) = Max(:,1);
            LC_MSMS_data(i).EIC(z).ScanNo(:,a) = Max(:,3);
            LC_MSMS_data(i).EIC(z).cmz(:,a) = Max(:,2);
            LC_MSMS_data(i).EIC(z).MetAbb(a) = EIC.abb(a);
            LC_MSMS_data(i).EIC(z).Met(a) = EIC.metName(a);
            LC_MSMS_data(i).EIC(z).mode(a) = EIC.mode(a);
            LC_MSMS_data(i).EIC(z).metMode(a) = EIC.metMode(a);

        end % loop all EICs
        clearvars -except files path i z Sample_inf_table ...
            hitList idx_pol LC_MSMS_data idx_MS1 idx_MS2 mzXMLStruct...
            EIC EIC_pos EIC_neg
        %% extract MS2 scan for where PrecMz has highest intensity
        idx_PrecMz = find(LC_MSMS_data(i).EIC(z).metMode==string(LC_MSMS_data(i).SampleInf.metMode(z)));
        Max = [LC_MSMS_data(i).EIC(z).Int(:,idx_PrecMz(1)),LC_MSMS_data(i).EIC(z).cmz(:,idx_PrecMz(1)),...
            LC_MSMS_data(i).EIC(z).ScanNo(:,idx_PrecMz(1)), LC_MSMS_data(i).EIC(z).RT(:,idx_PrecMz(1))];
    
        % Sort Max regarding intensity --> MS1 scan with highest intensity of precursor mz is in first row
        Max = sortrows(Max,1, 'descend'); % 1.col: Intensity of peak, 2.col: Mass from peak, 3.col: ScanNo in MS1, 4.col:RT in sec
        idxAll_MS1 = idx_MS1(Max(:,3)); % Scan No of highest Intensities (of prec ion) in all scans
        LC_MSMS_data(i).Scan.ScanNoMS1(z) = idxAll_MS1(1);
        % Extract MS1 spectrum for Scan where highest intensity of PrecIon is observed
        mz_MS1 = mzXMLStruct.scan(idxAll_MS1(1)).peaks.mz(1:2:end-1,:);
        int_MS1 = mzXMLStruct.scan(idxAll_MS1(1)).peaks.mz(2:2:end,:);

        % Is second highest MS1 Scan before OR after highest MS1 Scan --> Take the three MS2 Scans which are lying between the 2 highest MS1 scans
        if idxAll_MS1(1) == 1 % When highest scan is no 1 (rarely) --> Problem: measurement starts with 2 MS1 scnas!
            idx_MS2_following = [idxAll_MS1(1)+2;idxAll_MS1(1)+3; idxAll_MS1(1)+4];
        elseif idxAll_MS1(1)<idxAll_MS1(2) 
           idx_MS2_following=[idxAll_MS1(1)+1;idxAll_MS1(1)+2; idxAll_MS1(1)+3];
        else
           idx_MS2_following=[idxAll_MS1(1)-1;idxAll_MS1(1)-2; idxAll_MS1(1)-3];
        end
        LC_MSMS_data(i).Scan.ScanNoMS2(z) = {idx_MS2_following};
        % when highest/2nd highest scan is 1 or 2 --> then there is no peak for this metabolite
        if ~any(idx_MS2_following <1)
           % extract CEs from these three scans (should be 10, 20 and 40)
            for f = 1:size(idx_MS2_following,1)
                CE_temp(f) = mzXMLStruct.scan(idx_MS2_following(f)).collisionEnergy;
            end
            idx_CE10 = CE_temp==10;
            idx_CE20 = CE_temp==20;
            idx_CE40 = CE_temp==40;
    
            % Find peaks in the MS2 scans lying between hoghest MS1 scans for Precursor ion
            HeightFilter = 400; 
            PromFilter = 400; 
    
            % Peak Picking in the highest 3 Scans
            for l = 1:3 
                mz  = mzXMLStruct.scan(idx_MS2_following(l)).peaks.mz(1:2:end-1,:);
                int = mzXMLStruct.scan(idx_MS2_following(l)).peaks.mz(2:2:end,:);
                % peakPicking
                [pks,locs] = findpeaks(int,'MinPeakHeight',HeightFilter,'MinPeakProminence',PromFilter);
                pmz = mz(locs);
                peaks(l,1) = {[pmz,pks]};
                MS2(l,1) = {[mz,int]};
            end
            LC_MSMS_data(i).fragMS2(z).CE10 = peaks{idx_CE10};
            LC_MSMS_data(i).fragMS2(z).CE20 = peaks{idx_CE20};
            LC_MSMS_data(i).fragMS2(z).CE40 = peaks{idx_CE40};
            LC_MSMS_data(i).MS2spec(z).MS2 = MS2;
        else
            LC_MSMS_data(i).fragMS2(z).CE10 = {};
            LC_MSMS_data(i).fragMS2(z).CE20 = {};
            LC_MSMS_data(i).fragMS2(z).CE40 = {};
            LC_MSMS_data(i).MS2spec(z).MS2 = {};
        end
            LC_MSMS_data(i).MS1(z).mz = mz_MS1; % mz-values of MS1 Scan, where precIon has highest intensity
            LC_MSMS_data(i).MS1(z).Int = int_MS1; % intensities of MS1 Scan, where precIon has highest intensity
            LC_MSMS_data(i).RT(z,1) = Max(1,4); % Retention time of the compound (RT of MS1 Scan with highest intensity for precIon)

        clearvars -except files path i z Sample_inf_table ...
            hitList idx_pol LC_MSMS_data EIC EIC_pos EIC_neg

    end % END Sample loop
        clearvars -except files path i z Sample_inf_table ...
            hitList LC_MSMS_data EIC EIC_pos EIC_neg 

end % END pos/neg loop
end