%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: export_QC_Pred.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Results of Quality COntrol (mass accuracy and purity of isolated
% precursor ion are exported to an excel-sheet.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step8_export_QC_Pred(LC_MSMS_data, filename1)
Headers1 = LC_MSMS_data(1).SampleInf.Properties.VariableNames;
Headers = [Headers1, {'PrecMz','RT', 'IntPrecMz','deviation','Passed', 'Int', 'SidePeaks','fc',...
    'noFrag CE10','noFrag CE20','noFrag CE40','%Pred CE 10','%Pred CE 20','%Pred CE 40'}];
resultsfinal=array2table(zeros(1,size(Headers,2)), 'VariableNames',Headers);

%%
for i = 1:2 % pos/neg loop
    for k = 1:size(LC_MSMS_data(i).SampleInf,1) % Metabolite loop
        PrecMz = array2table(LC_MSMS_data(i).PrecMz(k), 'VariableNames',{'PrecMz'});
        RT = LC_MSMS_data(i).RT(k);
        fc = LC_MSMS_data(i).foldChange(k);

        % CE loop
        for a = 1:3 % CE loop
            if a == 1 
                Var = "CE10";
            elseif a ==2
                Var = "CE20";
            elseif a ==3
                Var = "CE40";
            end
            noFrag.(Var) = size(LC_MSMS_data(i).fragMS2(k).(Var),1);
            Pred = LC_MSMS_data(i).prediction(k).mz_pred.(Var);
            if ~isempty(Pred)
                percentPred.(Var) = size(find(~isnan(Pred)),1)/noFrag.(Var);
            else
                percentPred.(Var) ='NaN';
            end
        end

        if ~isempty(LC_MSMS_data(i).QC(k).IntPrecMz)
            results1 = [string(RT), string(LC_MSMS_data(i).QC(k).IntPrecMz),...
                string(LC_MSMS_data(i).QC(k).deviation), string(LC_MSMS_data(i).QC(k).Passed),...
                string(LC_MSMS_data(i).QC(k).IntPassed), string(LC_MSMS_data(i).QC(k).SidePeaksPassed),...
                fc, noFrag.CE10, noFrag.CE20, noFrag.CE40,percentPred.CE10,...
                percentPred.CE20, percentPred.CE40];
            results1 = array2table(results1,'VariableNames',{'RT', 'IntPrecMz','deviation','Passed',...
                'Int', 'SidePeaks', 'fc','noFrag CE10','noFrag CE20','noFrag CE40',...
                '%Pred CE 10','%Pred CE 20','%Pred CE 40'});
            results = [LC_MSMS_data(i).SampleInf(k,:), PrecMz];
            resultsfinal = [resultsfinal; [results, results1]];
        else
            results1 = [string(RT), "0", string(LC_MSMS_data(i).QC(k).deviation),...
                 string(LC_MSMS_data(i).QC(k).Passed), string(LC_MSMS_data(i).QC(k).IntPassed),...
                 string(LC_MSMS_data(i).QC(k).SidePeaksPassed),fc,...
                 noFrag.CE10, noFrag.CE20, noFrag.CE40, percentPred.CE10, percentPred.CE20, percentPred.CE40];
            results1 = array2table(results1,'VariableNames',{'RT', 'IntPrecMz','deviation','Passed', 'Int', 'SidePeaks', 'fc',...
                'noFrag CE10','noFrag CE20','noFrag CE40','%Pred CE 10','%Pred CE 20','%Pred CE 40'});
            results = [LC_MSMS_data(i).SampleInf(k,:), PrecMz];
            resultsfinal = [resultsfinal; [results, results1]];
        end
   

    end % end metabolite loop
    clearvars -except LC_MSMS_data resultsfinal file filename1
end % end pos/neg loop
% delete first row (is empty)
resultsfinal(1,:) = [];

% save sheet
filename = append(filename1,'_QC.xlsx');
writetable(resultsfinal, filename)
end