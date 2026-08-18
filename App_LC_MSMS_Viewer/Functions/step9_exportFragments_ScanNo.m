%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: step9_exportFragments_ScanNo.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Export of fold-changes, MS2 spectra, scan numbers, QC etc.
% in excel-sheet.
% Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"
% Version: 1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [LC_MSMS_data] = step9_exportFragments_ScanNo(LC_MSMS_data, filename1)

% Add Datafile name
for i = 1:2 % pos/neg loop
    LC_MSMS_data(i).SampleInf.DataFile = append(LC_MSMS_data(i).SampleInf.Metabolite,"_",...
        LC_MSMS_data(i).SampleInf.Gene,"_",LC_MSMS_data(i).SampleInf.PositionBatch,"_", ...
        LC_MSMS_data(i).SampleInf.Polarity,".d");
end % end pos/neg loop

% extract info
final1 = [];
for i = 1:2 % pos/neg loop
    mz = table();
    int = table();
    normInt = table();
    mzPred = table();
    noCarbon = table();
    SMILES = table();
    MS2scan = table();
    for a = 1:size(LC_MSMS_data(i).SampleInf,1) % sample loop
        SampleInf(a,:) = LC_MSMS_data(i).SampleInf(a,:);
        PrecMz(a,:) = LC_MSMS_data(i).PrecMz(a);
        RT(a,1) = LC_MSMS_data(i).RT(a);
        foldchange(a,1) = LC_MSMS_data(i).foldChange(a);
        IntPrecMz(a,1) = LC_MSMS_data(i).QC(a).IntPrecMz;
        deviation(a,1) = LC_MSMS_data(i).QC(a).deviation;
        Passed(a,1) = string(LC_MSMS_data(i).QC(a).Passed);
        IntPassed(a,1) = string(LC_MSMS_data(i).QC(a).IntPassed);
        SidePeaksPassed(a,1) = string(LC_MSMS_data(i).QC(a).SidePeaksPassed);
        ScanNoMS1(a,1) = LC_MSMS_data(i).Scan.ScanNoMS1(a);
        dummy_Scan_noMS2 = sortrows(LC_MSMS_data(i).Scan.ScanNoMS2{a},'ascend');
        for CE = 1:3 % CE loop
            if CE == 1
                Var = "CE10";
            elseif CE ==2
                Var = "CE20";
            elseif CE ==  3
                Var = "CE40";
            end
            frag_temp = LC_MSMS_data(i).fragMS2(a).(Var);
            normInt_temp = LC_MSMS_data(i).normFrag(a).(Var);
            mzPred_temp = LC_MSMS_data(i).prediction(a).mz_pred.(Var);
            noCarbon_temp = LC_MSMS_data(i).prediction(a).CarbonSmiles.(Var);
            SMILES_temp = LC_MSMS_data(i).prediction(a).smiles.(Var);
            if ~ isempty(frag_temp)
                mz.(Var)(a,1) = join(string(frag_temp(:,1)),'; ');
                int.(Var)(a,1) = join(string(frag_temp(:,2)),'; ');
                normInt.(Var)(a,1) = join(string(normInt_temp(:,2)),'; ');
                % find NaN in prediction
                idx_NaN = find(isnan(mzPred_temp));
                mzPred_temp = string(mzPred_temp);
                noCarbon_temp = string(noCarbon_temp);
                SMILES_temp = string(SMILES_temp);
                if ~isempty(idx_NaN)
                    mzPred_temp(idx_NaN) = "NaN";
                    noCarbon_temp(idx_NaN) = "NaN";
                    SMILES_temp(idx_NaN) = {'NaN'};
                end
                mzPred.(Var)(a,1) = join(string(mzPred_temp),'; ');
                noCarbon.(Var)(a,1) = join(string(noCarbon_temp),'; ');
                SMILES.(Var)(a,1) = join(SMILES_temp,'; ');
            else
                mz.(Var)(a,1) = "NaN";
                int.(Var)(a,1) = "NaN";
                normInt.(Var)(a,1) ="NaN";
                mzPred.(Var)(a,1) = "NaN";
                noCarbon.(Var)(a,1) = "NaN";
                SMILES.(Var)(a,1) = "NaN";
            end      
            
            MS2scan.(Var)(a,1) = dummy_Scan_noMS2(CE);

            clearvars -except i a SampleInf PrecMz RT foldchange IntPrecMz...
                deviation Passed IntPassed SidePeaksPassed CE LC_MSMS_data...
                mz int normInt mzPred noCarbon SMILES final1 MS2scan ScanNoMS1 dummy_Scan_noMS2 filename1
        end % end CE loop
        clearvars -except i a LC_MSMS_data SampleInf PrecMz RT foldchange IntPrecMz...
                deviation Passed IntPassed SidePeaksPassed mz int normInt...
                mzPred noCarbon SMILES final1 MS2scan ScanNoMS1 filename1
    end % end sample loop
    mz = renamevars(mz,["CE10","CE20","CE40"], ["mz_CE10","mz_CE20","mz_CE40"]);
    int = renamevars(int,["CE10","CE20","CE40"], ["int_CE10","int_CE20","int_CE40"]);
    normInt = renamevars(normInt,["CE10","CE20","CE40"], ["normInt_CE10","normInt_CE20","normInt_CE40"]);
    mzPred = renamevars(mzPred,["CE10","CE20","CE40"], ["mzPred_CE10","mzPred_CE20","mzPred_CE40"]);
    noCarbon = renamevars(noCarbon,["CE10","CE20","CE40"], ["noCarbon_CE10","noCarbon_CE20","noCarbon_CE40"]);
    SMILES = renamevars(SMILES,["CE10","CE20","CE40"], ["SMILES_CE10","SMILES_CE20","SMILES_CE40"]);
    MS2scan = renamevars(MS2scan, ["CE10","CE20","CE40"], ["ScanNo_CE10","ScanNo_CE20","ScanNo_CE40"]);

    final = [SampleInf,array2table([PrecMz, RT, foldchange, IntPrecMz, deviation],...
        'VariableNames', {'PrecMz','RT', 'fold-change', 'Int PrecMz', 'deviation'}),...
        array2table([Passed, IntPassed, SidePeaksPassed, ScanNoMS1],'VariableNames',...
        {'Passed', 'IntPassed','SidePeaksPassed','ScanNoMS1'}),MS2scan, mz, int, normInt, mzPred, noCarbon,...
        SMILES];
    final1 = [final1; final];
    clearvars -except i LC_MSMS_data final1 filename1
end % end pos/neg loop
filenameX = append(filename1,'.xlsx');
writetable(final1,filenameX)