%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Script: load_LC_MSMS_data.m
% Author: Johanna Rapp
% Date: 17-Aug-2026
% Description: Belongs to Matlab App "App_LC_MSMS_Viewer.mlapp"

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function MakePlot_Median(LC_MSMS_data, app, idx_selection)

metMode_temp = string(LC_MSMS_data(app.pol).SampleInf.metMode(idx_selection(1)));
idx_EIC = find(metMode_temp == LC_MSMS_data(app.pol).EIC(idx_selection(1)).metMode);
idx_EIC = idx_EIC(1); % some EICs occur more than once; because they hve the same vector take first one
idx_RT_temp = LC_MSMS_data(app.pol).RTidx(idx_selection(1));


%% Plot EIC 
xlabel(app.UIAxes,'Time (sec)')
ylabel(app.UIAxes,'Intensity (-)')
title(app.UIAxes, {string(LC_MSMS_data(app.pol).SampleInf.Gene(idx_selection(1))),...
    string(replace(append('EIC: ', LC_MSMS_data(app.pol).SampleInf.metMode(idx_selection(1))),'_', ' '))})

median_int= LC_MSMS_data(app.pol).EIC(idx_selection(1)).median1;
sd_int = LC_MSMS_data(app.pol).EIC(idx_selection(1)).std1;
RT_sample = LC_MSMS_data(app.pol).EIC(idx_selection(1)).RT_interpol;

% code for plotting median and sd
curve1 = median_int'+sd_int';
curve2 = median_int'-sd_int';
RT_sample = RT_sample';
RT2 = [RT_sample, fliplr(RT_sample)];
inBetween = [curve1, fliplr(curve2)];
hold(app.UIAxes,'on')
% fill(app.UIAxes, RT2, inBetween, [0.5 0.5 0.5],'EdgeColor', 'none', 'FaceAlpha', 0.4);
% 
% plot(app.UIAxes, RT_sample, curve1,'Color',[0.5 0.5 0.5])
% plot(app.UIAxes, RT_sample, curve2,'Color',[0.5 0.5 0.5])

plot(app.UIAxes, RT_sample',median_int, 'Color','black', 'LineStyle','--','LineWidth', 1)
plot(app.UIAxes, LC_MSMS_data(app.pol).EIC(idx_selection(1)).RT(:,idx_EIC),LC_MSMS_data(app.pol).EIC(idx_selection(1)).Int(:,idx_EIC),'Color', 'red','LineWidth', 1)
plot(app.UIAxes, LC_MSMS_data(app.pol).EIC(idx_selection(1)).RT(idx_RT_temp,idx_EIC), LC_MSMS_data(app.pol).EIC(idx_selection(1)).Int(idx_RT_temp,idx_EIC), '.', 'Color', 'r', 'MarkerSize',15)  

xline(app.UIAxes,LC_MSMS_data(app.pol).RT(idx_selection(1)), 'Color', 'red')
ylim(app.UIAxes,[0, inf]);
hold(app.UIAxes,'off')

%% mass accuracy of quadrupole
hold(app.UIAxes_2,'on')
prec_mz = LC_MSMS_data(app.pol).PrecMz(idx_selection(1));
window_plus=prec_mz+0.65; % isolation width=1.3 --> 1.3/2=0.65
window_minus=prec_mz-0.65; 

[~, mz_idx_plus] = min(abs(LC_MSMS_data(app.pol).MS1(idx_selection(1)).mz - window_plus));
[~, mz_idx_minus] = min(abs(LC_MSMS_data(app.pol).MS1(idx_selection(1)).mz - window_minus));
plot(app.UIAxes_2, LC_MSMS_data(app.pol).MS1(idx_selection(1)).mz(mz_idx_minus:mz_idx_plus), LC_MSMS_data(app.pol).MS1(idx_selection(1)).Int(mz_idx_minus:mz_idx_plus),  'Color', [0 0 0])
xline(app.UIAxes_2, prec_mz, 'k')
title(app.UIAxes_2, 'MS1 spectrum within isolation width')
xlabel(app.UIAxes_2,'{\it m/z} (-)')
ylabel(app.UIAxes_2, 'Intensity (-)')
xlim(app.UIAxes_2, [window_minus window_plus])
hold(app.UIAxes_2,'off')

%% MS2 spectra (all thre CEs)
for CE = 1:3 % CE loop
    if CE == 1
        Var = "CE10";
        ax = "UIAxes_3";
    elseif CE ==2
        Var = "CE20";
        ax = "UIAxes_4";
    elseif CE == 3
        Var = "CE40";
        ax = "UIAxes_5";
    end
    xlim(app.(ax), [50 prec_mz+50])
    xline(app.(ax), prec_mz, 'black', 'LineStyle', '--')
    xlabel(app.(ax), '{\it m/z} (-)')
    ylabel(app.(ax), 'Norm. intensity (-)')
    if ~isempty(LC_MSMS_data(app.pol).normFrag(idx_selection(1)).(Var))
        fragMS2 = LC_MSMS_data(app.pol).normFrag(idx_selection(1)).(Var);
        mz_temp = fragMS2(:,1);
        int_temp = fragMS2(:,2);
        Pred_temp = LC_MSMS_data(app.pol).prediction(idx_selection(1)).smiles.(Var);
        for j = 1:size(mz_temp,1)
            hold(app.(ax),'on')
            plot(app.(ax), [mz_temp(j,1) mz_temp(j,1)], [int_temp(j,1) 0], 'Color', 'black', 'LineWidth', 1.5)
            if ~isempty(Pred_temp)
                if string(Pred_temp(j,1))~="NaN"
                    text(app.(ax), mz_temp(j,1), int_temp(j,1), '*','VerticalAlignment','bottom', 'HorizontalAlignment','center', 'FontSize', 18)
                end
            end
        end
    end
    title(app.(ax), Var)
    hold(app.(ax),'off')
end % CE loop
end