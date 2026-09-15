function plot_sample(modeData, z, chromatogramAxes, isolationAxes, ms2Axes)
%PLOT_SAMPLE  Draw the five panels describing one sample.
%
%   PLOT_SAMPLE(MODEDATA, Z, CHROMATOGRAMAXES, ISOLATIONAXES, MS2AXES) draws
%   sample Z of one polarity:
%
%     * the chromatogram of its target metabolite in red (red dot shows 
%       timepoint where MS2 spectra are extracted, median intesity of
%       all samples in black dashed line
%     * the MS1 spectrum of MS1 scan with highets intensity of precursor ion
%       in a mass rangfe of +/ 0.65 m/z (isolation wisth of the quadrupole)
%     * the three normalised MS2 spectra, one per collision energy,
%     fragmnets matching predoicted fragments are labelled with a star
%     ("*")
%
%   See also LC_MSMS_Viewer_V2, STEP2_FOLD_CHANGES.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    modeData (1,1) struct
    z (1,1) double
    chromatogramAxes matlab.ui.control.UIAxes
    isolationAxes matlab.ui.control.UIAxes
    ms2Axes matlab.ui.control.UIAxes
end

params = modeData.params;
apex = modeData.apex(z);
column = apex.targetColumn;

%% --- chromatogram --------------------------------------------------
cla(chromatogramAxes);
hold(chromatogramAxes, 'on');

time = modeData.eic(z).retentionTime(:, column);
own  = modeData.eic(z).intensity(:, column);

plot(chromatogramAxes, modeData.eic(z).referenceTime, ...
    modeData.eic(z).medianIntensity, 'k--', 'LineWidth', 1);
plot(chromatogramAxes, time, own, 'r-', 'LineWidth', 1);
plot(chromatogramAxes, time(apex.rowInEic), own(apex.rowInEic), 'r.', 'MarkerSize', 15);
xline(chromatogramAxes, apex.retentionTime, 'r');

xlabel(chromatogramAxes, 'Time (s)');
ylabel(chromatogramAxes, 'Intensity');
title(chromatogramAxes, {char(modeData.samples.gene(z)), ...
    char("EIC: " + modeData.samples.abbreviation(z))}, 'Interpreter', 'none');
ylim(chromatogramAxes, [0 inf]);
hold(chromatogramAxes, 'off');

%% --- MS1 spectrum in m/z range of the isolation window --------------------------------
cla(isolationAxes);
hold(isolationAxes, 'on');

precursorMz = modeData.precursorMz(z);
window = isolation_window(modeData.ms1(z).mz, precursorMz, params.isolationHalfWidthDa);

plot(isolationAxes, double(modeData.ms1(z).mz(window)), ...
    double(modeData.ms1(z).intensity(window)), 'k-');
xline(isolationAxes, precursorMz, 'k');

xlabel(isolationAxes, '{\it m/z}');
ylabel(isolationAxes, 'Intensity');
title(isolationAxes, sprintf('MS1 in isolation window (%s)', ...
    modeData.purity(z).verdict));
xlim(isolationAxes, precursorMz + [-1 1] * params.isolationHalfWidthDa);
hold(isolationAxes, 'off');

%% --- MS2 spectra, one per collision energy --------------------------
energyFields = compose("CE%d", params.expectedCollisionEnergies);

for f = 1:numel(energyFields)
    axesHandle = ms2Axes(f);
    cla(axesHandle);
    hold(axesHandle, 'on');

    fragments = modeData.normalisedFragments(z).(energyFields{f});
    matched = modeData.prediction(z).(energyFields{f});

    for k = 1:size(fragments, 1)
        plot(axesHandle, fragments(k, 1) * [1 1], [0 fragments(k, 2)], ...
            'k-', 'LineWidth', 1.5);

        % A star marks a fragment the prediction accounts for.
        if k <= numel(matched.predictedMz) && ~isnan(matched.predictedMz(k))
            text(axesHandle, fragments(k, 1), fragments(k, 2), '*', ...
                'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'center', ...
                'FontSize', 18);
        end
    end

    xline(axesHandle, precursorMz, 'k--');
    xlabel(axesHandle, '{\it m/z}');
    ylabel(axesHandle, 'Norm. intensity');
    title(axesHandle, energyFields{f});
    xlim(axesHandle, [50 precursorMz + 50]);
    ylim(axesHandle, [0 1.15]);
    hold(axesHandle, 'off');
end

end
