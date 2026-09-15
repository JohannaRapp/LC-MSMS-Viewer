function lcData = step3_precursor_purity(lcData)
%STEP3_PRECURSOR_PURITY  Judge whether the MS2 spectra describe one compound.
%
%   LCDATA = STEP3_PRECURSOR_PURITY(LCDATA) analyzes the MS1 scan with the highest
%   peak intensity of the precursor ion. MS1 spectrum is analyzed within
%   the isolation width of the quadrupole (3 mDa --> 0.65 Da+/- of exact mass of
%   the precuror ion). Checks if there are another peaks besides the wanted
%   m/z peak within this isolation window. 

%   Thresholds are obtained from lc_default_params.
%   params.isolationHalfWidthDa
%       params.purityMinPeakHeight: Minimum peak height when finding peaks
%       within the isolation width of the quadrupole.
%       params.purityMinPrecursorIntensity: Minimum peak height of MS1
%       peak.
%       params.purityRelativeThreshold: A co-isolated peak higher than this 
%       fraction of the precursor peak fails the purity check: the MS2 
%       spectrum is then a mixture.
%   Each polarity gains a purity entry per sample with
%
%       peakMz, peakIntensity   peaks found inside the isolation window
%       relativeIntensity       those peaks relative to the precursor peak
%       precursorIntensity      height of the precursor peak itself
%       verdict                 "Yes", "SidePeaks", "Too Low" or "No peak found"
%       intensityVerdict        "Yes" or "Too Low", the intensity test alone
%       sidePeakVerdict         "Yes" or "SidePeaks", the purity test alone
%
%   The three check-ups are kept separately because a spectrum can fail both
%   tests, and the exports report them in separate columns.
%
%   See also STEP4_MASS_ACCURACY, LC_DEFAULT_PARAMS.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    lcData (1,2) struct
end

for mode = 1:2
    nSamples = height(lcData(mode).samples);
    if nSamples == 0
        continue
    end

    params = lcData(mode).params;
    purity = cell(nSamples, 1);

    for z = 1:nSamples
        precursorMz = lcData(mode).precursorMz(z);
        window = isolation_window(lcData(mode).ms1(z).mz, precursorMz, ...
            params.isolationHalfWidthDa);

        windowMz  = double(lcData(mode).ms1(z).mz(window));
        windowInt = double(lcData(mode).ms1(z).intensity(window));

        [peakIntensity, location] = findpeaks(windowInt, ...
            'MinPeakHeight', params.purityMinPeakHeight);
        peakMz = windowMz(location);

        entry = empty_purity();
        entry.peakMz = peakMz;
        entry.peakIntensity = peakIntensity;

        if isempty(peakMz)
            % Nothing above the noise floor: there is no precursor to judge.
            entry.verdict = "No peak found";
            entry.intensityVerdict = "Too Low";
            entry.sidePeakVerdict = "Yes";
            purity{z} = entry;
            continue
        end

        % The precursor is the peak nearest the isolated mass.
        [~, precursorPeak] = min(abs(peakMz - precursorMz));
        entry.precursorIntensity = peakIntensity(precursorPeak);
        entry.relativeIntensity = peakIntensity / peakIntensity(precursorPeak);

        tooFaint = entry.precursorIntensity < params.purityMinPrecursorIntensity;

        % Count the peaks reaching the threshold: the precursor is always one
        % of them, so more than one means something else was co-isolated.
        hasSidePeaks = sum(entry.relativeIntensity > params.purityRelativeThreshold) > 1;

        if tooFaint
            entry.verdict = "Too Low";
        elseif hasSidePeaks
            entry.verdict = "SidePeaks";
        else
            entry.verdict = "Yes";
        end

        entry.intensityVerdict = ternary(tooFaint, "Too Low", "Yes");
        entry.sidePeakVerdict = ternary(hasSidePeaks, "SidePeaks", "Yes");

        purity{z} = entry;
    end

    lcData(mode).purity = [purity{:}];
end

end

% ------------------------------------------------------------------------

function entry = empty_purity()
%EMPTY_PURITY  Template so every sample has the same fields.

entry = struct('peakMz', [], 'peakIntensity', [], 'relativeIntensity', [], ...
    'precursorIntensity', [], 'verdict', "", 'intensityVerdict', "", ...
    'sidePeakVerdict', "", 'massDeviationDa', NaN);

end

% ------------------------------------------------------------------------

function value = ternary(condition, whenTrue, whenFalse)

if condition
    value = whenTrue;
else
    value = whenFalse;
end

end
