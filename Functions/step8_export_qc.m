function outputFile = step8_export_qc(lcData, outputFile)
%STEP8_EXPORT_QC  
%
%   OUTPUTFILE = STEP8_EXPORT_QC(LCDATA, OUTPUTFILE) writes the summary
%   excel eheet. Results of Quality Control (mass accuracy and purity of isolated
%   precursor ion) are exported to an excel-sheet.
%
%   The columns beyond the sample description are
%
%       PrecMz              precursor the instrument isolated
%       PrecMzDelta_mDa     hit list mass minus that (mss deviation)
%       PrecMzCheck         whether the method targeted the right mass
%       CECheck             whether all collision energies were acquired
%       RT                  retention time of the peak, in seconds
%       Detected            whether the metabolite was seen at all
%       IntPrecMz           height of the precursor m/z
%       deviation           distance from the isolated mass to the highest
%                           peak in that window
%       Passed              overall purity check
%       Int, SidePeaks      the two purity tests separately
%       fc                  fold change against the median peak intensity of all samples
%       nFrag_CE10/20/40    number of fragments left after cleaning MS2
%                           spectra
%       pctPred_CE10/20/40  fraction of those matched to a prediction
%       PredSpectraFound    whether a CFM-ID file existed
%
%   See also STEP9_EXPORT_FRAGMENTS, STEP0_READ_HIT_LIST.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    lcData (1,2) struct
    outputFile (1,:) char
end

if isfile(outputFile)
    delete(outputFile);
end

perMode = cell(2, 1);

for mode = 1:2
    nSamples = height(lcData(mode).samples);
    if nSamples == 0
        continue
    end

    params = lcData(mode).params;
    energyFields = compose("CE%d", params.expectedCollisionEnergies);

    nFragments = nan(nSamples, numel(energyFields));
    fractionPredicted = nan(nSamples, numel(energyFields));

    for z = 1:nSamples
        for f = 1:numel(energyFields)
            fragments = lcData(mode).fragments(z).(energyFields{f});
            nFragments(z, f) = size(fragments, 1);

            matched = lcData(mode).prediction(z).(energyFields{f}).predictedMz;
            if nFragments(z, f) > 0
                fractionPredicted(z, f) = sum(~isnan(matched)) / nFragments(z, f);
            end
        end
    end

    summary = lcData(mode).samples;
    summary.PrecMz          = lcData(mode).precursorMz;
    summary.PrecMzDelta_mDa = 1000 * lcData(mode).precursorDeltaDa;
    summary.PrecMzCheck     = yes_no(lcData(mode).precursorOk, "Correct", "WRONG");
    summary.CECheck         = yes_no(lcData(mode).collisionOk, "Correct", "WRONG");
    summary.RT              = [lcData(mode).apex.retentionTime].';
    summary.Detected        = yes_no([lcData(mode).apex.detected].', "Yes", "No");
    summary.IntPrecMz       = numeric_column({lcData(mode).purity.precursorIntensity});
    summary.deviation       = [lcData(mode).purity.massDeviationDa].';
    summary.Passed          = [lcData(mode).purity.verdict].';
    summary.Int             = [lcData(mode).purity.intensityVerdict].';
    summary.SidePeaks       = [lcData(mode).purity.sidePeakVerdict].';
    summary.fc              = lcData(mode).foldChange;

    for f = 1:numel(energyFields)
        summary.("nFrag_" + energyFields{f}) = nFragments(:, f);
        summary.("pctPred_" + energyFields{f}) = fractionPredicted(:, f);
    end

    perMode{mode} = summary;
end

results = vertcat(perMode{:});

% The key and row index are internal bookkeeping, not results.
results = removevars(results, intersect({'key', 'hitRow', 'modeIndex'}, ...
    results.Properties.VariableNames));

writetable(results, outputFile);
fprintf('Exported QC summary to %s\n', outputFile);

end

% ------------------------------------------------------------------------

function text = yes_no(flags, whenTrue, whenFalse)
%YES_NO  Turn a logical column into a readable one.

text = repmat(string(whenFalse), numel(flags), 1);
text(logical(flags)) = whenTrue;

end

% ------------------------------------------------------------------------

function column = numeric_column(values)
%NUMERIC_COLUMN  Pack per-sample scalars, using NaN where one is missing.
%
%   A sample with no peak in the isolation window has an empty precursor
%   intensity.

column = nan(numel(values), 1);
for k = 1:numel(values)
    if ~isempty(values{k})
        column(k) = double(values{k});
    end
end

end
