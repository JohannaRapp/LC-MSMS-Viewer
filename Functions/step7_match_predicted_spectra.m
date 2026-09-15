function lcData = step7_match_predicted_spectra(lcData, predictedFolder)
%STEP7_MATCH_PREDICTED_SPECTRA  Compare measured fragments with CFM-ID.
%
%   LCDATA = STEP7_MATCH_PREDICTED_SPECTRA(LCDATA, PREDICTEDFOLDER) looks up
%   the in-silico fragmentation tree of each metabolite and asks, for every measured
%   fragment, whether the model predicted a structure at that m/z.
%
%   PREDICTEDFOLDER holds fragmentation trees as .txt file from CFM-ID 4.0 
%   named <abbreviation>_<polarity>.txt.
%
%   Isobaric metabolites share a hit list row, their abbreviations joined by
%   PARAMS.isobarSeparator, and their predictions are pooled: "26dapLL-26dapM"
%   reads both 26dapLL_pos.txt and 26dapM_pos.txt. Because the two cannot be
%   told apart by mass, a fragment matching either counts as a match.
%
%   Each polarity gains a prediction entry per sample with, per collision
%   energy, one row per measured fragment:
%
%       predictedMz    the matching predicted mass, NaN if none
%       smiles         its structure, "NaN" if none; if there are more than
%                      one matching smiles, only the first one is shown (see 
%                      nCandidates: shows how many matching smiles are found in total)
%       carbonCount    number of carbons in that structure, NaN if none
%       nCandidates    how many predicted structures share that mass
%
%   The samples table gains predictionFound, saying whether a prediction
%   file existed for each part of the metabolite name.
%
%   See also READ_PREDICTED_SPECTRUM, STEP8_EXPORT_QC.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    lcData (1,2) struct
    predictedFolder (1,:) char
end

if ~isfolder(predictedFolder)
    error('step7:folderNotFound', ...
        'Predicted spectra folder not found: %s', predictedFolder);
end

for mode = 1:2
    nSamples = height(lcData(mode).samples);
    if nSamples == 0
        continue
    end

    params = lcData(mode).params;
    energyFields = compose("CE%d", params.expectedCollisionEnergies);

    predictions = cell(nSamples, 1);
    foundText = strings(nSamples, 1);

    for z = 1:nSamples
        abbreviation = lcData(mode).samples.abbreviation(z);
        polarity = lcData(mode).samples.polarity(z);

        % One prediction file per isobar sharing this row.
        parts = strsplit(abbreviation, params.isobarSeparator);
        wanted = fullfile(predictedFolder, ...
            compose("%s_%s.txt", string(parts(:)), polarity));

        found = isfile(wanted);

        answer = repmat("No", numel(found), 1);
        answer(found) = "Yes";
        foundText(z) = strjoin(answer, "; ");

        perIsobar = cell(numel(wanted), 1);
        for p = 1:numel(wanted)
            perIsobar{p} = read_predicted_spectrum(char(wanted(p)));
        end
        predicted = vertcat(perIsobar{:});

        % Pool the isobars, then drop structures listed more than once.
        if ~isempty(predicted)
            [~, firstOccurrence] = unique(predicted.smiles, 'stable');
            predicted = predicted(firstOccurrence, :);
        end

        entry = struct();
        for f = 1:numel(energyFields)
            field = energyFields{f};
            measured = lcData(mode).fragments(z).(field);

            if isempty(measured)
                entry.(field) = empty_matches();
                continue
            end

            entry.(field) = match_fragments(measured(:, 1), predicted, ...
                params.predictionToleranceDa);
        end

        predictions{z} = entry;
    end

    lcData(mode).prediction = [predictions{:}];
    lcData(mode).samples.predictionFound = foundText;
end

end

% ------------------------------------------------------------------------

function matches = match_fragments(measuredMz, predicted, toleranceDa)
%MATCH_FRAGMENTS  Find a predicted structure for each measured fragment.

nMeasured = numel(measuredMz);

matches.predictedMz  = nan(nMeasured, 1);
matches.smiles       = repmat("NaN", nMeasured, 1);
matches.carbonCount  = nan(nMeasured, 1);
matches.nCandidates  = zeros(nMeasured, 1);

if isempty(predicted)
    return
end

for k = 1:nMeasured
    candidates = find(abs(predicted.mz - measuredMz(k)) <= toleranceDa);

    if isempty(candidates)
        continue
    end

    % Report the closest in mass, and say how many others were as close.
    [~, best] = min(abs(predicted.mz(candidates) - measuredMz(k)));
    chosen = candidates(best);

    matches.predictedMz(k) = predicted.mz(chosen);
    matches.smiles(k)      = predicted.smiles(chosen);
    matches.carbonCount(k) = count(upper(predicted.smiles(chosen)), "C");
    matches.nCandidates(k) = numel(candidates);
end

end

% ------------------------------------------------------------------------

function matches = empty_matches()
%EMPTY_MATCHES  Template used when a spectrum has no fragments left.

matches = struct('predictedMz', [], 'smiles', strings(0, 1), ...
    'carbonCount', [], 'nCandidates', []);

end
