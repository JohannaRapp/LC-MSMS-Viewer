function outputFile = step9_export_fragments(lcData, outputFile)
%STEP9_EXPORT_FRAGMENTS  The full result, fragment lists included.
%
%   OUTPUTFILE = STEP9_EXPORT_FRAGMENTS(LCDATA, OUTPUTFILE) writes the
%   detailed workbook: everything in the QC summary, plus the MS2 spectra
%   themselves and what each fragment was matched to.
%
%   Per collision energy it writes the scan number used and six
%   semicolon-separated lists, aligned with each other so that the nth entry
%   of every list describes the same fragment:
%
%       mz_CEnn           measured MS2 fragment masses
%       int_CEnn          their intensities
%       normInt_CEnn      their intensities relative to the highest peak
%       mzPred_CEnn       matching predicted mass, NaN where unmatched
%       SMILES_CEnn       structure of that prediction (if there is more than
%                         one matching smiles, only first one is reported; 
%                         number of matching smiles is shown in nCan_CEnn; 
%                         NaN where unmatched
%       noCarbon_CEnn     its carbon count, NaN where unmatched
%       nCand_CEnn        how many predicted structures shared that mass
%
%   See also STEP8_EXPORT_QC.
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
    energies = params.expectedCollisionEnergies;
    energyFields = compose("CE%d", energies);

    detail = lcData(mode).samples;
    detail.PrecMz     = lcData(mode).precursorMz;
    detail.RT         = [lcData(mode).apex.retentionTime].';
    detail.foldChange = lcData(mode).foldChange;
    detail.ScanNoMS1  = [lcData(mode).apex.scanIndex].';
    detail.IntPrecMz       = [lcData(mode).purity.precursorIntensity].'; %new
    detail.deviation       = [lcData(mode).purity.massDeviationDa].';%new
    detail.Passed          = [lcData(mode).purity.verdict].';%new
    detail.IntPassed             = [lcData(mode).purity.intensityVerdict].';%new
    detail.SidePeaksPassed       = [lcData(mode).purity.sidePeakVerdict].';%new

    for f = 1:numel(energyFields)
        field = energyFields{f};

        scanNumber = nan(nSamples, 1);
        mzText      = strings(nSamples, 1);
        intText     = strings(nSamples, 1);
        normText    = strings(nSamples, 1);
        predText    = strings(nSamples, 1);
        smilesText  = strings(nSamples, 1);
        carbonText  = strings(nSamples, 1);
        candText    = strings(nSamples, 1);

        for z = 1:nSamples
            scanNumber(z) = scan_for_energy(lcData(mode), z, energies(f));

            fragments  = lcData(mode).fragments(z).(field);
            normalised = lcData(mode).normalisedFragments(z).(field);
            matches    = lcData(mode).prediction(z).(field);

            if isempty(fragments)
                [mzText(z), intText(z), normText(z), predText(z), ...
                 smilesText(z), carbonText(z), candText(z)] = deal("NaN");
                continue
            end

            mzText(z)     = join_column(fragments(:, 1));
            intText(z)    = join_column(fragments(:, 2));
            normText(z)   = join_column(normalised(:, 2));
            predText(z)   = join_column(matches.predictedMz);
            smilesText(z) = strjoin(matches.smiles, "; ");
            carbonText(z) = join_column(matches.carbonCount);
            candText(z)   = join_column(matches.nCandidates);
        end

        detail.("ScanNo_" + field)  = scanNumber;
        detail.("mz_" + field)      = mzText;
        detail.("int_" + field)     = intText;
        detail.("normInt_" + field) = normText;
        detail.("mzPred_" + field)  = predText;
        detail.("SMILES_" + field)  = smilesText;
        detail.("noCarbon_" + field) = carbonText;
        detail.("nCand_" + field)   = candText;
    end

    perMode{mode} = detail;
end

results = vertcat(perMode{:});
results = removevars(results, intersect({'key', 'hitRow', 'modeIndex'}, ...
    results.Properties.VariableNames));

writetable(results, outputFile);
fprintf('Exported fragment detail to %s\n', outputFile);

end

% ------------------------------------------------------------------------

function scanNumber = scan_for_energy(modeData, z, energy)
%SCAN_FOR_ENERGY  Which scan was acquired at this collision energy.

scanNumber = NaN;
chosen = modeData.ms2Scans{z};

for k = 1:numel(chosen)
    if modeData.collisionEnergyOfScan{z}(k) == energy
        scanNumber = chosen(k);
        return
    end
end

end

% ------------------------------------------------------------------------

function text = join_column(values)
%JOIN_COLUMN  A numeric column as one semicolon-separated string.
stringValues = string(values(:));
stringValues(ismissing(stringValues)) = "NaN";
text = strjoin(stringValues, "; "); 

end
