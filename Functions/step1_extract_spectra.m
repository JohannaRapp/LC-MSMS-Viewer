function lcData = step1_extract_spectra(samples, eicTargets, rawFolder, params)
%STEP1_EXTRACT_SPECTRA  Read the raw files and pull out chromatograms and spectra.
%
%   LCDATA = STEP1_EXTRACT_SPECTRA(SAMPLES, EICTARGETS, RAWFOLDER, PARAMS)
%   reads every raw file listed in SAMPLES and reduces it to the pieces the
%   later steps work on. SAMPLES and EICTARGETS come from
%   STEP0_READ_HIT_LIST and RAWFOLDER is where the .mzXML files live.
%
%   LCDATA is a 1-by-2 struct array, index 1 positive polarity and index 2
%   negative, each entry carrying
%
%       samples          the rows of SAMPLES measured in this polarity
%       targets          the masses a chromatogram was extracted for
%       precursorMz      precursor the instrument actually isolated
%       precursorDeltaDa hit list mass minus that, per sample
%       precursorOk      whether the method targeted the requested mass
%       collisionOk      whether all expected collision energies were run
%       eic(z)           chromatograms of every target in this polarity:
%                        retentionTime, intensity, scanNumber and measuredMz,
%                        each nMs1Scans-by-nTargets
%       apex(z)          where this sample's own metabolite eluted:
%                        scanIndex, retentionTime, intensity, detected
%       ms1(z)           mz and intensity of the MS1 scan at the apex
%       ms2Scans(z)      indices of the MS2 scans taken at the apex
%       fragments(z)     picked MS2 peaks per collision energy, as CE10,
%                        CE20 and CE40, each an [m/z intensity] matrix
%       ms2Spectra(z)    the same three scans unpicked, for plotting
%
%   HOW THE MS2 SCANS ARE CHOSEN
%
%   The acquistion consists of one MS1 survexy scan follwoed by MS2 scans
%   at different collision energies (e.g. CE=10, CE=20, CE=40). MS1 scan
%   with highest intensity for precursor ion (set in the method) is used.
%   MS2 scans between the MS1 scans with the highest and second-highest
%   intenisty for the precurosr ion are used to extract MS2 spectra.
%
%   See also STEP0_READ_HIT_LIST, STEP2_FOLD_CHANGES, READ_MZXML.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    samples table
    eicTargets (1,2) struct
    rawFolder (1,:) char
    params (1,1) struct = lc_default_params()
end

polarityNames = ["pos", "neg"];
energyFields = compose("CE%d", params.expectedCollisionEnergies);

lcData = struct( ...
    'polarity',         {"positive", "negative"}, ...
    'samples',          {table(), table()}, ...
    'targets',          {eicTargets(1).targets, eicTargets(2).targets}, ...
    'precursorMz',      {[], []}, ...
    'precursorDeltaDa', {[], []}, ...
    'precursorOk',      {[], []}, ...
    'collisionOk',      {[], []}, ...
    'eic',              {struct([]), struct([])}, ...
    'apex',             {struct([]), struct([])}, ...
    'ms1',              {struct([]), struct([])}, ...
    'ms2Scans',         {{}, {}}, ...
    'collisionEnergyOfScan', {{}, {}}, ...
    'fragments',        {struct([]), struct([])}, ...
    'ms2Spectra',       {{}, {}}, ...
    'params',           {params, params}, ...
    'formatVersion',    {params.formatVersion, params.formatVersion});

for mode = 1:2
    inMode = samples.polarity == polarityNames(mode);
    modeSamples = samples(inMode, :);
    targets = eicTargets(mode).targets;

    nSamples = height(modeSamples);

    lcData(mode).samples = modeSamples;

    if nSamples == 0
        continue
    end

    lcData(mode).precursorMz      = nan(nSamples, 1);
    lcData(mode).precursorDeltaDa = nan(nSamples, 1);
    lcData(mode).precursorOk      = false(nSamples, 1);
    lcData(mode).collisionOk      = false(nSamples, 1);
    lcData(mode).ms2Scans             = cell(nSamples, 1);
    lcData(mode).ms2Spectra           = cell(nSamples, 1);
    lcData(mode).collisionEnergyOfScan = cell(nSamples, 1);

    % Each sample's chromatograms, apex and fragments are built as scalar
    % structs and joined at the end.
    eicPerSample       = cell(nSamples, 1);
    apexPerSample      = cell(nSamples, 1);
    ms1PerSample       = cell(nSamples, 1);
    fragmentsPerSample = cell(nSamples, 1);

    for z = 1:nSamples
        fileName = char(modeSamples.fileName(z));
        fprintf('load %s: file %d of %d in %s polarity\n', ...
            fileName, z, nSamples, lcData(mode).polarity);

        raw = read_mzxml(fullfile(rawFolder, fileName));

        msLevel = [raw.scan.msLevel];
        ms1Scans = find(msLevel == 1);
        ms2Scans = find(msLevel == 2);

        if isempty(ms1Scans) || isempty(ms2Scans)
            error('step1:missingMsLevel', ...
                '%s holds %d MS1 and %d MS2 scans; both are required.', ...
                fileName, numel(ms1Scans), numel(ms2Scans));
        end

        retentionTime = retention_times(raw, fileName);

        % --- what the instrument actually isolated ------------------
        precursorMz = raw.scan(ms2Scans(1)).precursorMz.value;
        deltaDa = modeSamples.targetMz(z) - precursorMz;

        lcData(mode).precursorMz(z)      = precursorMz;
        lcData(mode).precursorDeltaDa(z) = deltaDa;
        lcData(mode).precursorOk(z)      = abs(deltaDa) < params.precursorSetupToleranceDa;

        if ~lcData(mode).precursorOk(z)
            warning('step1:precursorMismatch', ...
                ['%s: the method isolated m/z %.4f but the hit list asks ' ...
                 'for %.4f, a difference of %.1f mDa.'], ...
                fileName, precursorMz, modeSamples.targetMz(z), 1000 * deltaDa);
        end

        lcData(mode).collisionOk(z) = isequal( ...
            unique([raw.scan.collisionEnergy]), params.expectedCollisionEnergies);

        % --- chromatograms for every target in this polarity ---------
        eic = extract_chromatograms(raw, ms1Scans, retentionTime, ...
            targets.targetMz, params.eicToleranceDa);
        eic.abbreviation = targets.abbreviation.';
        eicPerSample{z} = eic;

        % --- where this target metabolite of the sample eluted ---------------
        ownTarget = find(targets.abbreviation == modeSamples.abbreviation(z), 1);
        if isempty(ownTarget)
            error('step1:noOwnTarget', ...
                'No chromatogram target for %s, which %s is meant to contain.', ...
                modeSamples.abbreviation(z), fileName);
        end

        apex = find_apex(eic, ownTarget, ms1Scans);
        apexPerSample{z} = apex;

        % --- MS1 spectrum at the apex --------------------------------
        apexPeaks = raw.scan(apex.scanIndex).peaks.mz;
        ms1PerSample{z} = struct("mz", apexPeaks(1:2:end-1), ...
                                 "intensity", apexPeaks(2:2:end));

        % --- the three MS2 scans at the apex --------------------------
        chosen = choose_ms2_scans(apex, msLevel, params.nMs2ScansPerPeak);
        lcData(mode).ms2Scans{z} = chosen;

        % Record the energy each chosen scan was acquired at.
        lcData(mode).collisionEnergyOfScan{z} = [raw.scan(chosen).collisionEnergy];

        [fragments, spectra] = pick_fragments(raw, chosen, energyFields, params);
        fragmentsPerSample{z} = fragments;
        lcData(mode).ms2Spectra{z} = spectra;

        clear raw
    end

    lcData(mode).eic       = [eicPerSample{:}];
    lcData(mode).apex      = [apexPerSample{:}];
    lcData(mode).ms1       = [ms1PerSample{:}];
    lcData(mode).fragments = [fragmentsPerSample{:}];
end

end

% ========================================================================
% Helpers
% ========================================================================

function seconds = retention_times(raw, fileName)
%RETENTION_TIMES  Retention time of every scan, in seconds.
%
%   mzXML stores these as ISO 8601 durations, "PT3.995S".

nScans = numel(raw.scan);
seconds = nan(nScans, 1);

for k = 1:nScans
    text = char(string(raw.scan(k).retentionTime));
    value = regexp(text, '^PT([0-9.]+)S$', 'tokens', 'once');

    if isempty(value)
        error('step1:badRetentionTime', ...
            'Scan %d of %s has retention time "%s", which is not PT<seconds>S.', ...
            k, fileName, text);
    end

    seconds(k) = str2double(value{1});
end

end

% ------------------------------------------------------------------------

function eic = extract_chromatograms(raw, ms1Scans, retentionTime, targetMz, toleranceDa)
%EXTRACT_CHROMATOGRAMS  One chromatogram per target mass, from the MS1 scans.
%
%   For each MS1 scan and each target, the m/z-values within TOLERANCEDA of the
%   target mass are collected. If several matching m/z-values are found, the one with 
%   the highest intensity for the EIC is used. Scans with no m/z-value near 
%   the target contribute a zero.
%
%   The scan is unpacked once and all targets are tested against it.

nMs1 = numel(ms1Scans);
nTargets = numel(targetMz);

eic.retentionTime = zeros(nMs1, nTargets);
eic.intensity     = zeros(nMs1, nTargets);
eic.scanNumber    = zeros(nMs1, nTargets);
eic.measuredMz    = zeros(nMs1, nTargets);

for b = 1:nMs1
    interleaved = raw.scan(ms1Scans(b)).peaks.mz;
    mz  = interleaved(1:2:end-1);
    int = interleaved(2:2:end);

    eic.retentionTime(b, :) = retentionTime(ms1Scans(b));
    eic.scanNumber(b, :) = b;

    for a = 1:nTargets
        near = abs(mz - targetMz(a)) < toleranceDa;
        if ~any(near)
            continue
        end
        eic.intensity(b, a)  = max(int(near));
        eic.measuredMz(b, a) = mean(mz(near));
    end
end

end

% ------------------------------------------------------------------------

function apex = find_apex(eic, targetColumn, ms1Scans)
%FIND_APEX  The MS1 scan where the peak intensity of the target metabolite is highest.

intensity = eic.intensity(:, targetColumn);
[sorted, order] = sort(intensity, 'descend');

apex.targetColumn   = targetColumn;
apex.detected       = sorted(1) > 0;
apex.rowInEic       = order(1);
apex.secondRowInEic = order(min(2, numel(order)));
apex.scanIndex      = ms1Scans(order(1));
apex.secondScanIndex = ms1Scans(order(min(2, numel(order))));
apex.intensity      = sorted(1);
apex.retentionTime  = eic.retentionTime(order(1), targetColumn);
apex.measuredMz     = eic.measuredMz(order(1), targetColumn);

if ~apex.detected

    warning('step1:notDetected', ...
        ['The target was not found in any MS1 scan of this sample. Its ' ...
         'retention time and MS2 spectra are not meaningful.']);
end

end

% ------------------------------------------------------------------------

function chosen = choose_ms2_scans(apex, msLevel, nWanted)
%CHOOSE_MS2_SCANS  The MS2 scans acquired at the top of the peak.
%
%   They lie between the two most intense MS1 scans: after the apex when the
%   second-highest MS1 scan comes later, before it when it comes earlier.

if apex.secondScanIndex > apex.scanIndex
    step = 1;
else
    step = -1;
end

chosen = apex.scanIndex + step * (1:nWanted);

% Check both ends and check the level, so a run that starts or ends mid-cycle
% gives an empty result rather than MS1 scans treated as fragment spectra.
chosen = chosen(chosen >= 1 & chosen <= numel(msLevel));
chosen = chosen(msLevel(chosen) == 2);

chosen = sort(chosen);

end

% ------------------------------------------------------------------------

function [fragments, spectra] = pick_fragments(raw, chosen, energyFields, params)
%PICK_FRAGMENTS  Peak-pick the chosen MS2 scans, one field per collision energy.

fragments = struct();
for f = 1:numel(energyFields)
    fragments.(energyFields{f}) = zeros(0, 2);
end
spectra = cell(numel(energyFields), 1);

for k = 1:numel(chosen)
    interleaved = raw.scan(chosen(k)).peaks.mz;
    % m/z values and intensities alternate; vector starts with m/z
    mz  = double(interleaved(1:2:end-1));
    int = double(interleaved(2:2:end));

    energy = raw.scan(chosen(k)).collisionEnergy;
    slot = find(params.expectedCollisionEnergies == energy, 1);

    if isempty(slot)
        warning('step1:unexpectedCollisionEnergy', ...
            'Scan %d was acquired at %g eV, which is not one of %s; ignoring it.', ...
            chosen(k), energy, mat2str(params.expectedCollisionEnergies));
        continue
    end

    [height_, location] = findpeaks(int, ...
        'MinPeakHeight', params.ms2MinPeakHeight, ...
        'MinPeakProminence', params.ms2MinPeakProminence);

    fragments.(energyFields{slot}) = [mz(location), height_];
    spectra{slot} = [mz, int];
end

end
