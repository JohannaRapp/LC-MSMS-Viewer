function [hitList, samples, eicTargets] = step0_read_hit_list(hitListFile, rawFileNames, params)
%STEP0_READ_HIT_LIST  Read the hit list and match the raw files against it.
%
%   [HITLIST, SAMPLES, EICTARGETS] = STEP0_READ_HIT_LIST(HITLISTFILE,
%   RAWFILENAMES, PARAMS) reads the hit list.
%
%   HITLISTFILE is the .xlsx sheet listing the metabolites that were
%   targeted; RAWFILENAMES is a cell array of raw file names, without paths.
%
%   The hit list must provide the columns
%
%       Gene                    gene knocked down in the CRISPRi strain
%       Metabolite              full metabolite name (isobars merged)
%       MetaboliteAbbreviation  short name, also used for predicted spectra (isobars merged)
%       Polarity                "pos" or "neg"
%       Mode                    adduct, e.g. [M+H]+ or [M-H]-
%       Mass                    m/z of that adduct, i.e. the precursor mass
%
%   Any other columns, such as MonoMass or the fold changes carried over
%   from the discovery screen, are kept but not used.
%
%   HITLIST is that table with a "key" column added. SAMPLES has one row per
%   raw file, in the order the file names were given, carrying what the file name says (abbreviation, gene, position
%   and polarity) joined to what the hit list says (metabolite, adduct,
%   expected precursor mass). EICTARGETS is a 1-by-2 struct array, index 1
%   positive and index 2 negative, listing the distinct masses to extract a
%   chromatogram for in that polarity.
%
%   Every raw file must match exactly one hit list row, and the matching is
%   on abbreviation, gene and polarity together. 
%
%   See also STEP1_EXTRACT_SPECTRA, LC_DEFAULT_PARAMS.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    hitListFile (1,:) char
    rawFileNames (1,:) cell
    params (1,1) struct = lc_default_params()
end

%% ------------------------------------------------------------------
%  Read and check the hit list
%  ------------------------------------------------------------------

if ~isfile(hitListFile)
    error('step0:hitListNotFound', 'Hit list not found: %s', hitListFile);
end

hitList = readtable(hitListFile, 'VariableNamingRule', 'preserve');

% Give the columns valid MATLAB names without the warning readtable emits
% when it does this itself, so that "Kegg ID" becomes KeggID predictably.
hitList.Properties.VariableNames = ...
    matlab.lang.makeValidName(hitList.Properties.VariableNames);

required = {'Gene', 'Metabolite', 'MetaboliteAbbreviation', 'Polarity', 'Mode', 'Mass'};
missing = required(~ismember(required, hitList.Properties.VariableNames));
if ~isempty(missing)
    error('step0:hitListColumns', ...
        'The hit list is missing the column(s) %s.\nIt has: %s', ...
        strjoin(missing, ', '), strjoin(hitList.Properties.VariableNames, ', '));
end

hitList.Polarity = lower(string(hitList.Polarity));
badPolarity = ~ismember(hitList.Polarity, ["pos", "neg"]);
if any(badPolarity)
    error('step0:hitListPolarity', ...
        'Hit list row(s) %s have a Polarity that is neither "pos" nor "neg".', ...
        mat2str(find(badPolarity).'));
end

% One key for linking a raw file to its hit list row, used everywhere.
hitList.key = hit_key(hitList.MetaboliteAbbreviation, hitList.Gene, hitList.Polarity);

duplicateKeys = duplicates_of(hitList.key);
if ~isempty(duplicateKeys)
    error('step0:hitListDuplicates', ...
        'The hit list has more than one row for %s.', strjoin(duplicateKeys, ', '));
end

%% ------------------------------------------------------------------
%  Read what the file names say, and join to the hit list
%  ------------------------------------------------------------------

nFiles = numel(rawFileNames);
fileName      = strings(nFiles, 1);
abbreviation  = strings(nFiles, 1);
gene          = strings(nFiles, 1);
positionBatch = strings(nFiles, 1);
polarity      = strings(nFiles, 1);

for k = 1:nFiles
    parts = regexp(rawFileNames{k}, params.fileNamePattern, ...
        'names', 'once', 'ignorecase');

    if isempty(parts)
        error('step0:badFileName', ...
            ['File name "%s" cannot be interpreted.\nExpected ' ...
             '<abbreviation>_<gene>_<positionBatch>_<polarity>.mzXML with ' ...
             'polarity "pos" or "neg", for example ' ...
             '"3psme_aroC_P4C6msAV958_neg.mzXML".'], rawFileNames{k});
    end

    fileName(k)      = rawFileNames{k};
    abbreviation(k)  = parts.abbreviation;
    gene(k)          = parts.gene;
    positionBatch(k) = parts.positionBatch;
    polarity(k)      = lower(parts.polarity);
end

samples = table(fileName, abbreviation, gene, positionBatch, polarity);
samples.key = hit_key(abbreviation, gene, polarity);

[isKnown, hitRow] = ismember(samples.key, hitList.key);

if ~all(isKnown)
    error('step0:fileNotInHitList', ...
        ['No hit list row for %s.\nThe file name must give the same ' ...
         'abbreviation, gene and polarity as the hit list.'], ...
        strjoin(cellstr(samples.fileName(~isKnown)), ', '));
end

samples.metabolite = string(hitList.Metabolite(hitRow));
samples.mode       = string(hitList.Mode(hitRow));
samples.targetMz   = hitList.Mass(hitRow);
samples.hitRow     = hitRow;

% Index into the 1-by-2 data structure: 1 is positive, 2 is negative.
samples.modeIndex = 2 - double(samples.polarity == "pos");

%% ------------------------------------------------------------------
%  Distinct chromatogram targets per polarity
%  ------------------------------------------------------------------

% Every sample gets a chromatogram for every metabolite targeted in the same
% polarity, which is what makes the fold change in step 2 possible: the same
% mass is followed across all samples of that polarity.

polarityNames = ["pos", "neg"];
eicTargets = struct('polarity', {"positive", "negative"}, 'targets', {table(), table()});

for mode = 1:2
    inMode = samples.polarity == polarityNames(mode);

    targets = table(samples.abbreviation(inMode), samples.metabolite(inMode), ...
        samples.mode(inMode), samples.targetMz(inMode), ...
        'VariableNames', {'abbreviation', 'metabolite', 'mode', 'targetMz'});

    [~, firstOccurrence] = unique(targets.abbreviation, 'stable');
    eicTargets(mode).targets = targets(firstOccurrence, :);
end

end

% ------------------------------------------------------------------------

function key = hit_key(abbreviation, gene, polarity)
%HIT_KEY  The one identifier linking a raw file to its hit list row.

key = lower(string(abbreviation) + "_" + string(gene) + "_" + string(polarity));

end

% ------------------------------------------------------------------------

function repeated = duplicates_of(values)
%DUPLICATES_OF  Values that occur more than once, each reported once.

[unique_, ~, group] = unique(values);
repeated = cellstr(unique_(accumarray(group, 1) > 1));
repeated = repeated(:).';

end
