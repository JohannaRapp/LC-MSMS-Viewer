function fragments = read_predicted_spectrum(fileName)
%READ_PREDICTED_SPECTRUM  Read the fragment table from a CFM-ID .txt file.
%
%   FRAGMENTS = READ_PREDICTED_SPECTRUM(FILENAME) returns a table with one
%   row per predicted fragment structure:
%
%       index    the number CFM-ID gives the fragment
%       mz       its mass over charge ratio (m/z)
%       smiles   its structure
%
%   A CFM-ID file holds a header, then one block of predicted peaks per
%   collision energy, then a blank line, then the fragment table this
%   function reads. That last table is the useful one: it lists every
%   fragment structure the model produced together with its mass, pooled
%   across the energies.
%
%   See also STEP7_MATCH_PREDICTED_SPECTRA.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    fileName (1,:) char
end

fragments = table('Size', [0 3], ...
    'VariableTypes', {'double', 'double', 'string'}, ...
    'VariableNames', {'index', 'mz', 'smiles'});

if ~isfile(fileName)
    return
end

lines = string(readlines(fileName));

% Trim surrounding whitespace, including the carriage return a Windows
% written file leaves at the end of every line.
lines = strip(lines);

blank = find(lines == "");

if isempty(blank)
    error('read_predicted_spectrum:noFragmentTable', ...
        ['%s has no blank line, so the fragment table cannot be located. ' ...
         'Expected a CFM-ID file: header, energy blocks, blank line, then ' ...
         'the fragment table.'], fileName);
end

% The table runs from just after the first blank line to just before the
% next one, or to the end of the file if there is no second blank.
firstRow = blank(1) + 1;
laterBlanks = blank(blank > firstRow);
if isempty(laterBlanks)
    lastRow = numel(lines);
else
    lastRow = laterBlanks(1) - 1;
end

if lastRow < firstRow
    return
end

block = lines(firstRow:lastRow);
block = block(block ~= "");

% Each row is "<index> <mz> <smiles>". The SMILES never contains a space,
% so splitting on the first two spaces is unambiguous.
parts = regexp(block, '^(\S+)\s+(\S+)\s+(.*)$', 'tokens', 'once');
valid = ~cellfun(@isempty, parts);

if ~any(valid)
    return
end

parts = vertcat(parts{valid});

fragments = table( ...
    str2double(parts(:, 1)), ...
    str2double(parts(:, 2)), ...
    strip(parts(:, 3)), ...
    'VariableNames', {'index', 'mz', 'smiles'});

end
