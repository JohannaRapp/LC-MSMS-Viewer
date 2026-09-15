function lcData = step5_clean_ms2_fragments(lcData)
%STEP5_CLEAN_MS2_FRAGMENTS  Clean up MS2 spectra.
%
%   LCDATA = STEP5_CLEAN_MS2_FRAGMENTS(LCDATA) removes two kinds of entry
%   from every picked MS2 spectrum:
%
%     * peaks at the precursor mass, which are precursor that are not fragmented in the
%       collision cell, and
%     * peaks having a larger mass than the precursor.
%
%   Removed fragments are kept in removedFragments, per collision energy, so
%   the original spectrum can be traced back.
%
%   The precursor comparison uses PARAMS.precursorMatchToleranceDa.
%
%   See also STEP6_NORMALISE_FRAGMENTS.
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
    energyFields = compose("CE%d", params.expectedCollisionEnergies);
    removed = cell(nSamples, 1);

    for z = 1:nSamples
        precursorMz = lcData(mode).precursorMz(z);
        removedHere = struct();

        for f = 1:numel(energyFields)
            field = energyFields{f};
            fragments = lcData(mode).fragments(z).(field);

            if isempty(fragments)
                removedHere.(field) = zeros(0, 2);
                continue
            end

            mz = fragments(:, 1);
            isPrecursor = abs(mz - precursorMz) < params.precursorMatchToleranceDa;
            isTooHeavy  = mz > precursorMz;
            drop = isPrecursor | isTooHeavy;

            removedHere.(field) = fragments(drop, :);
            lcData(mode).fragments(z).(field) = fragments(~drop, :);
        end

        removed{z} = removedHere;
    end

    lcData(mode).removedFragments = [removed{:}];
end

end
