function lcData = step6_normalise_fragments(lcData)
%STEP6_NORMALISE_FRAGMENTS  Scale fragments to the base peak and drop noise.
%
%   LCDATA = STEP6_NORMALISE_FRAGMENTS(LCDATA) normalizes peak intensity of 
%   every MS2 fragment to the peak intensity of the fragments with the highest 
%   intensity, and discards those below PARAMS.fragmentMinRelativeIntensity.
%
%   Normalising MS2 peak intensities are plotted in the GUI.
%
%   See also STEP5_CLEAN_MS2_FRAGMENTS, STEP7_MATCH_PREDICTED_SPECTRA.
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
    normalised = cell(nSamples, 1);

    for z = 1:nSamples
        normalisedHere = struct();

        for f = 1:numel(energyFields)
            field = energyFields{f};
            fragments = lcData(mode).fragments(z).(field);

            if isempty(fragments)
                normalisedHere.(field) = zeros(0, 2);
                continue
            end

            mz = fragments(:, 1);
            intensity = fragments(:, 2);
            relative = intensity / max(intensity);

            keep = relative >= params.fragmentMinRelativeIntensity;

            normalisedHere.(field) = [mz(keep), relative(keep)];
            lcData(mode).fragments(z).(field) = [mz(keep), intensity(keep)];
        end

        normalised{z} = normalisedHere;
    end

    lcData(mode).normalisedFragments = [normalised{:}];
end

end
