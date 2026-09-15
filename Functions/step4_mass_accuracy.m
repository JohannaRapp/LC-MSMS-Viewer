function lcData = step4_mass_accuracy(lcData)
%STEP4_MASS_ACCURACY  Mass deviation of isolated mass from the expected
%one.
%
%   LCDATA = STEP4_MASS_ACCURACY(LCDATA) measures the distance between the
%   expected m/z.value and the most intense peak present in that isolation window,
%   and records it as purity(z).massDeviationDa.
%
%   A large deviation means the highest peak in the window is not the
%   compound the method was aiming at. Read together with the purity
%   check from step 3 it separates two different failures: a clean but misplaced
%   isolation, and a correctly placed one contaminated by a neighbour.
%
%   See also STEP3_PRECURSOR_PURITY.
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

    for z = 1:nSamples
        precursorMz = lcData(mode).precursorMz(z);
        window = isolation_window(lcData(mode).ms1(z).mz, precursorMz, ...
            params.isolationHalfWidthDa);

        windowMz  = double(lcData(mode).ms1(z).mz(window));
        windowInt = double(lcData(mode).ms1(z).intensity(window));

        [~, tallest] = max(windowInt);
        lcData(mode).purity(z).massDeviationDa = abs(precursorMz - windowMz(tallest));
    end
end

end
