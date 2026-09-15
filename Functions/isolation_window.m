function window = isolation_window(mz, precursorMz, halfWidthDa)
%ISOLATION_WINDOW  Index range of the quadrupole isolation window.
%
%   WINDOW = ISOLATION_WINDOW(MZ, PRECURSORMZ, HALFWIDTHDA) returns the
%   indices of MZ lying within HALFWIDTHDA of PRECURSORMZ, as the range the
%   quadrupole let through on its way to the collision cell.
%
%   See also STEP3_PRECURSOR_PURITY, STEP4_MASS_ACCURACY.
%
%   Part of the MATLAB app "LC_MSMS_Viewer_V2.mlapp".
%   Authors: Hannes Link, Johanna Rapp, 10.09.2026.

arguments
    mz (:,1)
    precursorMz (1,1) double
    halfWidthDa (1,1) double
end

mz = double(mz);

[~, firstIndex] = min(abs(mz - (precursorMz - halfWidthDa)));
[~, lastIndex]  = min(abs(mz - (precursorMz + halfWidthDa)));

window = firstIndex:lastIndex;

end
