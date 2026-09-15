function lcData = step2_fold_changes(lcData)
%STEP2_FOLD_CHANGES  Compare peak intensity of each targeted metabolite with the 
%   peak intensity of this metabolite in all other samples.
%
%   LCDATA = STEP2_FOLD_CHANGES(LCDATA) calculate fold chnages of all
%   target metabolites compared to the median peak intensity of all
%   samples.
%
%   The comparison is made at the retention time where the metabolite eluted
%   in the sample where the metabolite was targeted. Every other sample of 
%   the same polarity has a chromatogram for the same mass, so their intensities 
%   at that moment can be read off directly; retention times drift slightly 
%   between runs, so the others are interpolated onto this sample's time axis first. 
%   The fold change is the sample's peak height divided by the median peak height of all of them
%   at that same timepoint.
%
%   Each polarity gains
%
%       foldChange       one value per sample
%       eic(z).reference retention time axis the comparison was made on
%       eic(z).interpolatedIntensity  every sample on that axis
%       eic(z).medianIntensity        median across samples, per time point
%       eic(z).sdIntensity            standard deviation, per time point
%
%   The median and standard deviation are what the app draws behind the
%   chromatogram.
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

    lcData(mode).foldChange = nan(nSamples, 1);

    for z = 1:nSamples
        apex = lcData(mode).apex(z);
        column = apex.targetColumn;

        % The sample's own chromatogram is the reference axis.
        referenceTime = lcData(mode).eic(z).retentionTime(:, column);
        ownIntensity  = lcData(mode).eic(z).intensity(:, column);

        interpolated = nan(numel(referenceTime), nSamples);

        for other = 1:nSamples
            otherTime      = lcData(mode).eic(other).retentionTime(:, column);
            otherIntensity = lcData(mode).eic(other).intensity(:, column);
            interpolated(:, other) = interp1(otherTime, otherIntensity, ...
                referenceTime, 'linear');
        end

        medianIntensity = median(interpolated, 2, 'omitmissing');
        sdIntensity     = std(interpolated, [], 2, 'omitmissing');

        atApex = apex.rowInEic;
        lcData(mode).foldChange(z) = ownIntensity(atApex) / medianIntensity(atApex);

        lcData(mode).eic(z).referenceTime          = referenceTime;
        lcData(mode).eic(z).interpolatedIntensity  = interpolated;
        lcData(mode).eic(z).medianIntensity        = medianIntensity;
        lcData(mode).eic(z).sdIntensity            = sdIntensity;
    end
end

end
