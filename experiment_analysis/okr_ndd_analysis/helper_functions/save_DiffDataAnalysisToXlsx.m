function save_DiffDataAnalysisToXlsx(timepoints, diffdata, params)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

ntimepoints = length(timepoints);
ndiffdata = length(diffdata);
nmetricrows = ntimepoints + ndiffdata;
rownum = 0;

% Initialize metrics "Type", "BlockNums", and "BaseAmp" columns for table
metrics.Type = cell(nmetricrows, 1);
metrics.BlockNums = cell(nmetricrows, 1);
metrics.DrumAmp = zeros(nmetricrows, 1);
metrics.DrumAmp_SEM = zeros(nmetricrows, 1);
metrics.DrumPhase = zeros(nmetricrows, 1);
metrics.DrumPhase_SEM = zeros(nmetricrows, 1);
metrics.RelGain = zeros(nmetricrows, 1);
metrics.RelGain_SEM = zeros(nmetricrows, 1);
metrics.RelPhase = zeros(nmetricrows, 1);
metrics.RelPhase_SEM = zeros(nmetricrows, 1);

% Nasal-Temporal (NT) half-cycle peak times, centroid, lag, and skew metrics
metrics.NT_PeakTime1 = zeros(nmetricrows, 1);
metrics.NT_PeakTime2 = zeros(nmetricrows, 1);
metrics.NT_Mean = zeros(nmetricrows, 1);
metrics.NT_MeanLag = zeros(nmetricrows, 1);
metrics.NT_MeanSkew = zeros(nmetricrows, 1);
metrics.NT_Median = zeros(nmetricrows, 1);
metrics.NT_MedianLag = zeros(nmetricrows, 1);
metrics.NT_MedianSkew = zeros(nmetricrows, 1);
metrics.NT_BowleySkew = zeros(nmetricrows, 1);
metrics.NT_MeanMedianSkew = zeros(nmetricrows, 1);

% Temporal-Nasal (TN) half-cycle peak times, centroid, lag, and skew metrics
metrics.TN_PeakTime1 = zeros(nmetricrows, 1);
metrics.TN_PeakTime2 = zeros(nmetricrows, 1);
metrics.TN_Mean = zeros(nmetricrows, 1);
metrics.TN_MeanLag = zeros(nmetricrows, 1);
metrics.TN_MeanSkew = zeros(nmetricrows, 1);
metrics.TN_Median = zeros(nmetricrows, 1);
metrics.TN_MedianLag = zeros(nmetricrows, 1);
metrics.TN_MedianSkew = zeros(nmetricrows, 1);
metrics.TN_BowleySkew = zeros(nmetricrows, 1);
metrics.TN_MeanMedianSkew = zeros(nmetricrows, 1);

% Initialize diffdata struct for table
data.CycleTimes = (0:length(timepoints(1).pooled_drumvel_cyclemean_cvd.cycleMean)-1)' / params.fs;
data.StimVel = timepoints(1).pooled_drumvel_cyclemean_cvd.cycleMean';

for ii = 1:ntimepoints
    rownum = rownum + 1;
    data.(timepoints(ii).timePoint) = timepoints(ii).pooled_good_cyclemean_cvd.cycleMean';
    data.([timepoints(ii).timePoint,'_SEM']) = timepoints(ii).pooled_good_cyclemean_cvd.cycleSEM';
    data.([timepoints(ii).timePoint,'_Fit']) = timepoints(ii).pooled_good_cyclemean_cvd.fitMean';
    resii = timepoints(ii).pooled_good_cyclemean_cm;  % Retrieve cycle metric results

    DrumAmp = timepoints(ii).pooled_drumvel_cyclemean_cvd.amplitudeMean;
    DrumAmp_SEM = timepoints(ii).pooled_drumvel_cyclemean_cvd.amplitudeSEM;
    EyeAmp = timepoints(ii).pooled_good_cyclemean_cvd.amplitudeMean;
    EyeAmp_SEM = timepoints(ii).pooled_good_cyclemean_cvd.amplitudeSEM;
    [RelGain,RelGain_SEM] = calc_eyeRelGainSEM(EyeAmp, EyeAmp_SEM, ...
            DrumAmp, DrumAmp_SEM);

    DrumPhase = timepoints(ii).pooled_drumvel_cyclemean_cvd.phaseMean_deg;
    DrumPhase_SEM = timepoints(ii).pooled_drumvel_cyclemean_cvd.phaseSEM_deg;
    EyePhase = timepoints(ii).pooled_good_cyclemean_cvd.phaseMean_deg;
    EyePhase_SEM = timepoints(ii).pooled_good_cyclemean_cvd.phaseSEM_deg;
    [RelPhase,RelPhase_SEM] = calc_eyeRelPhaseSEM(EyePhase, EyePhase_SEM, ...
            DrumPhase, DrumPhase_SEM);

    metrics.Type{rownum} = timepoints(ii).timePoint;
    metrics.BlockNums{rownum} = strrep(timepoints(ii).blockNumbers, ', ', ';');

    metrics.DrumAmp(rownum) = DrumAmp;
    metrics.DrumAmp_SEM(rownum) = DrumAmp_SEM;
    metrics.DrumPhase(rownum) = DrumPhase;
    metrics.DrumPhase_SEM(rownum) = DrumPhase_SEM;
    metrics.RelGain(rownum) = RelGain;
    metrics.RelGain_SEM(rownum) = RelGain_SEM;
    metrics.RelPhase(rownum) = RelPhase;
    metrics.RelPhase_SEM(rownum) = RelPhase_SEM;

    metrics.NT_PeakTime1(rownum) = resii.eye.peak1TimeMs(1);
    metrics.NT_PeakTime2(rownum) = resii.eye.peak1TimeMs(2);
    metrics.NT_Mean(rownum) = resii.eye.centroidMean(1);
    metrics.NT_MeanLag(rownum) = resii.centroidMeanDiff(1);
    metrics.NT_MeanSkew(rownum) = resii.eye.skewFromCentroidMean(1);
    metrics.NT_Median(rownum) = resii.eye.centroidMedian(1);
    metrics.NT_MedianLag(rownum) = resii.centroidMedianDiff(1);
    metrics.NT_MedianSkew(rownum) = resii.eye.skewFromCentroidMedian(1);
    metrics.NT_BowleySkew(rownum) = resii.eye.skewQuantile(1);
    metrics.NT_MeanMedianSkew(rownum) = resii.eye.skewMeanMedian(1);
    
    metrics.TN_PeakTime1(rownum) = resii.eye.peak2TimeMs(1);
    metrics.TN_PeakTime2(rownum) = resii.eye.peak2TimeMs(2);
    metrics.TN_Mean(rownum) = resii.eye.centroidMean(2);
    metrics.TN_MeanLag(rownum) = resii.centroidMeanDiff(2);
    metrics.TN_MeanSkew(rownum) = resii.eye.skewFromCentroidMean(2);
    metrics.TN_Median(rownum) = resii.eye.centroidMedian(2);
    metrics.TN_MedianLag(rownum) = resii.centroidMedianDiff(2);
    metrics.TN_MedianSkew(rownum) = resii.eye.skewFromCentroidMedian(2);
    metrics.TN_BowleySkew(rownum) = resii.eye.skewQuantile(2);
    metrics.TN_MeanMedianSkew(rownum) = resii.eye.skewMeanMedian(2);
end


for ii = 1:ndiffdata
    rownum = rownum + 1;
    fname = ['Diff_', strrep(diffdata(ii).timePointDiff,'-','_')];
    data.(fname) = diffdata(ii).pooled_good_cyclemean_diff';
    data.([fname,'_SEM']) = diffdata(ii).pooled_good_cyclemean_diff_SEM';
    resii = diffdata(ii).pooled_good_cyclemean_diff_cm;  % Retrieve cycle metric results

    metrics.Type{rownum} = fname;
    metrics.BlockNums{rownum} = "NA";

    metrics.DrumAmp(rownum) = nan;
    metrics.DrumAmp_SEM(rownum) = nan;
    metrics.DrumPhase(rownum) = nan;
    metrics.DrumPhase_SEM(rownum) = nan;
    metrics.RelGain(rownum) = nan;
    metrics.RelGain_SEM(rownum) = nan;
    metrics.RelPhase(rownum) = nan;
    metrics.RelPhase_SEM(rownum) = nan;
    
    metrics.NT_PeakTime1(rownum) = resii.eye.peak1TimeMs(1);
    metrics.NT_PeakTime2(rownum) = resii.eye.peak1TimeMs(2);
    metrics.NT_Mean(rownum) = resii.eye.centroidMean(1);
    metrics.NT_MeanLag(rownum) = resii.centroidMeanDiff(1);
    metrics.NT_MeanSkew(rownum) = resii.eye.skewFromCentroidMean(1);
    metrics.NT_Median(rownum) = resii.eye.centroidMedian(1);
    metrics.NT_MedianLag(rownum) = resii.centroidMedianDiff(1);
    metrics.NT_MedianSkew(rownum) = resii.eye.skewFromCentroidMedian(1);
    metrics.NT_BowleySkew(rownum) = resii.eye.skewQuantile(1);
    metrics.NT_MeanMedianSkew(rownum) = resii.eye.skewMeanMedian(1);
    
    metrics.TN_PeakTime1(rownum) = resii.eye.peak2TimeMs(1);
    metrics.TN_PeakTime2(rownum) = resii.eye.peak2TimeMs(2);
    metrics.TN_Mean(rownum) = resii.eye.centroidMean(2);
    metrics.TN_MeanLag(rownum) = resii.centroidMeanDiff(2);
    metrics.TN_MeanSkew(rownum) = resii.eye.skewFromCentroidMean(2);
    metrics.TN_Median(rownum) = resii.eye.centroidMedian(2);
    metrics.TN_MedianLag(rownum) = resii.centroidMedianDiff(2);
    metrics.TN_MedianSkew(rownum) = resii.eye.skewFromCentroidMedian(2);
    metrics.TN_BowleySkew(rownum) = resii.eye.skewQuantile(2);
    metrics.TN_MeanMedianSkew(rownum) = resii.eye.skewMeanMedian(2);
end

% Convert the structs to tables
DT = struct2table(data);
MT = struct2table(metrics);

% Write this data to new sheets in excel sheet
[folder,file,~] = fileparts(params.save_filepath);
xlsx_filepath = fullfile(folder, [file, '.xlsx']);
writetable(DT, xlsx_filepath, 'Sheet','DiffData');
writetable(MT, xlsx_filepath, 'Sheet','DiffInfo');


end