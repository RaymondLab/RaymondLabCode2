%% Load data
clear; loadAnalysisInfo_APP;

%% Parameters
lambda = 6;
freq = 1;
cutoff = 30;
samplerate = mag1.samplerate;

%% Preprocessing
posRaw = mag1.pos_data_aligned;
% posRaw = mag2.pos_data_aligned;
% posRaw = vid.pos_data_upsampled_aligned;

posFilt = butterworthfilter(posRaw, cutoff, samplerate);

% Differentiate position to calculate velocity
veltau = .01;
velRaw = movingslopeCausal(posRaw, round(samplerate*veltau))*samplerate;
velFilt = movingslopeCausal(posFilt, round(samplerate*veltau))*samplerate;

%% Fit Error Method
% Remove experiment frequency by calculating the error from initial fit guess
segLength = length(posFilt);
timeVec = 0:(1/samplerate):(segLength-1)/samplerate;

y1 = sin(2*pi*freq*timeVec(:));
y2 = cos(2*pi*freq*timeVec(:));
constant = ones(segLength,1);

vars = [y1, y2, constant];
keep = abs(velFilt) < 5*std(abs(velFilt)) + mean(abs(velFilt));
b = regress(velFilt(keep), vars(keep,:));
fit1 = vars *b;
times = 1:length(velFilt);

%% Wavelet Denoising
% mag1Vel_denoised = wden(mag1Vel, 'sqtwolog', 'h', 'mln', 6, 'sym8');
% [imf, residual] = emd(mag1Vel, 'MaxNumIMF',5);
% mag1Vel_denoised = imf(:,5);

%% Fast Iterative Filtering (FIF/EMD)
opts = Settings_FIF_v3('delta',0.1, ...
                       'NIMFs',11, ...
                       'Xi',4, ...
                       'alpha',24);
[imf, logM] = FIF_v2_12(velRaw, opts);
velImf = imf(5,:)';


%% Comparative Analysis
velFitErr = (velFilt - fit1);
velImfErr = (velFilt - velImf);

velRefMAD = lambda*movmad(velFilt, 1000);
velFitMAD = lambda*movmad(velFitErr, 1000);
velImfMAD = lambda*movmad(velImfErr, 1000);

velAbsMovMed = abs(velFilt - movmedian(velFilt, 1000));
velRefLocs = velAbsMovMed > velRefMAD;
velFitLocs = velAbsMovMed > velFitMAD;
velImfLocs = velAbsMovMed > velImfMAD;

velRefTimes = times(velRefLocs);
velFitTimes = times(velFitLocs);
velImfTimes = times(velImfLocs);

fig = figure('units', 'normalized', 'outerposition', [0.1, 0.1, 0.9, 0.9]);
h = tiledlayout(3, 1, "TileSpacing","compact", "Padding","compact");

% plot(velFilt - fit1, 'b');
% hold on;
% plot(velFilt - velImf', 'r');
% xlim([0, 5000]);
% ylim([-10, 10]);
% grid on;
% hold off;

ax1 = nexttile(1,[2,1]);
plot(times, posFilt, 'k');
hold on;
plot(velRefTimes, ones(length(velRefTimes)).*1.90, '.g');
plot(velFitTimes, ones(length(velFitTimes)).*1.87, '.b');
plot(velImfTimes, ones(length(velImfTimes)).*1.84, '.r');
xlim([0, round(length(velFilt)/4)]);
ylim([1.3, 2.0])
hold off;

ax2 = nexttile;
plot(times, velAbsMovMed,'k');
hold on;
plot(times, velRefMAD, 'g');
plot(times, velFitMAD, 'b');
plot(times, velImfMAD, 'r');
xlim([0, round(length(velFilt)/4)]);
ylim([0, 7]);
grid on;
hold off;

linkaxes([ax1, ax2], 'x');