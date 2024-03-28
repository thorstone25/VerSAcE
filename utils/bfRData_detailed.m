%% Load FSA experimental data, visualize (channels, spectra), then beamform.
outdir = '/home/benfrey/Vantage-4.8.4/tmp';

% Load in mat files
% (a) The MAT file containing RData
% (b) The MAT file containing the QUPS configuration
% (c) The MAT file containing the VSX post-acquisition configuration
tic;
RData = load('/home/benfrey/Vantage-4.8.4/tmp/240322_122422_9876.mat').RData;
load('/home/benfrey/Vantage-4.8.4/MatFiles/qups-conf.mat', 'us', 'chd', 'QUPS_BF_PARAMS');
vs = load('/home/benfrey/Vantage-4.8.4/tmp/qups-vsx-post.mat');

% Take a look at the TGC information
lambda = us.seq.c0/(vs.Trans.frequency*1e6); % [m]
rangeM = vs.TGC.rangeMax * lambda; % [m]
figure();
xData = 0:rangeM/(512-1):rangeM;
plot(xData, vs.TGC.Waveform, LineWidth=3); % The waveform range has 512 elements.
xlim([0, xData(end)])
ylim([0, 1050]);
xticklabels(xticks*1e3)
xlabel("Axial [mm]")
ylabel("Gain [units?]")
title("Applied TGC amplitude over (axial?) depth");
grid on;
saveas(gcf, sprintf('%s/tgc_curve.png', outdir))

% Below is code from the bfQUPS function by Thurston
chd0 = chd;
[prms, chd0] = QUPS_RT_update_params(QUPS_BF_PARAMS, vs, chd0);
chd = QUPS_RT_load_data(RData, chd0, vs, prms);
disp("Data loaded in " + toc() + " seconds."); % loading time

% Take a look at the FFT of the channel data
elem = 1;
tmp = abs(fft(chd.data(:, elem, elem)));
figure; plot(db(tmp/max(tmp(:)))); title(sprintf('FFT of chd, elem %s', string(elem))); xlabel('Frequency [Hz?]'); ylabel('Amplitude'); xlim([0, length(tmp)])

% pre-process
tic;
chd = hilbert(filter(chd, prms.D));

% Take a look at the FFT of the filtered channel data
mid=1
tmp = abs(fft(chd.data(:, mid, mid)));
hold on; plot(db(tmp/max(tmp(:)))); hold off; legend(["Chd", "Filtered Chd"], "Location", "southwest")
saveas(gcf, sprintf('%s/fft_chd.png', outdir))

% get image
b = DAS(us, chd);
disp("Beamformed in " + toc() + " seconds."); % bf time
b = double(gather(abs(b)));
bnorm = max(b, [], 'all');
bimg = 20 * log10(b / bnorm);

% plot image
figure();
zmmRange = int32([us.scan.zb(1):10e-3:us.scan.zb(2)]*1e3);
xmmRange = int32([us.scan.xb(1):10e-3:us.scan.xb(2)]*1e3);
z_ticks_positions = linspace(1, size(b, 1), numel(zmmRange));
x_ticks_positions = linspace(1, size(b, 2), numel(xmmRange));
imagesc(bimg, 'XData', linspace(1, size(b, 2)), 'YData', linspace(1, size(b, 1)));
set(gca, 'XTick', x_ticks_positions, 'XTickLabel', xmmRange);
set(gca, 'YTick', z_ticks_positions, 'YTickLabel', zmmRange);
mminVal = min(bimg, [], 'all');
caxis([-70, 0])
colormap('bone');
cb = colorbar();
cb.Label.String = '[dB]';
ylabel("Axial [mm]")
xlabel("Lateral [mm]")
axis image
title('bfQUPS Example - CIRS-049A Type III Target (large circle)')