% This file is meant to test matlab functions and how different data could
% change stuff

% Generate two signals with coherence at 20 Hz
% and compute wavelet coherence between 1-50 Hz

clear all;
close all;
clc;

%% Test coherence with different amplitudes
% Parameters
fs = 1000;              % Sampling frequency (Hz)
T = 10;                 % Total duration (seconds)
t = 0:1/fs:T-1/fs;      % Time vector
N = length(t);

% Target frequency for coherence
f_coherent = 20;        % Hz

% Linspace for interpolation
freq_lin = linspace(1, 200, 200);

% Generate the coherent component at 20 Hz (same amplitude)
amplitude = 1;
coherent_signal = amplitude * sin(2*pi*f_coherent*t);

% Generate two signals with the coherent component in the middle
% Add random noise to make signals realistic but keep coherence at 20 Hz

% Define the middle section (e.g., 30-70% of the signal)
middle_start = round(0.3 * N);
middle_end = round(0.7 * N);

% Initialize signals
signal1 = randn(1, N) * 0.3;  % Background noise
signal2 = randn(1, N) * 0.3;  % Background noise

% Add coherent component in the middle section
signal1(middle_start:middle_end) = signal1(middle_start:middle_end) + ...
    coherent_signal(middle_start:middle_end);
signal2(middle_start:middle_end) = signal2(middle_start:middle_end) + ...
    coherent_signal(middle_start:middle_end);

% Test flipping signal1
signal1 = -1 * signal1;

% Compute wavelet coherence
[wcoh, ~, freq] = wcoherence(signal1, signal2, fs);

% Plot the results
figure('Position', [100, 100, 1200, 800]);

% Plot Signal 1
subplot(4, 1, 1);
plot(t, signal1, 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal 1');
grid on;
xlim([0 T]);

% Plot Signal 2
subplot(4, 1, 2);
plot(t, signal2, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal 2');
grid on;
xlim([0 T]);

% Plot full wavelet coherence
subplot(4, 1, 3);
coh_interp = interp1(freq, abs(wcoh), freq_lin);
imagesc(t, freq_lin, coh_interp);
%surface(t, freq, abs(wcoh));
axis xy;
ylabel('Frequency (Hz)');
xlabel('Time (s)');
title('Wavelet Coherence (Full Range)');
colormap(jet);
colorbar;
caxis([0 1]);
ylim([1 50]);

% Highlight the 20 Hz line
hold on;
plot([0 T], [20 20], 'w--', 'LineWidth', 2);
hold off;

% Plot coherence at 20 Hz over time
subplot(4, 1, 4);
[~, idx_20Hz] = min(abs(freq - 20));
coherence_at_20Hz = wcoh(idx_20Hz, :);
plot(t, coherence_at_20Hz, 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Coherence');
title('Coherence at 20 Hz over Time');
grid on;
xlim([0 T]);
ylim([0 1]);

% Add shaded region for middle section
hold on;
patch([t(middle_start) t(middle_end) t(middle_end) t(middle_start)], ...
    [0 0 1 1], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold off;
legend('Coherence', 'Middle Section (High Coherence)');

% Display summary statistics
fprintf('\n=== Coherence Analysis Summary ===\n');
fprintf('Sampling Frequency: %d Hz\n', fs);
fprintf('Signal Duration: %.1f seconds\n', T);
fprintf('Target Coherent Frequency: %d Hz\n', f_coherent);
fprintf('\nCoherence at 20 Hz:\n');
fprintf('  Mean coherence (full signal): %.3f\n', mean(coherence_at_20Hz));
fprintf('  Mean coherence (middle section): %.3f\n', ...
    mean(coherence_at_20Hz(middle_start:middle_end)));
fprintf('  Max coherence: %.3f\n', max(coherence_at_20Hz));

% Change amplitude to one of the signals to check for difference in
% coherences
first_half_idxs = 1:round(length(signal1)/2);
signal1(first_half_idxs) = 0.5*signal1(first_half_idxs);

axs = [];

% Compute wavelet coherence
[wcoh, ~, freq] = wcoherence(signal1, signal2, fs);

% Plot the results
figure('Position', [100, 100, 1200, 800]);

% Plot Signal 1
subplot(4, 1, 1);
plot(t, signal1, 'b', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal 1');
grid on;
xlim([0 T]);
axs(end + 1) = gca;

% Plot Signal 2
subplot(4, 1, 2);
plot(t, signal2, 'r', 'LineWidth', 1);
xlabel('Time (s)');
ylabel('Amplitude');
title('Signal 2');
grid on;
xlim([0 T]);
axs(end + 1) = gca;

% Plot full wavelet coherence
subplot(4, 1, 3);
coh_interp = interp1(freq, abs(wcoh), freq_lin);
imagesc(t, freq_lin, coh_interp);
%surface(t, freq, abs(wcoh));
axis xy;
ylabel('Frequency (Hz)');
xlabel('Time (s)');
title('Wavelet Coherence with different amplitude');
colormap(jet);
colorbar;
caxis([0 1]);
ylim([1 50]);
axs(end + 1) = gca;

% Highlight the 20 Hz line
hold on;
plot([0 T], [20 20], 'w--', 'LineWidth', 2);
hold off;

% Plot coherence at 20 Hz over time
subplot(4, 1, 4);
[~, idx_20Hz] = min(abs(freq - 20));
coherence_at_20Hz = wcoh(idx_20Hz, :);
plot(t, coherence_at_20Hz, 'k', 'LineWidth', 2);
xlabel('Time (s)');
ylabel('Coherence');
title('Coherence at 20 Hz over Time');
grid on;
xlim([0 T]);
ylim([0 1]);
axs(end + 1) = gca;

% Add shaded region for middle section
hold on;
patch([t(middle_start) t(middle_end) t(middle_end) t(middle_start)], ...
    [0 0 1 1], 'g', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
hold off;
legend('Coherence', 'Middle Section (High Coherence)');

linkaxes(axs, 'x');

% Save the results
%save('coherence_results.mat', 'signal1', 'signal2', 'wcoh', 'freq', 't', 'fs');
%fprintf('\nResults saved to coherence_results.mat\n');

%%