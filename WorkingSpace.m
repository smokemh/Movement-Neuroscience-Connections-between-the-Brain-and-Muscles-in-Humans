% Course : Movement Neuroscience: Connections between the Brain and Muscles in Humans
% Description : Various surface EMG data postprocessing and analysis steps.
% FAU Erlangen
%% Matlab script settings
clc;
clear all;
close all;
format compact;
rng('default');
warning('off', 'MATLAB:dispatcher:InexactMatch');
format long;
customColors = [0.2 0.4 0.6; 
                0.8 0.1 0.3;  
                0.1 0.6 0.2];
customFontName = 'Arial';
customFontSize = 12;
customFontWeight = 'bold';
customLineWidth = 2;
legendLocation = 'northeast';
customFigureWidth = 800;
customFigureHeight = 500;
screenSize = get(0, 'ScreenSize');
centerX = (screenSize(3) - customFigureWidth) / 2;
centerY = (screenSize(4) - customFigureHeight) / 2;
customFigure = [centerX, centerY, customFigureWidth, customFigureHeight];
set(groot, 'defaultAxesColorOrder', customColors);
set(groot, 'defaultLineLineWidth', customLineWidth);
set(groot, 'defaultAxesFontName', customFontName);
set(groot, 'defaultAxesFontSize', customFontSize);
set(groot, 'defaultAxesFontWeight', customFontWeight);
set(groot, 'defaultLegendLocation', legendLocation);
set(groot, 'defaultAxesXGrid', 'on');
set(groot, 'defaultAxesYGrid', 'on');
set(groot, 'defaultAxesZGrid', 'on');
%% Pre-Processing the data and visulizing the preprocessed data
slow_contraction = load('Slow_Contraction.mat');
rapid_contraction = load('Rapid_Contractions.mat');

% Plot for slow contraction reference signal
figure('Position', customFigure); 
subplot(2,1,1)
plot(slow_contraction.ref_signal);
title('Slow contraction reference signal');
xlabel('Time');
ylabel('Force in mV');
legend('Slow Reference signal');
subplot(2,1,2)
plot(slow_contraction.SIG{1});
title('Slow contraction SIG');
ylabel('Voltage in V');
legend('Slow SIG signal');

%Plot for rapid contraction reference signal
figure('Position', customFigure);
subplot(2,1,1)
plot(rapid_contraction.ref_signal);
title('Rapid contraction sampling frequence');
xlabel('Time' );
ylabel('Force in mV)');
legend('Rapid Reference Signal');
subplot(2,1,2)
plot(rapid_contraction.SIG{1});
title('Rapid contraction SIG');
ylabel('Voltage in V');
legend('Rapid SIG signal');
%% 1.1 : Tranforming reference signal to Newton
C = 0.02;
g = 9.81;
slow_contraction.signal=slow_contraction.ref_signal/C*g;

%Plot for converted signal with respect to orignal
figure('Position', customFigure);
subplot(2,1,1);
plot(slow_contraction.ref_signal);
title('Slow reference signal');
xlabel('Time');
ylabel('Force in mV');
legend('reference signal');
subplot(2,1,2);
plot(slow_contraction.signal);
title('Slow transformed signal');
xlabel('Time');
ylabel('Force in N');
legend('Transformed Signal');
%% 1.2 : Getting time axis
dt = 1 / slow_contraction.fsamp;
T = (numel(slow_contraction.signal) - 1) * dt;
t_s = 0:dt:T;
figure('Position', customFigure);
plot(t_s,slow_contraction.signal);
title('Slow Signal in Time');
xlabel('Time in s');
ylabel('Force in N');
%% 1.3 : Applying low Pass Filter
% Design Butterworth low-pass filter
cutoff_frequency = 20;
order = 5;
[b_butter, a_butter] = butter(order, cutoff_frequency / (slow_contraction.fsamp / 2));
slow_contraction.filterd_signal = filtfilt(b_butter, a_butter, slow_contraction.signal);

% Plot the original and filtered signals
figure('Position', customFigure);
plot(t_s, slow_contraction.signal);
hold on;
plot(t_s, slow_contraction.filterd_signal);
hold off;
title('Unfiltered v.s Filtered');
xlabel('Time in s');
ylabel('Force in N');
legend('Unfiltered Signal', 'Filtered Signal');
%% 1.4 : EMG data with signal
figure ('Position', customFigure);
plot(t_s, slow_contraction.filterd_signal);
xlabel('Time in s');
ylabel('Force in N');
yyaxis right;
plot(t_s, slow_contraction.SIG{1,1});
ylabel('Voltage in V');
title('EMG data with signal');

% Q: How do force and EMG signal relate?
% A: When force increases EMG also increase
%% 2.1 : Computing CV
Sig = slow_contraction.filterd_signal;
[~, p_start] = max(Sig);
p_X = Sig(p_start);
Sig(max(1, p_start - 1000):min(end, p_start + 1000)) = 0;
p_end = find(Sig == max(Sig), 1);
p_Y = Sig(p_end);
fprintf('Plateau start value = %d\n', p_start);
fprintf('Plateau end value = %d\n', p_end);

% Ploting the Signal 
figure('Position',customFigure);
plot(slow_contraction.filterd_signal, 'DisplayName', 'Reference Signal');
hold on;
plot(p_start, p_X, 'o', 'DisplayName', 'Plateau Start');
plot(p_end, p_Y, 'o', 'DisplayName', 'Plateau End');
text(p_start, p_X, ['  Plateau Start = ' num2str(p_start)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
text(p_end, p_Y, ['  Plateau End = ' num2str(p_end)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
title('Filtered Signal with Identified Plateau');
xlabel('Time in s');
ylabel('Force in N');
legend('Show');

%Compute CV
Sig_p = slow_contraction.filterd_signal(p_start:p_end);
CV = (std(Sig_p)/mean(Sig_p)) * 100;
fprintf('CV = %f\n', CV);
%% 2.2 : Q: What Does it mean physiologically? 
%Ans : In Report
%% 3.1 : Pre-processing rapid contraction
rapid_contraction.signal=rapid_contraction.ref_signal/C*g;

%Plot for converted signal with respect to orignal
figure('Position', customFigure);
subplot(2,1,1);
plot(rapid_contraction.ref_signal);
title('Rapid reference signal');
xlabel('Time');
ylabel('Force in mV');
legend('reference signal');
subplot(2,1,2);
plot(rapid_contraction.signal);
title('Rapid transformed signal');
xlabel('Time');
ylabel('Force in N');
legend('Transformed Signal');

dt = 1 / rapid_contraction.fsamp;
T = (numel(rapid_contraction.signal) - 1) * dt;
t_r = 0:dt:T;
figure('Position', customFigure);
plot(t_r,rapid_contraction.signal);
title('Rapid Signal in Time');
xlabel('Time in s');
ylabel('Force in N');

% Design Butterworth low-pass filter
cutoff_frequency = 20;
order = 5;
[b_butter, a_butter] = butter(order, cutoff_frequency / (rapid_contraction.fsamp / 2));
rapid_contraction.filterd_signal = filtfilt(b_butter, a_butter, rapid_contraction.signal);

% Plot the original and filtered signals
figure('Position', customFigure);
plot(t_r, rapid_contraction.signal);
hold on;
plot(t_r, rapid_contraction.filterd_signal);
hold off;
title('Unfiltered v.s Filtered');
xlabel('Time in s');
ylabel('Force in N');
legend('Unfiltered Signal', 'Filtered Signal');

% For EMG data
figure ('Position', customFigure);
plot(t_r, rapid_contraction.filterd_signal);
xlabel('Time in s');
ylabel('Force in N');
yyaxis right;
plot(t_r, rapid_contraction.SIG{1,1});
ylabel('Voltage in V');
title('EMG data with signal');
%% 3.2 Computation of RFD
td = 0.1 * max(rapid_contraction.filterd_signal);
onset_p = find(rapid_contraction.filterd_signal > td);
minDistance = 100;
onset_p = onset_p([true, diff(onset_p) > minDistance]);
disp('Detected Onset Indices:');
disp(onset_p);

% Plot to see if values are correct
figure('Position', customFigure);
t_samp = (0:(length(rapid_contraction.filterd_signal) - 1)) / rapid_contraction.fsamp;
plot(t_samp, rapid_contraction.filterd_signal, 'DisplayName', 'Filtered Signal');
hold on;
plot(t_samp(onset_p), rapid_contraction.filterd_signal(onset_p), 'o', 'DisplayName', 'Detected Onset');
title('Filtered signal with onset points');
xlabel('Time in s');
ylabel('Force in N');
legend('show');
for i = 1:length(onset_p)
    text(t_samp(onset_p(i)), rapid_contraction.filterd_signal(onset_p(i)), num2str(onset_p(i)), 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

% Calculating  instantaneous RFD
S_window = 0.5 * rapid_contraction.fsamp;
F_window = arrayfun(@(start_idx) rapid_contraction.filterd_signal(start_idx:start_idx + S_window - 1), onset_p, 'UniformOutput', false);
instant_RFD = cellfun(@(force) diff(force) / (1 / rapid_contraction.fsamp), F_window, 'UniformOutput', false);

% Visualize contractions with their respective RFD values
%Plot 1: 5 Contraction Visulaize with RFD
figure('Position', customFigure);
t_vector = (1:size(instant_RFD{1}, 2)) / rapid_contraction.fsamp;
for i = 1:5
    subplot(5, 1, i)
    plot(t_vector, instant_RFD{i})
    ylabel('N/s')
    yyaxis right
    plot(t_vector, F_window{i}(1:end-1),'--')
    ylabel('Force in N')
    xlabel('Time in s')
    title("Contraction Number:  " + num2str(i))
end
sgtitle('Contractions with Respective RFD Values', 'FontSize', 14);
%Plot 2: Histogram with 1 contraction
data_hist = instant_RFD{1};
figure('Position', customFigure);
subplot(2, 1, 1);
plot(t_vector, instant_RFD{1})
ylabel('N/s')
yyaxis right
plot(t_vector, F_window{1}(1:end-1),'--')
ylabel('Force in N')
xlabel('Time in s')
title("Contraction Number: 1")
subplot(2, 1, 2);
histogram(data_hist, 20);
xlabel('RFD Bins');
ylabel('Frequency');
title('Histogram of RFD');
sgtitle('Contractions with Respective RFD Values and Histogram', 'FontSize', 14);
%% 3.3 What is the physiological relevance? 

% Ans: In Report

%% 4.1: SIG to 2D Array
slow_contraction.EMG = mean(cat(3, slow_contraction.SIG{:}), 3);
rapid_contraction.EMG = mean(cat(3, rapid_contraction.SIG{:}), 3);
%% 4.2 : Compute the RMS 200 ms GIVEN
win_leng=(200/1000);
win = round(win_leng * slow_contraction.fsamp);
movWindow = ones(1, win) / win;
convrms_slow = sqrt(conv(slow_contraction.EMG.^2, movWindow, 'same'));
convrms_rapid = sqrt(conv(rapid_contraction.EMG.^2, movWindow, 'same'));
dt_slow = 1 / slow_contraction.fsamp;
t_slow = 0:dt_slow:(size(slow_contraction.signal, 2) - 1) * dt_slow;
dt_rapid = 1 / rapid_contraction.fsamp;
t_rapid = 0:dt_rapid:(size(rapid_contraction.signal, 2) - 1) * dt_rapid;

% To Plot
figure ('Position', customFigure);
plot(t_slow, slow_contraction.EMG);
hold on;
plot(t_slow, convrms_slow);
title('Slow Contraction 200ms');
xlabel('Time in s)');
figure ("Position", customFigure);
plot(t_rapid, rapid_contraction.EMG);
hold on;
plot(t_rapid, convrms_rapid);
title('Rapid Contraction 200ms');
xlabel('Time in s');

%% 4.3 : Correlation and Scatter Plot
mean_conv = mean(convrms_slow);
mean_signal = mean(slow_contraction.filterd_signal);
numerator = sum((convrms_slow - mean_conv) .* (slow_contraction.filterd_signal - mean_signal));
denominator_conv = sqrt(sum((convrms_slow - mean_conv).^2));
denominator_signal = sqrt(sum((slow_contraction.filterd_signal - mean_signal).^2));
R_s = numerator / (denominator_conv * denominator_signal);
fprintf('Coefficient R in Slow Data = %.5f\n', R_s);

%ScatterPlot
figure('Position', customFigure);
scatter(convrms_slow/mean(convrms_slow), slow_contraction.filterd_signal/mean(slow_contraction.filterd_signal))
hold on
plot([0 4], [0 4]); 
xlim ([0 2]);
ylim([0 2])
title('Scatterplot of Slow Data')
xlabel('RMS')
ylabel('Signal')

% Correlation for rapid data
mean_conv_rapid = mean(convrms_rapid);
mean_signal_rapid = mean(rapid_contraction.filterd_signal);
numerator_rapid = sum((convrms_rapid - mean_conv_rapid) .* (rapid_contraction.filterd_signal - mean_signal_rapid));
denominator_conv_rapid = sqrt(sum((convrms_rapid - mean_conv_rapid).^2));
denominator_signal_rapid = sqrt(sum((rapid_contraction.filterd_signal - mean_signal_rapid).^2));
R_rapid = numerator_rapid / (denominator_conv_rapid * denominator_signal_rapid);
fprintf('Coefficient R in Rapid Data = %.5f\n', R_rapid);

% ScatterPlot 
figure('Position', customFigure);
scatter(convrms_rapid/mean(convrms_rapid), rapid_contraction.filterd_signal/mean(rapid_contraction.filterd_signal))
hold on
plot([0 6], [0 6]);  
xlim([0 3.5])
ylim([0 3.5])
title('Scatterplot Rapid Data')
xlabel('RMS')
ylabel('Signal')




