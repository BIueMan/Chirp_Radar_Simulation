clear all
close all
addpath('objects');
addpath('func');
addpath('test_func');

%% set data
C = physconst('LightSpeed');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_transmit_start = -1e-5;
time_record_end = 7e-4;
f = 1e11;
T = 1/f;

f_0 = 76e9;
f_1 = f_0 + 3.4e9;
T_sample = 7.031256e-6;
time = 0: T_sample: time_record_end;
chirp_len = 125e-6;
c = (f_1-f_0)/chirp_len;
effective_fft_size = 1024*128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set source and target
source_locs = [0,0,0];
source = Element(source_locs, @chirp, {f_0, chirp_len, f_1, 'complex'});
source.SavedSignals = source.get_signal(time);
args_dict = containers.Map({'Pt', 'Gt', 'sigma', 'Ae'}, ...
                           { 1  ,    1,  1     ,  1  });
target_properties = 1;

target_locs = [0,0.1,0];
target = Element(target_locs, NaN, NaN);


%% get signal
[r_signal, true_dist_travel] = source.to_locs1_to_locs2(target.Locs, source.Locs, time, args_dict);
r_signal = reshape(r_signal, [1, length(r_signal)]);
s_signal = source.get_signal(time);

%% find dist
% aprrox dist
dist_travel = Find.dist_travel(s_signal, r_signal, c, T_sample, effective_fft_size);

dist_travel/2
true_dist_travel/2

%% plot signals and mixed
% figure
% subplot(2,1,1);
% stft(source.get_signal(time))
% subplot(2,1,2);
% stft(squeeze(signal(1,1,:)).')
% figure
% stft(source.get_signal(time).*conj(squeeze(signal(1,1,:)).'))

%% test range estimation
range = 0.0001:0.0001:0.07;
target_locs = [zeros(size(range)); zeros(size(range)); range].';
% % % test_dist_changes(time, source, target_locs, effective_fft_size);
%% test SNR
repet = 100;
SNR = 50;
%[range_true, var_error] = test_radar_SNR(time, source, target_locs, effective_fft_size, SNR, repet);

SNRdb = 10:5:80;
SNRlin = 10.^(SNRdb./10);
%SNR = 10:10:400;
range = 0.02:0.01:0.07;
target_locs = [zeros(size(range)); zeros(size(range)); range].';
[range_true, var_error] = test_SNR(time, source, target_locs, effective_fft_size, SNRlin, repet);

% %% SNR
% range = 0.001:0.001:2;
% target_locs = [zeros(size(range)); zeros(size(range)); range].';
% N = 1;
% test_radar_SNR(time, source, target_locs, args_dict, N, 1);


% %% SNR equation
% Pr_func = @(r) (args_dict('Pt')*args_dict('Gt')*args_dict('sigma')*args_dict('Ae')) ./...
%                ((4*pi)^2 .* r.^4);
% %% get SNR vs R plot
% noise = 1;
% snr_func = @(r) (args_dict('Pt')*args_dict('Gt')*args_dict('sigma')*args_dict('Ae')) ./ ...
%                 ((4*pi)^2 .* r.^4 .* noise);
% r = 0.1:0.01:200;
% figure
% plot(r, snr_func(r))
% title('SNR vs dist')

%% add noise to the signal
