clear all
close all
addpath('objects');
addpath('func');
addpath('test_func');

%% set data
C = physconst('LightSpeed');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_transmit_start = -1e-8;
time_record_end = 7e-4;
T_sample = 7.031256e-6;
upsample = 1e-6;
T = T_sample * upsample;
time = [time_transmit_start: T: 0, 0: T: time_record_end];

f_0 = 76e9;
f_1 = f_0 + 4e9;
amp = 1;
chirp_len = 125e-6;
c = (f_1-f_0)/chirp_len;
effective_fft_size = 1024*128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% time_transmit_start = 0; %-1e-5;
% time_record_end = 2e-6;
% T = 1e-10;
% time = time_transmit_start: T: time_record_end;
% 
% f_0 = 3e9;
% f_1 = f_0 + 0.4e9;
% amp = 1;
% chirp_len = 2e-8;
% c = (f_1-f_0)/chirp_len;
% effective_fft_size = 1024*128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set source and target
source_locs = [0,0,0];
source = Element(source_locs, @chirp, {f_0, chirp_len, f_1, 'complex'});
source.SavedSignals = source.get_signal(time);
args_dict = containers.Map({'Pt', 'Gt', 'sigma', 'Ae'}, ...
                           { 1  ,    1,  1     ,  1  });
target_properties = 1;

target_locs = [0,0,0.2];
target = Element(target_locs, NaN, NaN);
%% calculate signal travel, from the first element to the traget and beck
% source to target
[signal, dist] = source.to_locs(target_locs, time, target_properties);
target.SavedSignals = squeeze(signal(1,:,:)).';
s2t_dist = dist(1,:,:);
% target to source
[signal, dist] = target.to_locs(source_locs, time, target_properties);
r_signal = reshape(signal, [1,size(signal, 3)]);
t2s_dist = squeeze(dist);

%% get source signal from the first element, and extract dist
% get source signal, remove negative time, and downsampling
[~,idx] = min(abs(time));
time_positive = time(idx:end);
time_positive(1) = 0;
r_signal = r_signal(:,idx:end);

% resample
downsample_factor = T/T_sample;
[p,q] = rat(downsample_factor);
time_salmple = 0:T_sample:time_record_end;
t_signal_sample = source.get_signal(time_salmple);
r_signal_sample = r_signal(mod(time_positive, T_sample) == 0);

%% find dist
% true dist
true_dist_travel = s2t_dist + t2s_dist;
% aprrox dist
dist = Find.dist_travel(t_signal_sample, r_signal_sample, c, T_sample, effective_fft_size)/2
true_dist = true_dist_travel(1,1)/2
error = abs(true_dist - dist)

[signal_2, dist_2] = source.to_locs1_to_locs2(target.Locs, source.Locs, time_salmple, args_dict);
dist_2 = dist_2/2
signal_2 = squeeze(signal_2).';
dist_2_est = Find.dist_travel(t_signal_sample, signal_2, c, T_sample, effective_fft_size)/2
error_2 = abs(dist_2_est-dist_2)

%% plot signals and mixed
% figure
% subplot(2,1,1);
% stft(source.get_signal(time))
% subplot(2,1,2);
% stft(squeeze(signal(1,1,:)).')
% figure
% stft(source.get_signal(time).*conj(squeeze(signal(1,1,:)).'))

%% test range estimation
range = 0.001:0.001:2;
target_locs = [zeros(size(range)); zeros(size(range)); range].';
test_dist_changes(time_salmple, source, target_locs, effective_fft_size);

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
