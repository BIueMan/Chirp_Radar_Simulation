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
T_sample = 7e-6;
T = T_sample;
time = [time_transmit_start: T: time_record_end];
[~, idx] = min(abs(time));
time = time - (time(idx) < 0) * time(idx + 1) - (time(idx) >= 0) * time(idx);

f_0 = 76e9;
f_1 = f_0 + 4e9;
amp = 1;
chirp_len = 125e-6;
c = (f_1-f_0)/chirp_len;
effective_fft_size = 1024*128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set source and target
source_locs = [0,0,0;
               0,0,0.01];
source = Element(source_locs(1,:), @chirp, {f_0, chirp_len, f_1, 'complex'});
source.SavedSignals = source.get_signal(time);
args_dict = containers.Map({'Pt', 'Gt', 'sigma', 'Ae'}, ...
                           { 1  ,    1,  1     ,  1  });
target_properties = 1;

% crate set of angles
true_angles = double(0:0.1:1) * pi;
r = 0.1;
target_locs = r * [zeros(1,length(true_angles));sin(true_angles);cos(true_angles)].';
target = Element(target_locs, NaN, NaN);

% plot points
scatter3(target_locs(:,1), target_locs(:,2), target_locs(:,3))
hold on
scatter3(source_locs(:,1), source_locs(:,2), source_locs(:,3))
hold off



%% get signal
% [r_signal, true_dist_travel] = source.to_locs(target.Locs, time, target_properties);
[r_signal, true_dist_travel] = source.to_locs1_to_locs2(target.Locs, source_locs, time, args_dict);
s_signal = source.get_signal(time);

%% find angle
est_angles = zeros(size(true_angles));
for i = 1:size(r_signal, 1)
    r_signal_i = reshape(r_signal(i,:,:), [size(r_signal, 2), size(r_signal, 3)]);
    find = Find();
    est_angles(i) = find.angle_2d(source_locs(1,:), source_locs(2,:), r_signal_i(1,:), r_signal_i(2,:), c, T, effective_fft_size);
end

angles_true_est = [true_angles; est_angles]
