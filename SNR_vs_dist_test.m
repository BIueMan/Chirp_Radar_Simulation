clear
close all

addpath('objects');
addpath('func');
addpath('test_func');

%% set data
C = physconst('LightSpeed');
%% signal params
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time_transmit_start = 0; %-1e-5;
time_record_end = 7e-4;
T = 7.031256e-6;
%T = 1e-11;
time = time_transmit_start: T: time_record_end;

f_0 = 76e9;
f_1 = f_0 + 4e9;
amp = 1;
chirp_len = 125e-6;
c = (f_1-f_0)/chirp_len;
effective_fft_size = 1024*128;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% create noise
SNR = 200:-0.1:0;