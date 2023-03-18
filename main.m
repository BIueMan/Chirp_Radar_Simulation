clear
close all
addpath('objects');
addpath('func');
addpath('test_func');

%% data
global C
C = physconst('LightSpeed');
amp = 1;
error = 1;
num_of_reppet_per_test = 5;
aqurecy_permission = 4;

R = 0.1:0.1:200;
SNR = 10.^( (1:-0.1:-100)./20 );
R_vs_SNR_dB = zeros(length(SNR),1);

snr_ii = 1;
progress_bar = waitbar(0,'progress\_bar for R: ');
for r_ii = 1:length(R)
  r = R(r_ii);
  % continue until SNR distory the dist
  while(true)
    aqurecy = 0;
    snr = SNR(snr_ii);
    for jj = 1:num_of_reppet_per_test
      dist = test_dist_vs_SNR(r, snr);
      if (abs(dist-r) < error)
        aqurecy = aqurecy + 1;
      end
    end
    % if SNR is to big then brack
    if (aqurecy < aqurecy_permission || snr_ii >= length(SNR))
      R_vs_SNR_dB(r_ii) = db(snr);
      waitbar(r_ii/length(R), progress_bar);
      break
    else      
      snr_ii = snr_ii+1;
    end
  end
end
close(progress_bar);

figure
plot(R, R_vs_SNR_dB)
title('R vs SNR')
xlabel('R [m]')
ylabel('SNR [dB]')
