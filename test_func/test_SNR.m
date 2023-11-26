function [range_true, var_error] = test_SNR(time, source, target_loc, effective_fft_size, SNRs, repet, doplot)
  assert(isa(source, 'Element'), 'source must be an Element class');
  if nargin < 7
    doplot = true;
  end
  args_dict = containers.Map({'Pt', 'Gt', 'sigma', 'Ae'}, ...
                           { 1  ,  1,  1  ,  1});
  T = time(2) - time(1);
  f_1 = cell2mat(source.args_func(3));
  f_0 = cell2mat(source.args_func(1));
  chirp_len = cell2mat(source.args_func(2));
  c = (f_1 - f_0)/chirp_len;
  
  Error = zeros(size(target_loc, 1), repet, length(SNRs));
  % first get signal
  [signal, true_dist_travel] = source.to_locs1_to_locs2(target_loc, source.Locs, time, args_dict);
  r_signal = reshape(signal, [size(signal,1), size(signal,3)]);
  range_true = (true_dist_travel/2).';
  
  t_signal = source.get_signal(time);
  A = max(abs(t_signal))*max(abs(r_signal.')).';
  N = size(r_signal, 2);
  
  %% for every SNR, repet. add SNR and test it
  for snr_ii = 1:length(SNRs)
      SNR = SNRs(snr_ii)
      sigma = (A.^2./SNR).^0.5;
      %% for every repet
      for rep_ii = 1:repet
        % add noise
        noise_real = sigma .* randn(size(r_signal));
        noise_imag = sigma .* randn(size(r_signal));
        %r_signal = r_signal + complex(noise_real, noise_imag);

        range_est = Find.dist_travel(t_signal, r_signal, c, T, effective_fft_size, complex(noise_real, noise_imag))/2;

        range_true = (true_dist_travel/2);
        Error(:, rep_ii, snr_ii) = range_true - range_est;
      end
  end
  
  %% extract varianse and mean
  mean_error = squeeze(mean(Error, 2));
  var_error = squeeze(std(Error, 0, 2));
  if size(mean_error,2)==1 % Check if the result is a matrix
    mean_error = mean_error.'; % Transpose the matrix to get a column vector
    var_error = var_error.';
  end
  
  %% plots
  if doplot
    CRLB_cal = @(SNR, B, N) (sqrt((3*physconst('LightSpeed')^2) ./ ((2*pi)^2 * N*B^2 .* SNR)));
    CRLB = CRLB_cal(SNRs, (f_1-f_0), N);
    
    figure
%     subplot(2, 1, 1);
%     plot(SNRs, CRLB, '-o', 'DisplayName', 'CRLB')
%     hold on
%     for ii = 1:length(range_true)
%       plot(SNRs, var_error(ii, :), '-*', 'DisplayName', sprintf('range = %d', range_true(ii)))
%     end
%     hold off
%     title(sprintf('varians - %d repet', repet));
%     ylabel('STD')
%     xlabel('SNR')
%     legend('show', 'Location', 'best')
%     
%     subplot(2, 1, 2);
    semilogy(10*log10(SNRs), 1e3*CRLB, '-o', 'DisplayName', 'CRLB')
    hold on
    for ii = 1:length(range_true)
      semilogy(10*log10(SNRs), 1e3*var_error(ii, :), '-*', 'DisplayName', sprintf('range = %d cm', 1e2*range_true(ii)))
    end
    hold off
    grid on
    title(sprintf('loglog varians - %d repet', repet));
    ylabel('STD[mm]')
    xlabel('SNR[db]')
    legend('show', 'Location', 'best')
    
%     figure
%     hold on
%     for ii = 1:length(range_true)
%       plot(SNRs, mean_error(ii,:), 'DisplayName', sprintf('range = %d', range_true(ii)))
%     end
%     hold off
%     title(sprintf('mean - %d repet', repet));
%     ylabel('mean')
%     xlabel('range')
%     legend('show', 'Location', 'best')
  end
end

