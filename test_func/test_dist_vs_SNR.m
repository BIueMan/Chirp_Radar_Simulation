function dist = test_dist_vs_SNR(r, snr, to_print)
  if nargin < 3
    to_print = false;
  end
  %% def emenets
  % def time
  time_transmit_start = -1e-5;
% time_record_start = 0;
  time_record_end = 2e-5;
  T = 1e-9;
  t = time_transmit_start: T: time_record_end;

  % def locs
  loc_source = [0, 0, 0];
  source = Element(loc_source);

  loc_target = [0, r, 0];
  target = Element(loc_target);

  %% def signal
  % def source 
  f_0 = 3e9;
  f_1 = 4e9;
  amp = 1;
  chirp_len = 5e-6;
  c = (f_1-f_0)/chirp_len;
  source_signal = amp * chirp(t,f_0, chirp_len, f_1, 'complex');
  source.saved_signal = source_signal;

  % cal returned signal
  target.saved_signal = source.to_loc(target.locs, t);
  return_signal = target.to_loc(source.locs, t);

  % remove the negative time
  [~ ,cut_of] = min(abs(t));
  t = t(cut_of:end);
  source_signal = source_signal(cut_of:end);
  return_signal = return_signal(cut_of:end);
  
  %% add noise
  true_dist = r;
  amp_return = amp * ((4*pi*true_dist.^2).^-1).^2;
  amp_noise = amp_return/snr;
  noise = amp_noise * (randn(size(return_signal)) + randn(size(return_signal))*1j);
  return_signal = return_signal + noise;


  %% find dist
  mixed = source_signal.*conj(return_signal);

  if(to_print)
    figure
    plot(t, real(mixed))
    title("time - real\_mixed")

    figure
    plot(db(abs(fft(mixed))))
    title("fft - mixed")
  end

  % calulate dist
  dist = Find.dist(source_signal, return_signal, c, T);
end

