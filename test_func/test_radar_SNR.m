function [range_true, var_error] = test_radar_SNR(time, source, target_locs, effective_fft_size, SNR, repet, doplot)
  assert(isa(source, 'Element'), 'source must be an Element class');
  if nargin < 7
    doplot = true;
  end
  args_dict = containers.Map({'Pt', 'Gt', 'sigma', 'Ae'}, ...
                           { 1  ,  0.1,  0.01  ,  0.1});
  T = time(2) - time(1);
  f_1 = cell2mat(source.args_func(3));
  f_0 = cell2mat(source.args_func(1));
  chirp_len = cell2mat(source.args_func(2));
  c = (f_1 - f_0)/chirp_len;
  
  Error = zeros(repet, size(target_locs, 1));
  %% for every repet, add SNR to every target range
  for jj = 1:repet
    % for every target locs, get dist
    [signal, true_dist_travel] = source.to_locs1_to_locs2(target_locs, source.Locs, time, args_dict);
    signal = squeeze(signal);
    % add SNR
    noisey_signal = awgn(signal, SNR); % this line is broken. it dont add the SNR like we want

    t_signal = source.get_signal(time);
    r_signal = noisey_signal;

    range_est = arrayfun(@(ii) Find.dist_travel(t_signal, r_signal(ii,:), c, T, effective_fft_size)/2, 1:size(target_locs,1));

    range_true = (true_dist_travel/2).';
    Error(jj,:) = range_true - range_est;
  end
  
  %% extract varianse and mean
  mean_error = mean(Error);
  var_error = var(Error);
  
  %% plots
  if doplot    
    figure
    plot(range_true, range_est)
    title('one repet')
    ylabel('varianse')
    xlabel('range')
    
    figure
    plot(range_true, var_error)
    title(sprintf('varians - %d repet', repet));
    ylabel('varianse')
    xlabel('range')
    
    figure
    plot(range_true, mean_error)
    title(sprintf('mean - %d repet', repet));
    ylabel('mean')
    xlabel('range')
    
    % plot signals
    num_elements = 5;
    % compute the indices of the elements to extract
    plot_indices = round(linspace(1,size(signal,1),num_elements));
    figure
    for i = 1:length(plot_indices)
      subplot(length(plot_indices), 1, i);
      hold on
      plot(real(signal(plot_indices(i), :)));
      plot(real(noisey_signal(plot_indices(i), :)));
      title('range=', range_true(plot_indices(i)))
      hold off
    end
  end
end

