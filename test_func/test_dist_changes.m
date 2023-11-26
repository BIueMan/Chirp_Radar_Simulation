function [range_true, range_est] = test_dist_changes(time, source, target_locs, effective_fft_size, doplot) % continuous
  assert(isa(source, 'Element'), 'source must be an Element class');
  if nargin < 5
    doplot = true;
  end
  args_dict = containers.Map({'Pt', 'Gt', 'sigma', 'Ae'}, ...
                           { 1  ,  0.1,  0.01  ,  0.1});
  T = time(2) - time(1);
  f_1 = cell2mat(source.args_func(3));
  f_0 = cell2mat(source.args_func(1));
  chirp_len = cell2mat(source.args_func(2));
  c = (f_1 - f_0)/chirp_len;
  
  %% for every target locs, get dist
  [signal, true_dist_travel] = source.to_locs1_to_locs2(target_locs, source.Locs, time, args_dict);

  t_signal = source.get_signal(time);
  r_signal = squeeze(signal);
  range_est = Find.dist_travel(t_signal, r_signal, c, T, effective_fft_size)/2;
% %   range_est = arrayfun(@(ii) Find.dist_travel(t_signal, r_signal(ii,:), c, T, effective_fft_size)/2, 1:size(target_locs,1));
%   for ii=1:size(target_locs,1)
%     % find dist
%     range_est(ii) = Find.dist_travel(t_signal, r_signal(ii,:), c, T, effective_fft_size)/2;
%   end
  range_true = true_dist_travel/2;
  error = abs(range_true - range_est);
  C = physconst('LightSpeed');
  range_max = C/(2*c*T);
  
  %% plots
  if doplot
    figure
    hold on
    plot(range_true, range_true, '-b', 'DisplayName', 'True Range');
    plot(range_true, range_est, '-r', 'DisplayName', 'Estimated Range');
    line([range_max, range_max], [0, range_true(end)], 'DisplayName', 'Calculated Max Estimated Range', 'LineStyle',  ':')
    hold off
    xlabel('True Range')
    ylabel('Range')
    legend show
    title('Comparison Estimated vs True Ranges')

    figure
    plot(range_true, error)
    xlabel('True Range')
    ylabel('Error')
    title('Error of Estimated Ranges')
  end
end

