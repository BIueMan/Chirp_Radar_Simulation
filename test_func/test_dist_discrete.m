function [range_true, range_est] = test_dist_discrete(time, T_sample, source, target, effective_fft_size, doplot)
    assert(isa(source, 'Element'), 'source must be an Element class');
    if nargin < 6
      doplot = true;
    end
    
    target_properties = 1;

    T = time(2) - time(1);
    f_1 = cell2mat(source.args_func(3));
    f_0 = cell2mat(source.args_func(1));
    chirp_len = cell2mat(source.args_func(2));
    c = (f_1 - f_0)/chirp_len;
    
    %% for every target locs, get dist
    %%% from source to target
    [signal, dist] = source.to_locs(target.Locs, time, target_properties);
    target.SavedSignals = squeeze(signal(1,:,:));
    clear signal
    s2t_dist = dist(1,:,:);
    %%% from target to source
    [signal, dist] = target.to_locs(source.Locs, time, target_properties);
    clear target.SavedSignals
    r_signal = squeeze(signal);
    t2s_dist = squeeze(dist);
    clear signal

    %%% remove negative time
    [~,idx] = min(abs(time));
    time_positive = time(idx:end);
    time_positive(1) = 0;
    r_signal = r_signal(:,idx:end);

    %%% resample
    time_salmple = 0:T_sample:time(end);
    t_signal_sample = source.get_signal(time_salmple);
    r_signal_sample = r_signal(:, mod(time_positive, T_sample) == 0);

    %% find dist
    %%% true dist
    true_dist_travel = s2t_dist + t2s_dist.';
    %%% aprrox dist
    range_est = (Find.dist_travel(t_signal_sample, r_signal_sample, c, T_sample, effective_fft_size)/2).'
    range_true = true_dist_travel/2
    error = abs(range_true - range_est)

    %% calculate max range
    C = physconst('LightSpeed');
    range_max = C/(2*c*T_sample);
    
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
  
  