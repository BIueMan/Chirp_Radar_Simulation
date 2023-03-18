% %%%%%%%%%%%%%%%%
%
% func - make sure func is diffined good at t<0. in our calculation we shift the time to t<0,
%        so if you want the signal to start from t=0. diffine it with "u(t)"
% %%%%%%%%%%%%%%%%
classdef Element_old
    %%%%
    % obj.loc - a matrix of size(m,3)
    % obj.saved_signal - matrix of sizze(m,len(t))
    % when:
    %  	m - is the elements locations
    %   len(t) - the lenth of the time vector.
    %%%%
    properties
        loc
        saved_signal
        func_pointer
        args_func
    end
    methods
      % %%%%%%%%%
      % loc          - is the cordinant of the source, size(3)
      % func_pointer - is the pointer for the signal the the source emit
      % args_func    - are the argumant the func use, exape time
      % %%%%%%%%%
        function obj = Element(loc, func_pointer, args_func)
            obj.loc = loc;
            obj.func_pointer = func_pointer;
            obj.args_func = args_func;
        end
        
        function h = h_to_loc(obj, t, target_locs)
          C = 3e8;
          h = zeros(size(obj.loc), size(target_locs,1), length(t));
          % for every targel loc, we will find h
          for ii = 1:size(obj.loc,1)
            dist = sum( (obj.loc(ii,:) - target_locs).^2, 2 ).^0.5;

            power = (4 * pi * dist.^2).^-1;

            % create time shift, and remove any time_shift that after t(end)
            shift = round( dist / (C * (t(2)-t(1))) ) + 1;
            shift(shift > length(t)) = 0;

            h(ii,:) = accumarray(shift, power, [length(t), 1]);
          end

          return
        end
        
        function signal = to_loc(obj, targets_loc, t)
          % check if saved signal is good
          if sum(size(obj.saved_signal) ~= [size(obj.loc, 1), length(t)])
             error('this.saved_signal must be at size [size(source_loc, 1), length(t)]')
          end
          
          % h is the delta shift for every source*target combination
          h = h_to_loc(obj, t, targets_loc);
          % shift every combination
          signal = zeros(size(h));
          for ii = 1:size(h,1)
              h_ii = reshape(h(ii, :, :), [size(h, 2), size(h, 3)]);
              signal(ii,:,:) = reshape(obj.saved_signal(ii,:) * h_ii, [1, size(h, 2), size(h, 3)]);
          end

          return
        end
        
        function signal = get_signal(obj, t)
          signal = obj.func_pointer(t, obj.args_func{:} );
        end
        
        % locs1 - a cordinant vector size(3,n)
        % locs2 - a cordinant vector size(3,m)
        function [signal, dist_travel] = to_locs1_to_loc2(obj, locs1, locs2, t)
          C = 3e8;
          
          dist1 = sum( (obj.loc - locs1).^2 , 2).^0.5;   % a size(1,n)
          dist2 = reshape(...
                   sqrt(sum(bsxfun(@minus, locs1, permute(locs2, [3 2 1])).^2, 2)), ...
                   size(locs1, 1), size(locs2, 1));

          
          dist_travel = dist1 + dist2;
          power = (4 * pi * dist_travel.^2).^-1;
          shift = dist_travel ./ C;
          
          signal = zeros(size(locs1,1), size(locs2,1), length(t));
          for ii = 1:size(signal,1)
            for jj = 1:size(signal,2)
              signal(ii,jj,:) = power(ii, jj) * obj.func_pointer(t-shift(ii,jj), obj.args_func{:} );
            end
          end
          
          return
        end

        function plot_dot(obj, color, marker)
          for i = 1:size(obj.loc, 2)
              plot3(obj.loc(1, i), obj.loc(2, i), obj.loc(3, i), 'LineWidth', 2, 'Color', color, 'Marker', marker);
          end
        end
      end
end