% %%%%%%%%%
% loc          - is the cordinant of the source, size(3)
% func_pointer - is the pointer for the signal the the source emit
% args_func    - are the argumant the func use, exape time
% %%%%%%%%%
classdef Element
  properties
      Locs
      func_pointer
      args_func
      
      SavedSignals
  end
  methods
    % constractor
    function obj = Element(Locs, func_pointer, args_func)
      obj.Locs = Locs;
      obj.func_pointer = func_pointer;
      obj.args_func = args_func;

      obj.SavedSignals;
    end
    
    %%%% get signal at the source 
    function Signal = get_signal(obj, t)
      Signal = obj.func_pointer(t, obj.args_func{:});
    end
    
    function [Signals, Dist] = to_locs(obj, Locs, t, TargetsProperties)
      C = physconst('LightSpeed');
      n = size(obj.Locs, 1);
      m = size(Locs, 1);
      T = t(2)-t(1);
      
      Dist = calc_distance_matrix(obj.Locs, Locs);   % from this to Loc, size(n,m, 3)
      Pr = (TargetsProperties./(4*pi.*Dist.^2));
      Shift = round(Dist ./ (C.*T));
      % tmp_shift = Dist ./ C;
      
      Signals = zeros(n, m, length(t));
      for index = 1:n*m
        [ii,jj] = ind2sub([n, m], index);
        if size(obj.SavedSignals, 1) == 1
          Signals(ii,jj,:) = Pr(ii, jj) * circshift(obj.SavedSignals, [0, Shift(ii,jj)]);
          % tmp = Pr(ii, jj) * obj.func_pointer(t-tmp_shift, obj.args_func{:} );
          % Signals(ii,jj,:) = Pr(ii, jj) * interp1(t, obj.SavedSignals, Dist ./ (C.*T), 'nearest');
        else
          Signals(ii,jj,:) = Pr(ii, jj) * circshift(obj.SavedSignals(ii,:), [0, Shift(ii,jj)]);
        end
      end
    end

    %%%% get signal that hit loc1, and return to loc2. from source
    % locs1 - a cordinant vector size(n,3)
    % locs2 - a cordinant vector size(m,3)
    % args_dict -
    %   Pt - power of the source signal
    %   Gt - is the radiation intensity of the antenna in a given direction
    %        over that of an isotropic (Gt=4*pi*A/lambda^2) 
    %   sigma - radar cross section units [meters^2]
    %   Ae - effective area of receiving antenna
    function [signal, dist_travel] = to_locs1_to_locs2(obj, locs1, locs2, t, args_dict)
      C = physconst('LightSpeed');
      n = size(locs1, 1);
      m = size(locs2, 1);

      dist1 = sum( (obj.Locs - locs1).^2 , 2).^0.5;   % from source to loc1, size(1,n)
      dist2 = reshape(...
               sqrt(sum(bsxfun(@minus, locs1, permute(locs2, [3 2 1])).^2, 2)), ...
               n, m);      % from loc1 to loc2,   size(n,m)

      dist_travel = dist1 + dist2;
      Pr = (args_dict('Gt').*args_dict('Gt')./(4*pi.*dist1.^2)) .* ...
              (args_dict('sigma').*args_dict('Ae')./(4*pi.*dist2.^2));
      shift = dist_travel ./ C;

      signal = zeros(n, m, length(t));
      for index = 1:n*m
        [ii,jj] = ind2sub([n, m], index);
        signal(ii,jj,:) = Pr(ii, jj) * obj.func_pointer(t-shift(ii,jj), obj.args_func{:} );
      end

      return
    end
    %%%% plot dot in a 3D space
    function plot(obj, color, marker)
      plot3(obj.loc(1), obj.loc(2), obj.loc(3), 'LineWidth', 2, 'Color', color, 'Marker', marker);
    end
  end
end

function dist_matrix = calc_distance_matrix(locs1, locs2)
  % Calculate Euclidean distance between every pair of points in loc1 and loc2

  % Get number of points in each vector
  n = size(locs1, 1);
  m = size(locs2, 1);

  % Repeat each row of loc1 and loc2 m and n times, respectively, to create
  % matrices of size (n,m,3) and (n,m,3) that contain all pairs of points
  loc1_mat = repmat(reshape(locs1, n, 1, 3), [1, m, 1]);
  loc2_mat = repmat(reshape(locs2, 1, m, 3), [n, 1, 1]);

  % Calculate Euclidean distance between each pair of points using the norm function
  diff_mat = loc1_mat - loc2_mat;
  dist_matrix = sqrt(sum(diff_mat .^ 2, 3));
end


