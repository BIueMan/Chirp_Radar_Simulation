classdef Find
    methods     ( Static = true )
        % find the dist using chirp_source signal, and a returning one
        %   c = (f_1 - f_0)/chirp_len
        function dist_travel = dist_travel(source_signal, return_signal, c, time_step, effective_fft_size)
            C = physconst('LightSpeed');
            
            mixed_signal = source_signal.*conj(return_signal);
            fft_mixed = abs(fft(mixed_signal, effective_fft_size, 2));

            % use 3 point around the maxima to aproximat the true maxima
            function x_vertex = max_parable_aproximat(fft_mixed)
                %xf_max = find(fft_mixed == max(fft_mixed,[],2));
                [~, xf_max] = max(fft_mixed,[],2);
                x = [xf_max-1, xf_max, xf_max+1];
                x(x(:,1)<1, :) = x(x(:,1)<1, :)+1;
                x(x(:,3)>size(x, 1), :) = x(x(:,3)>size(x, 1), :)-1;
                y = reshape(fft_mixed(x).', [], 3);
                [x_vertex, ~] = parabolix_max_approx(x, y);
            end
            x_vertex = max_parable_aproximat(fft_mixed);

            Fs = 1/time_step;
            f_true = (x_vertex/length(fft_mixed)) * Fs;

            dist_travel = (f_true*C)/c; % it's the dist_trave, it's 2*R from source to target and beck
        end
        
        % asoming loc_0 is the location of both the source signal and the return_signal_0.
        function angle = angle_2d(loc_0, loc_1, return_signal_0, return_signal_1, c, time_step, effective_fft_size)
          h1 = Find.dist_travel(return_signal_0, return_signal_1, c, time_step, effective_fft_size);
          h2 = Find.dist_travel(return_signal_1, return_signal_0, c, time_step, effective_fft_size);
          d = norm(loc_0 - loc_1);
          if h1<h2
              if h1>d
                angle = pi;
              else
                angle = asin(h1/d) + pi/2;
              end
          else
              if h2>d
                angle = 0;
              else
                angle = acos(h2/d);% (pi-acos(h2/d)) + pi;
              end
          end
        end
        
        % return a function handler for the angle in 3D
        % circle_function(t), t = [0, 2pi]
        function circle_function = angle_function_3d_circle(loc_0, loc_1, source_signal, return_signal_0, return_signal_1, c, time_step, effective_fft_size)
            angle_2d = Find.angle_2d(loc_0, loc_1, source_signal, return_signal_0, return_signal_1, c, time_step, effective_fft_size);
            [theta, phi, ~] = cart2sph(loc_1 - loc_0);     % TODO: test if [theta, phi, ~] OR [phi, theta, ~]
            % return a function handler for the angle in 3D
            function circle_function = get_circle_3d(angle_2d, theta, phi, r)
                % get rotation matrixs
                rot_mat_z = [cos(phi) -sin(phi) 0;...
                             sin(phi)  cos(phi) 0;...
                             0         0        1];
                theta = -theta;
                rot_mat_y = [cos(theta)  0      sin(theta);...
                             0           1      0; ...
                             -sin(theta) 0      cos(theta)];
                % return a function to the 3D-circle
                circle_function = @(t) r * rot_mat_y * rot_mat_z *          ...
                                           [cos(angle_2d) * ones(size(t));  ...
                                            sin(angle_2d) * cos(t);         ...
                                            sin(angle_2d) * sin(t)];
            end
            circle_function = get_circle_3d(angle_2d, theta, phi, r);
        end
        
    end
end

