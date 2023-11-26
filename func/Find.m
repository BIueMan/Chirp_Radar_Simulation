classdef Find
    methods     ( Static = true )
        % find the dist using chirp_source signal, and a returning one
        %   c = (f_1 - f_0)/chirp_len
        function dist_travel = dist_travel(source_signal, return_signal, c, time_step, effective_fft_size, noise)
            persistent print_fft_count;
            if isempty(print_fft_count)
                print_fft_count = -1; % Initialize on the first call
            end
            print_fft_count = print_fft_count + 1;
            if nargin < 6
                noise = 0;
            end
            C = physconst('LightSpeed');
            
            mixed_signal = source_signal.*conj(return_signal) + noise;
            fft_mixed = abs(fft(mixed_signal, effective_fft_size, 2));

            % use 3 point around the maxima to aproximat the true maxima
            function x_vertex = max_parable_aproximat(fft_mixed) %% add errors
                %xf_max = find(fft_mixed == max(fft_mixed,[],2));
                xf_max = argmax(fft_mixed, 2);
                % make sure we dont fall on edge
                if (sum(xf_max == 1) || sum(xf_max == length(fft_mixed)))
                  warning('argmax is on the egde of the vector. could not use parabolic max approx')
                end

                x = [xf_max-1, xf_max, xf_max+1];
                y = fft_mixed(x);
                [x_vertex, ~] = parabolix_max_approx(x, y);
            end
%             x_vertex = max_parable_aproximat(fft_mixed);
            [~, x_vertex] = max(fft_mixed, [], 2);

            Fs = 1/time_step;
            f_true = (x_vertex/length(fft_mixed)) * Fs;

            dist_travel = (f_true*C)/c; % it's the dist_trave, it's 2*R from source to target and beck
            
            % plot
            if any(print_fft_count == [0 5]) % Changed this line             
              figure
              % Subplot 1: Plot the absolute value of mixed_signal
              subplot(2, 1, 1);
              plot(real(mixed_signal));
              title('Absolute Value of Mixed Signal');
              ylabel('Amplitude');

              % Subplot 2: Plot the shifted FFT of mixed_signal
              subplot(2, 1, 2);
              plot(fftshift(abs(fft_mixed)));
              title('Shifted FFT of Mixed Signal');
              ylabel('Magnitude');
            end
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

