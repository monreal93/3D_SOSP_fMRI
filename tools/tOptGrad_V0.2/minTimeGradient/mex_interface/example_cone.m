disp('############################################');
disp('####   Design a cone trajectory         ####');
disp('####                                    ####');
disp('############################################');
disp(' ');

r = linspace(0,5, 500)';
th = linspace(0,2*pi, 500)';
C = r.*exp(2*pi*1i*th);
C = [real(C) imag(C) r];

[C_rv, time_rv, g_rv, s_rv, k_rv] = minTimeGradient(C,0);          % Rotationally variant solution
[C_riv, time_riv, g_riv, s_riv, k_riv] = minTimeGradient(C,1, 0);  % Rotationally invariant solution
L = max(length(s_riv), length(s_rv));

figure, subplot(2,2,1), plot3(C_rv(:,1), C_rv(:,2), C_rv(:,3)); title('k-space'); axis([-6 6 -6 6]);
subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Invariant)')
hold on, subplot(2,2,2), plot(g_riv(:,2), 'r');
hold on, subplot(2,2,2), plot(g_riv(:,3), 'g');
legend('gx', 'gy', 'gz', 'Location', 'NorthEast');
subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, '--'), 
hold on, subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');  axis([0 L 0 6]);
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2+ s_rv(:,3).^2).^0.5, '--'); title('slew-rate magnitude');  axis([0 L 0 20]);
hold on, subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2+ s_riv(:,3).^2).^0.5, 'r'); 
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthWest');