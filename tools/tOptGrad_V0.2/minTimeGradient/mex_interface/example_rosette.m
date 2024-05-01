disp('############################################');
disp('#### Design a rosette trajectory        ####');
disp('####                                    ####');
disp('############################################');
disp(' ');

Gmx = 4;
Smx = 15;
T = 17/Gmx;
Kmx = 6;
w1 = 0.147*2*pi*Gmx;
w2 = 0.087/1.02*2*pi*Gmx;
t = 0e-3:4e-3:T;
C = Kmx*sin(w1*t').*exp(i*w2*t');
C = [real(C) imag(C) t'];

[C_rv, time_rv, g_rv, s_rv, k_rv] = minTimeGradient(C,0);          % Rotationally variant solution
[C_riv, time_riv, g_riv, s_riv, k_riv] = minTimeGradient(C,1, 0);  % Rotationally invariant solution
L = max(length(s_riv), length(s_rv));

figure, subplot(2,2,1), plot(C_rv(:,1), C_rv(:,2)); title('k-space'); axis([-6 6 -6 6]);
subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Invariant)')
hold on, subplot(2,2,2), plot(g_riv(:,2), 'r');
legend('gx', 'gy', 'Location', 'NorthEast');
subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, '--'), 
hold on, subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');  axis([0 L 0 6]);
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2).^0.5, '--'); title('slew-rate magnitude');  axis([0 L 0 20]);
hold on, subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, 'r'); 
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthWest');