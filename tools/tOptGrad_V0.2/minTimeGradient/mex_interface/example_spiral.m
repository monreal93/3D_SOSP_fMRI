%% Spiral

disp('######################################');
disp('#### Design a dual density spiral ####');
disp('####                              ####');
disp('######################################');
disp(' ');

[k_rv,g_rv,s_rv,time_rv,Ck_rv] = vdSpiralDesign(1, 16, 0.83,[55,55,10,10],[0,0.2,0.3,1],4,15,4e-3,'cubic');
[k_riv,g_riv,s_riv,time_riv,Ck_riv] = vdSpiralDesign(0, 16, 0.83,[55,55,10,10],[0,0.2,0.3,1],4,15,4e-3,'cubic');

L = max(length(s_riv), length(s_rv));

figure, subplot(2,2,1), plot(k_rv(:,1), k_rv(:,2)); title('k-space'); axis([-6 6 -6 6]);
subplot(2,2,2), plot(g_riv(:,1)); axis([0,L,-4.5,4.5]); title('gradient waveforms (R. Invariant)')
hold on, subplot(2,2,2), plot(g_riv(:,2), 'r');
legend('gx', 'gy', 'Location', 'NorthEast');
subplot(2,2,3), plot((g_rv(:,1).^2 + g_rv(:,2).^2).^0.5, '--'), 
hold on, subplot(2,2,3), plot((g_riv(:,1).^2 + g_riv(:,2).^2).^0.5, 'r');  axis([0 L 0 6]);
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthEast'); title('gradient magnitude')
subplot(2,2,4), plot((s_rv(:,1).^2 + s_rv(:,2).^2).^0.5, '--'); title('slew-rate magnitude');  axis([0 L 0 20]);
hold on, subplot(2,2,4), plot((s_riv(:,1).^2 + s_riv(:,2).^2).^0.5, 'r'); 
legend('rotationally invariant', 'rotationally variant', 'Location', 'SouthWest');