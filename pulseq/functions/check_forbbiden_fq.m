
function f_pwr = check_forbbiden_fq(g_hz,scanner,plot_graph)
% g_hz = spiral.grad_x.seg(1).shot(1).kslice(1).grad(1).waveform;
% g_hz = gx(1).waveform;

    g = g_hz/42.58e3;
    
    if scanner == 7
        f1 = 1100;
        bw1 = 300;
        f2 = 550;
        bw2 = 100;
    elseif scanner == 9
        f1 = 3000;
        bw1 = 100;
    elseif scanner == 7i 
        f1 = 575;
        bw1 = 70;
        f2 = 1060;
        bw2 = 70;
        f3 = 1410;
        bw3 = 70;
        f4 = 1920;
        bw4 = 120;
    end
    oversample = 4;
    
    L=length(g)*oversample;
    L=L-rem(L,2);
    dt=(10e-6);
    ft=fftshift(abs(fft(g,L)));
    ft=ft(L/2+1:end);
    f=linspace(0,1/dt/2,L/2);

    % Computing the total power
    ft_tmp = interp1(linspace(1,max(f),length(ft)),ft,linspace(1,max(f),max(f)));
    pwr = trapz(ft_tmp);

    % Power in forbidden fq...
    f1_pwr = trapz(ft_tmp(f1-bw1/2:f1+bw1/2));
    % Perc in forbidden fq
    f_pwr = (f1_pwr)/pwr;
    if scanner == 7
        f2_pwr = trapz(ft_tmp(f2-bw2/2:f2+bw2/2));
            % Perc in forbidden fq
            f_pwr = (f1_pwr+f2_pwr)/pwr;
    end
    
    if plot_graph
        figure,hold on
        plot([f1-bw1/2 f1-bw1/2],[0,10],'.-r')
        plot([f1+bw1/2 f1+bw1/2],[0,10],'.-r')
        if scanner == 7 || scanner == 7i 
            plot([f2-bw2/2 f2-bw2/2],[0,10],'.-r')
            plot([f2+bw2/2 f2+bw2/2],[0,10],'.-r')
        end
        if scanner == 7i
            plot([f3-bw3/2 f3-bw3/2],[0,10],'.-r')
            plot([f3+bw3/2 f3+bw3/2],[0,10],'.-r')
            plot([f4-bw4/2 f4-bw4/2],[0,10],'.-r')
            plot([f4+bw4/2 f4+bw4/2],[0,10],'.-r')
        end

        axis([0 4000 0 1]);
        plot(f,(ft)/max(ft),'k','LineWidth',1)
    end

end