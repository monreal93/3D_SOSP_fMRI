function [B1,phase,rf_complex,Gz,fa] = tr_foci(params)
    % Some parameters
    rf_ampl = 2.1865;
    larmor_freq = 4.2576e+04;
    rf_dur = params.vaso.foci_dur;

% RF Slab-Selective
    Gs=1;
    Amax = 3.32;
    B1max =  0.327152;
    w = 0.30; r1=0.64; r2=0.27; r3=0.59; r4=0.00; r5=1.00;
    % AMM: For ASL, use 300, for VASO 150
%     mu=7.71*(300/100); % 7.71* BW in percentage
    mu = 7.71*(params.vaso.foci_bw/100);
    beta=3.90; t1=0.25; t2=0.40;

    c=Amax*(1-r1);
    Amin=r2*c;
    b1=r3*(c-Amin)/((w-1)^2);
    b2=r4*(1-r3)*(c-Amin)/((w-1)^4);
    b3=r5*(1-r4)*(1-r3)*(c-Amin)/((w-1)^6);
    b4=(1-r5)*(1-r4)*(1-r3)*(c-Amin)/((w-1)^8);
    
    % Timing for RF pulse
    dt = 2/((rf_dur*1e6)-1);
    t=-rf_dur/rf_dur:dt:(rf_dur/rf_dur);

    % Timing for gradient
    grad_raster = 10;
    grad_ramp = 400;
    flat_samples = (rf_dur./grad_raster.*1e6);
    ramp_samples = grad_ramp./grad_raster;
    grad_samples = round((ramp_samples*2)+flat_samples);
    t_grad = (-rf_dur/rf_dur)+(1/grad_samples/4):1/(grad_samples/2):(rf_dur/rf_dur)-+(1/grad_samples/4);
    
    
    T= @(t) (t1*(t.^5)+t2*(t.^3)+t)/(t1+t2+1);
    x=beta*T(t);
    sech=2./(exp(x)+exp(-x));
    mytanh= @(x) ((exp(x)-exp(-x))./(exp(x)+exp(-x)));
    
    A1 = @(t) Amax*(1-((r1*(t+1))/w));
    A2 = @(x) Amin+(b1*(x.^2))+(b2*(x.^4))+(b3*(x.^6))+(b4*(x.^8));
    grad1 = @(t)(A1(t)*Gs/A1(-1));
    freq1 = @(t) (-A1(t).*beta.*mytanh(beta*T(t)));
    grad2 = @(t) (A2(t)*Gs/A1(-1));
    freq2 = @(t) (-A2(t).*beta.*mytanh(beta*T(t)));
    
    dummy = 0;
   
    % Checking the timing for RF segments
    tt1=t(t<=(w-1));
    tt2=t(t>(w-1)); tt2=tt2(tt2<(1-w)); % flat area of pulse
    tt3=t(t>(w-1)); tt3=tt3(tt3>(1-w));

    % Checking the timing for Gradient segments  
    t_grad_flat = -tt2(1)/tt2(1):1/(flat_samples-1)*2:tt2(end)/tt2(end);
    tt1_grad=t_grad_flat(t_grad_flat<=(w-1));
    tt2_grad=t_grad_flat(t_grad_flat>(w-1)); tt2_grad=tt2_grad(tt2_grad<(1-w)); % flat area of pulse
%     tt3_grad=t(t_grad_flat>(w-1)); tt3_grad=tt3_grad(tt3_grad>(1-w));
    
    % Segment 1
    ampl_A1=A1(tt1);
    grad_A1 = grad1(tt1);
    freq_A1 = freq1(tt1);
    gz_1 = grad1(tt1_grad);
    
    for i=1:length(tt1)
        dummy = dummy-((grad_A1(i)+(mu.*freq_A1(i)))./length(t));
%         dummy = dummy-(grad_A1(i)+(freq_A1(i)));
%         phase1(i) = dummy- renzo_dummyphase; 
        phase1(i) = dummy; % AMM: Removing dummyphase..
    end
    
    % Segment 2
    ampl_A2=A2(tt2);
    grad2_A1 = grad2(tt2);
    freq2_A1 = freq2(tt2);
    gz_2 = grad2(tt2_grad);
    
    for i=1:length(tt2)
        dummy = dummy-((grad2_A1(i)+(mu.*freq2_A1(i)))./length(t));
%         phase2(i) = dummy- renzo_dummyphase;
        phase2(i) = dummy; % AMM: Removing dummyphase..
   end
    
    % Segment 3
    ampl_A3 = flip(ampl_A1);
    grad3_A1 = grad2(tt3);
    freq3_A1 = freq2(tt3);
    gz_3 = flip(gz_1);
    
    % Gz Rampup
    tmp = 0:ramp_samples-1;
    gz_ramp = tmp./ramp_samples;
    
    for i=1:length(tt3)
        dummy = dummy-((grad3_A1(i)+(mu.*freq3_A1(i)))./length(t));
%         phase3(i) = dummy- renzo_dummyphase;
        phase3(i) = dummy; % AMM: Removing dummyphase..
    end
    
    A = [ampl_A1 ampl_A2 ampl_A3];
    F = [phase1 phase2 phase3];
    Gz = [gz_ramp gz_1 gz_2 gz_3 flip(gz_ramp)];
    Gz = (Gz.*rf_ampl.*larmor_freq).';
    F_wrap = wrapTo2Pi(-F);
    
    B1=A.*sech./B1max;
%     B1 = padarray(B1,[0 1],'both');  % AMM: trying to pad the zero...
%     freq=-A.*mu.*beta.*mytanh;

    grad = A*Gs;
%     freq=-A.*beta.*mytanh(x);  % Without mu
   phase = F_wrap;
%    phase = padarray(phase,[0 1],'both');  % AMM: trying to pad the zero...
   rf_complex = B1.*exp(1i.*phase);
   
   % Here I want to calculate the FA taking into account the amplitude
   m_dBW_kHz = 1.0E-3 * mu * beta* Amax / pi()* 1/(1.0E-6 *length(t));
   ampl = 2*pi*m_dBW_kHz/(267522209. * params.vaso.foci_ampl);
   fa = (2*pi)/(ampl*360*(larmor_freq*1e3*2*pi)*trapz(0:1e-6:params.vaso.foci_dur-1e-6,abs(rf_complex)));
   fa = fa*180/pi;   % FA in degrees
   
end