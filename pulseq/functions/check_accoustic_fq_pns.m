function check_accoustic_fq_pns(seq,params, grad_file)
    lims = params.gen.lims;
    % t     time of samples in microseconds, y     sample values, Fs    sampling rate in Hz
    Fs = 1/lims.adcRasterTime/100;
    gradients = seq.gradient_waveforms1()/lims.gamma*1000;   % Full sequence
    gradients = gradients(1:2,:);               % taking only gx and gy

    % Full sequence (all readouts)
    time = linspace(0,seq.duration(),size(gradients,2)).*1e6; % Full sequence
    
%     % Only 1 readout
%     gradients = gradients(:,1:floor(length(gradients)/params.gen.n_ov(3)));   % 1 readout
%     time = time(1:length(gradients));   % 1 readout

    gaxes = ['X' 'Y'];
    for i=1:length(gaxes)   
        [~, ~, ~, ~, bpass(i)] = gradFreqPlot_pulseq(time,gradients(i,:),Fs,gaxes(i),params.gen.field_strength);
    end
    
    % For safety reasons if test fails, the script will stop and
    % sequence won't be written....
    if ~bpass(1) || ~bpass(2)
        error(sprintf('########################################## \n Sequence failed resonance test,Please modify Spiral/EPI gradient Amplitde and Slew Rate and try again \n##########################################'))
    end

    if ~isempty(grad_file)
        seq.calcPNS(grad_file);
    end
    
end