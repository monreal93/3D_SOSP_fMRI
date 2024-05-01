function  [gx,gy_blips,gx_pre,gy_pre,gy_blip_up,gy_blip_down,adc,params] = create_epi_adc_pulseq(params)
    lims = params.gen.lims;

    
    gy_pre_area = -(params.gen.del_k(2)*params.gen.n(2)/2/params.epi.ry)*params.epi.pf;
    if params.epi.pf ~= 1; gy_pre_area = gy_pre_area+(gy_pre_area/4*-1); end
    gy_pre = mr.makeTrapezoid('y',lims,'maxGrad',15e-3*lims.gamma,'Area',gy_pre_area);
    gy_blip = mr.makeTrapezoid('y',lims,'Area',params.gen.del_k(2));

    tot_bw = params.gen.n(1)*params.epi.bw_px;
    adcDwell = floor((1./tot_bw/lims.adcRasterTime))*lims.adcRasterTime;

    % ADC
    adcSamples = params.gen.n(1);
    % In SIEMENS number of ADC samples should be divisible by 4
    adcSamples = floor(adcSamples/4)*4;  

    % Making sure we are aligned to adcRastertime
    params.epi.bw = 1./adcDwell;
    adc = mr.makeAdc(adcSamples,'Dwell',adcDwell,'Delay',mr.calcDuration(gy_blip)/2); % Ramp sampling

    % Getting the correct number to split the ADC
    % The ADC obj has to be splitted into N equal parts, with duration multiple of 10us
    adc_total_samples = adcSamples*round(params.gen.n(2)/params.epi.ry*(params.epi.pf))*params.gen.n(3);
    % if params.gen.dork; adc_total_samples = adc_total_samples + adc_post.numSamples; end
    for i = 1:500
        if mod(adc_total_samples,i) == 0 && mod(adc_total_samples/i*adcDwell,10e-9) == 0 && (adc_total_samples/i) < 8192 && mod(adc_total_samples/i,4) == 0
            adcSplit = adc_total_samples/i;
            break
        end
    end
    params.gen.adc_split = adcSplit;
    params.gen.adc_segments = i;
    
    tmp = floor((adc.duration+(adc.delay*2))./lims.gradRasterTime)*lims.gradRasterTime;
    gx = mr.makeTrapezoid('x',lims,'Area',params.gen.n(1)*params.gen.del_k(1),'Duration',tmp);
    gx_pre = mr.makeTrapezoid('x',lims,'maxGrad',30e-3*lims.gamma,'Area',-gx.area/2);
    gx.amplitude = -gx.amplitude;

    % Pulseq (and siemens) define the samples to happen in the center of the dwell period
%     time_to_center=adc.dwell*((adcSamples-1)/2+0.5);
    adc.delay = round(adc.delay/lims.rfRasterTime)*lims.rfRasterTime;

    % Let's split the gy_blip
    gy_blip_up = mr.makeExtendedTrapezoid('y',lims,'times',[0 gy_blip.riseTime],'amplitude',[0 gy_blip.amplitude]);
    gy_blip_down = mr.makeExtendedTrapezoid('y',lims,'times',[0 gy_blip.riseTime],'amplitude',[gy_blip.amplitude 0]);
    [gy_blip_up,gy_blip_down,~] = mr.align('right',gy_blip_up,'left',gy_blip_down,gx);
    gy_blips = mr.addGradients({gy_blip_down, gy_blip_up},lims);

end