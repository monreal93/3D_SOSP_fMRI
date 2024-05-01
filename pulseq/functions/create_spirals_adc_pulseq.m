function [gx,gy,gx_pre,gy_pre,adc,params] =  create_spirals_adc_pulseq(spiral_grad_shape,adcSamples,adcDwell,params)
        
        lims = params.gen.lims;

        % Dummy trapezoid for pre-allocation
        gx_pre = mr.makeTrapezoid('x',lims,'maxGrad',25e-3*lims.gamma,'Area',0);
        gy_pre = mr.makeTrapezoid('x',lims,'maxGrad',25e-3*lims.gamma,'Area',0);
        
        for j=1:params.gen.n(3)
            for i=1:params.spi.interl

                % Readout gradients
                % Temp: for Pulseq version 1.4.1
                gx(j,i) = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:,i,j), lims, 'first', 0, 'last', 0);
                gy(j,i) = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:,i,j), lims, 'first', 0, 'last', 0);
                
%                 % For Pulseq version 1.4.0
%                 gx(j,i) = mr.makeArbitraryGrad('x',spiral_grad_shape(1,:,i,j), lims);
%                 gy(j,i) = mr.makeArbitraryGrad('y',spiral_grad_shape(2,:,i,j), lims);

                % ADC 
                adc = mr.makeAdc(adcSamples,lims,'Dwell',adcDwell);

% %                 Pulseq (and siemens) define the samples to happen in the center of the dwell period
%                 time_to_center = adc.dwell*((adcSamples-1)/2+0.5);
%                 adc.delay = round((gx.riseTime+gx.flatTime/2-time_to_center)/lims.rfRasterTime)*lims.rfRasterTime;

                % Getting the correct number to split the ADC
                % The ADC obj has to be splitted into N equal parts, with duration multiple of 10us
                adc_total_samples = adcSamples;
                for k = 1:50
                    if mod(adc_total_samples,k) == 0 && mod(adc_total_samples/k*adcDwell,10e-9) == 0 && (adc_total_samples/k) < 8192 && mod(adc_total_samples/k,4) == 0
                        adcSplit = adc_total_samples/k;
                        break
                    end
                end
                params.gen.adc_split = adcSplit;
                params.gen.adc_segments = k;

                if params.spi.type == 0
                    if params.gen.echos > 1
                        tmp = cumsum(gx(j,i).waveform);
                        tmp = tmp./lims.gamma*100;   
                        area = ((tmp(1)-tmp(end))*4)+(((tmp(1)-tmp(end))*4)*6.5/100);
                        gx_pre(j,i) = mr.makeTrapezoid('x',lims,'Area',area);
                        tmp = cumsum(gy(j,i).waveform);
                        tmp = tmp./lims.gamma*100;   
                        area = ((tmp(1)-tmp(end))*4)+(((tmp(1)-tmp(end))*4)*6.5/100);
                        gy_pre(j,i) = mr.makeTrapezoid('y',lims,'Area',area);
                    end
                elseif params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
                    % Pre-phasing gradients, I adjust for 6.5% trajectory
                    % error to be closer to zero at center of k-space
                    % for rxy=3=6.5, rxy=4=
                    err = 6.5; % (6.5)
                    tmp = cumsum(gx(j,i).waveform);
                    tmp = tmp./lims.gamma*100;  
                    if params.spi.type == 1
                        area = ((tmp(1)-tmp(end))*4)+(((tmp(1)-tmp(end))*4)*err/100);
                    elseif params.spi.type == 3 || params.spi.type == 4
                        area = ((tmp(1)-tmp(end/2))*4)+(((tmp(1)-tmp(end/2))*4)*err/100);
                    end
                    gx_pre(j,i) = mr.makeTrapezoid('x',lims,'maxGrad',25e-3*lims.gamma,'Area',area);
                    tmp = cumsum(gy(j,i).waveform);
                    tmp = tmp./lims.gamma*100;   
                    if params.spi.type == 1
                        area = ((tmp(1)-tmp(end))*4)+(((tmp(1)-tmp(end))*4)*err/100);
                    elseif params.spi.type == 3 || params.spi.type == 4
                        area = ((tmp(1)-tmp(end/2))*4)+(((tmp(1)-tmp(end/2))*4)*err/100);
                    end
                    gy_pre(j,i) = mr.makeTrapezoid('y',lims,'maxGrad',25e-3*lims.gamma,'Area',area);

                end
            end
        end

end