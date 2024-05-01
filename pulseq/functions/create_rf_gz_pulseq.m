function [rf0,gz,gzReph,params] =  create_rf_gz_pulseq(params,tr_tmp)

    % RF pulses and FA depending on sequence type (VASO,ABC) and field strength,
    % Values for: Duration, apodization and timeBwProduct, can be modified to
    % avoid SAR or other constraints
    lims = params.gen.lims;
    for i=1:params.gen.n(3)
        % Array for flip angles
        if params.gen.vfa
            if i > 1
                % Variable flip angle
                e1 = exp(-tr_tmp/params.gen.ernst_t1);
                params.gen.fa(i) = asin((sin(params.gen.fa(1))*tan(params.gen.fa(i-1))) ...
                    /((e1*sin(params.gen.fa(1)))+((1-e1)*tan(params.gen.fa(i-1)))));

                % Cut-off at prescribed value
                if params.gen.fa(i) > params.gen.vfa_cutoff /180*pi
                    params.gen.fa(i) = params.gen.vfa_cutoff /180*pi;
                end
            end
        else
            % Fixed flip angle
            params.gen.fa(i) = params.gen.fa(1);
        end
        
        % Cut-off at 90 if >90
        if ~isreal(params.gen.fa(i))
            params.gen.fa(i) = 90/180*pi;
        end
        
        % Array of RF pulses, depending on sequence type:
        if params.gen.seq == 2
            % ABC
            [rf0(i), gz(i)] = mr.makeSincPulse(params.gen.fa(i),'system',lims,'Duration',1e-3,...
                'SliceThickness',params.gen.fov(3),'apodization',0.5);
        else
            % VASO/BOLD/Fieldmap
            if params.gen.field_strength == 7 || params.gen.field_strength == 7i
                % 7T and 7Ti
                [rf0(i), gz(i)] = mr.makeSincPulse(params.gen.fa(i),'system',lims,'Duration',2.56e-3,...
                    'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
            elseif params.gen.field_strength == 9
                % 9.4T
                [rf0(i), gz(i)] = mr.makeSincPulse(params.gen.fa(i),'system',lims,'Duration',2.56e-3,...
                    'SliceThickness',params.gen.fov(3),'apodization',0.5,'timeBwProduct',25);
            end
        end
        gzReph(i) = mr.makeTrapezoid('z',lims,'Area',-gz(i).area/2);
    end

end