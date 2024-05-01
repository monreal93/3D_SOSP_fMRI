function params = prepare_flip_angle(tr_tmp,gx, params)

    % Flip angle 
    % tr_tmp Rough estimate of TR, to calculate Ernst Angle
    if params.gen.ro_type == 'c'
        params.epi.n_lines = ceil(params.gen.n(2)/params.epi.ry)-((round(params.gen.n(2)/params.epi.ry)-round(params.gen.n(2)/params.epi.ry*params.epi.pf)));
        tr_tmp = (mr.calcDuration(gx)*(params.epi.n_lines+3))+mr.calcDuration(gx)+params.gen.te; % +3 of navigators
        tr_tmp = tr_tmp+3.5e-3; % 3.5e-3 (approx rf + rephasing time)
    end

    if params.gen.fa == 0
        % Ernst angle
        params.gen.fa(1) = acos(exp(-(tr_tmp)/(params.gen.ernst_t1)));
    else
        % Custom angle
        params.gen.fa(1) = params.gen.fa*pi/180;
    end


end