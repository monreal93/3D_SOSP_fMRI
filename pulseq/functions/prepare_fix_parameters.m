function [params,ro_blocks] = prepare_fix_parameters(params)

    params.gen.del_k = (1./params.gen.fov);
    if params.gen.ro_type == 's'; params.gen.del_k(1:2) = params.gen.del_k(1:2).*params.spi.rxy; end
    if params.gen.ro_type == 'c'; params.spi.interl = 1; end
    if params.gen.ro_type == 'c'; params.gen.del_k(2) = params.gen.del_k(2).*params.epi.ry; end
    if params.mt.bold; params.mt.mt = 1; end        % If MT BOLD corrected, use MT 
    params.gen.del_k(3) = params.gen.del_k(3).*params.gen.kz;
    params.gen.n = round(params.gen.fov./params.gen.res);
    if and(params.gen.seq == 1,params.vaso.bold_ref == 1) || and(params.gen.seq == 2,params.mt.bold)
        ro_blocks = 2; 
    elseif params.gen.seq == 3
        ro_blocks = length(params.epi.te);
    else
        ro_blocks = 1; 
    end
    % Partition or phase oversampling
    params.gen.n_ov = params.gen.n;
    if params.gen.ph_oversampling > 0
        params.gen.del_k(3) = params.gen.del_k(3)/(1+(params.gen.ph_oversampling/100)); 
        params.gen.n(3) = round(round(params.gen.n(3)*(1+(params.gen.ph_oversampling/100)))/2)*2;
    end
    % Fat sat angle
    if params.gen.fs_angle == 0
        if params.gen.field_strength == 7  || params.gen.field_strength == 7i 
            params.gen.fs_angle = 80; % 110/80
        elseif params.gen.field_strength == 9
            params.gen.fs_angle = 40; % 110/80
        end
    end
    
    % Trying to make n multiple of 4,update res
    tmp = mod(params.gen.n,4);
    params.gen.n(1:2) = params.gen.n(1:2)+tmp(1:2);
    % ToDo: Check... Not sure if I do need this...
    params.gen.n(3) = params.gen.n(3)/params.gen.kz/params.gen.pf;
    % Making sure fovz/rz are integers and even
    if round((params.gen.fov(3)/params.gen.kz)*1000,9)/round((params.gen.fov(3)/params.gen.kz)*1000) ~= 1
        params.gen.fov(3) = round(params.gen.fov(3)*1000/2/params.gen.kz)*2*params.gen.kz*1e-3;
    end
    % if params.gen.seq == 3; params.gen.ro_type = 'c'; end   % if fieldmap, cartesian
    if params.gen.seq == 2; params.gen.ro_type = 's'; end   % if ABC, Spiral
    if params.mt.mt == 0; params.mt.bold = 0; end           % if no MT, no ref BOLD  
    if params.gen.ro_type == 'c'
        params.gen.seg = params.epi.seg;
    elseif params.gen.ro_type == 's'
        params.gen.seg = params.spi.interl;
    end

end