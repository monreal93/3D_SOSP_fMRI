function [rf_fs,gx_fs,gy_fs,gz_fs] = prepare_fat_sat(params, sat_ppm)
    lims = params.gen.lims;

    sat_freq = sat_ppm*1e-6*lims.B0*lims.gamma;
    
    if params.gen.field_strength == 7 || params.gen.field_strength == 7i 
        rf_fs = mr.makeGaussPulse(params.gen.fs_angle*pi/180,'system',lims,'Duration',5e-3,...
            'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
    elseif params.gen.field_strength == 9
        rf_fs = mr.makeGaussPulse(params.gen.fs_angle*pi/180,'system',lims,'Duration',8e-3,...
            'bandwidth',abs(sat_freq),'freqOffset',sat_freq);
    end
    
    % AMM: ToDo: Check correct values for spoiling gradient
    gx_fs = mr.makeTrapezoid('x',lims,'MaxGrad',18e-3*lims.gamma,'Area',-1/1e-3,'Duration',1.5e-3);
    gy_fs = mr.makeTrapezoid('y',lims,'MaxGrad',18e-3*lims.gamma,'Area',1/1e-3,'Duration',1.5e-3);
    gz_fs = mr.makeTrapezoid('z',lims,'MaxGrad',18e-3*lims.gamma,'Area',1/1e-3,'Duration',1.5e-3);
    gx_fs.delay = mr.calcDuration(rf_fs);
    gy_fs.delay = mr.calcDuration(rf_fs);
    gz_fs.delay = mr.calcDuration(rf_fs);

end