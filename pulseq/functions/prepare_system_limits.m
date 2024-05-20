function params = prepare_system_limits(params)

if params.gen.field_strength == 7
    % 7T, SC72 gradient
    lims = mr.opts('MaxGrad',65,'GradUnit','mT/m',...
        'MaxSlew',200,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 100e-6,'adcDeadTime', 10e-6, 'B0',6.98);  % To read it in VM I need rfDeadtime = 180e-6, 100e-6 for scanner
    params.gen.lims.forbidden_fq = [1100,550];
    params.gen.lims.forbidden_band = [300,100];
elseif params.gen.field_strength == 9
    % 9.4T, AC84 gradient
    lims = mr.opts('MaxGrad',80,'GradUnit','mT/m',...
        'MaxSlew',333,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 150e-6,'adcDeadTime', 10e-6, 'B0',9.38);
    params.gen.lims.forbidden_fq = [3000];
    params.gen.lims.forbidden_band = [100];
elseif params.gen.field_strength == 7i
    % 7T, impuse gradient
    lims = mr.opts('MaxGrad',198,'GradUnit','mT/m',...
        'MaxSlew',910,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 100e-6,'adcDeadTime', 10e-6, 'B0',6.98);
    params.gen.lims.forbidden_fq = [575,1060,1410,1920];
    params.gen.lims.forbidden_band = [70,70,70,120];
end
params.gen.lims = lims;

% Set gradient files for PNS check
if params.gen.field_strength == 7
    grad_file = './tools/pns_prediction/gradient_files/MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc';
elseif params.gen.field_strength == 7i
    grad_file = './tools/pns_prediction/gradient_files/MP_GPA_K2298_2250V_1250A_AC207_Base.asc';
else
    grad_file = '';
end
params.gen.grad_file = grad_file;

end

