clear endtime simdir starttime sumres_xyz;
ctopdir = '/home/amonreal/Documents/PhD/tools/check_grad_idea_Des';
simdir = uigetdir(ctopdir, 'pick the simulation results directory');
if simdir==0
    error('need directory input');
end
if ~exist(simdir,'dir')
    error('directory does not exist');
end
cd(simdir);

gaxes = ['X' 'Y' 'Z'];
gfilehead = 'SimulationProtocol_GR';
gfiletail = '.dsv';

for ii=1:numel(gaxes)
    if ~exist([simdir '/' gfilehead gaxes(ii) gfiletail],'file')
        error([simdir '/' gfilehead gaxes(ii) gfiletail 'does not exist']);
    end
end %endfor gaxex

t_tot = getDSV([gfilehead gaxes(1) gfiletail]);
tsn = inputdlg( {'start time (s)' 'end time (s)'}, 'choose sim time range', [1 35], {'0.00000', num2str(t_tot.*1e-6,'%.5f')});
if isempty(tsn)
    error('need start and end time inputs');
else  
    starttime = str2double(tsn{1});
    endtime = str2double(tsn{2});
end
if starttime < 0 || starttime >= endtime || endtime <= 0 || round(endtime*1e6) > t_tot || round(starttime*1e6) >= t_tot 
    error('errors in start or end time')
end
clear tsn t_tot;
disp(['checking from ' num2str(starttime,'%.6f') ' s to ' num2str(endtime,'%.6f') ' s.']);

disp(' ');
choice = questdlg('Save figures?','Save figures?','Yes','No','Yes');
switch choice
    case 'Yes'
        bsavefig = true;
    otherwise
        bsavefig = false;
end
clear choice;

sumres_xyz = [];
for ii=1:numel(gaxes)
    disp(['loading ' gfilehead gaxes(ii) gfiletail '...']);
    [t,y,Fs]=getDSV(['SimulationProtocol_GR' gaxes(ii) '.dsv'], endtime);
    tfilt = (t >= starttime*1e6);
    t = t(tfilt);
    y = y(tfilt);
    disp(['loaded data from ' num2str(t(1)*1e-6,'%.6f') ' s to ' num2str(t(end)*1e-6,'%.6f') ' s.']);
    [~, ~, sumres, ~] = gradFreqPlot(t,y,Fs,gaxes(ii));
    set(gcf,'Position',[36 8 577 781]);
    sumres_xyz = [sumres_xyz sumres];
    if bsavefig
        saveas(gcf,['SimulationProtocol_GR' gaxes(ii) '_' num2str(t(1)) '_us_to_' num2str(t(end)) '_us.svg' ]);
    end
    clear t y Fs sumres tfilt;
    disp(' ');
end %endfor gaxex
cd(ctopdir);
clear gaxes gfilehead gfiletail ii bsavefig ctopdir;