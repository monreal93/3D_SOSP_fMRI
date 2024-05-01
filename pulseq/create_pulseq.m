clear all; clc
% Change directory to path where sosp_vaso git has been cloned
cd /home/amonreal/Documents/PhD/PhD_2024/sosp_vaso/
addpath(genpath("./pulseq/functions"))

%% Adding paths 
addpath(genpath("/home/amonreal/Documents/PhD/tools/pulseq_1.4.1/"))                                    % Pulseq toolbox
addpath(genpath("/home/amonreal/Documents/PhD/tools/tOptGrad_V0.2/minTimeGradient/mex_interface/"))     % Minminue time gradieng Lusitg
addpath(genpath("/home/amonreal/Documents/PhD/tools/pns_prediction/"))                                  % PNS prediction
addpath(genpath("/home/amonreal/Documents/PhD/tools/check_grad_idea_Des/"))                             % Check Forbidden Fq
warning('OFF', 'mr:restoreShape')

%% Define parameters
folder_name = 'tmp';        % Day I am scanning
seq_name = 'sample';                % use sv/abc/sb_n (n for the diff scans at each day)
params.gen.seq = 4;                 % 1-VASO 2-ABC 3-Multi-Echo 4-BOLD
params.gen.field_strength = 7;      % Field Strength (7=7T,7i=7T-impuse_grad,9=9.4T,11=11.7T)
params.gen.plot = 0;                % Plot frequencies, PNS, etc...

% General parameters
params.gen.fov = [192 192 24].*1e-3;% FOV
params.gen.res = [0.8 0.8 1].*1e-3; % Nominal resolution
params.gen.fa = 0;                  % FA in degrees. Set to 0, to use Ernst Angle
params.gen.vfa = 0;                 % Variable FA (WIP)
params.gen.vfa_cutoff = 33;         % Variable FA cut-off (WIP) in degrees
params.gen.ernst_t1 = 1425e-3;      % T1 to calculate Ernst Angle (s) 7T=(WM-1220e-3)(GM-1800e-3)(blood=2587e-3) 9T=(1425e-3)(2100e-3)
params.gen.te = 0e-3;               % Set to 0 to shortest TE possible
params.gen.tr_delay = 0e-3;         % Delay between acquisitions in sec
params.gen.ro_type = 's';           % 's'-Spiral, 'c'-Cartesiaen (WIP)
params.gen.kz = 1;                  % Acceleration in Kz
params.gen.pf = 1;                  % Partial fourier in Kz
params.gen.fat_sat = 0;             % Fat saturation (1=yes,0=no)
params.gen.fs_angle = 0;            % Fat sat angle (0=default)
params.gen.skope = 0;               % Add skope sync scan and triggers, 0=N0, 1=sep scan, 2=concurrent(center partition), 3=1&2
params.gen.skope_sync = 0;          % Skope pre-scans, only added in skope seq
params.gen.dork = 0;                % extra adc's for DORK correction (WIP)
params.gen.kz_enc = 0;              % k-space partition encoding 0=linear,1=center-out For Cartesian now only linear encoding
params.gen.interl_enc = 0;          % Interleave encoding, 0=1-interl 1-plane -> 2-interl 1-plane, 1=1-inter 1-plane -> 1-interl 2-plane
params.gen.ph_oversampling = 0;     % Partition phase oversampling in %, to avoid partition phase-wrap
params.gen.echos = 1;               % Echos per RF pulse (WIP)
         
% Spiral parameters
params.spi.type = 0;                % spiral type 0=spiral-Out , 1=spiral-In, 3=In-Out (WIP), 4=In-Out kspace interleavead (WIP)
params.spi.in_out_order = 1;        % 0=In-Out same k-space path (separate vol.), 1=In-Out k-space path shift (WIP)
params.spi.rotate = 'none';         % Spiral rotation ('none','golden','180','120')
params.spi.increment = 'linear';    % Spiral increment mode (for now only 'linear') (WIP)
params.spi.max_grad  = 35;          % Peak gradient amplitude for spiral (mT/m)
params.spi.max_sr = 155;            % Max gradient slew rate for spiral (mT/m/ms)
params.spi.interl = 1;              % Spiral interleaves
params.spi.vd = 1.3;                % Variability density (accepts negative values)
params.spi.rxy = 3;               % In-plane (radial) undersampling
params.spi.rxy_az = 1;              % In-plane (azimuthal) undersampling (WIP)
params.spi.bw = 0e3;                % Spiral BW in Hz 0=Nyquist, highest=1000e3

% MT pulse parameters
params.mt.mt = 0;                   % Add MT pulse, 0 for reference scan without MT
params.mt.mt_prep = 0;              % MT (pre) dummy scans
params.mt.rf_spoil = 1;             % RF spoliing
params.mt.alpha = 225;              % default=225
params.mt.delta = -650;             % default=-650
params.mt.trf = 0.004;
params.mt.mt_rep = 100e-3;          % How often to play the MT pulse
params.mt.bold = 0;                 % Get BOLD reference acq after MT one (WIP)

% EPI parameters
params.epi.ry = 3;                  % In-plane undersampling
params.epi.pf = 6/8;                % In-plane PF, (1,7/8,6/8)
params.epi.te = [33.6 36.6 38.6 40]*1e-3+0.07; % Echo times for field map
params.epi.seg = 1;                 % EPI Segments
params.epi.tr_delay = 0;            % Delay after each TR, needed to reduce SAR
params.epi.bw_px = 1096;            % BW in Hz

% VASO parameters
params.vaso.foci = 1;               % FOCI inversion?
params.vaso.bold_ref = 1;           % BOLD reference volume
params.vaso.foci_ampl = 270;        % FOCI amplitude (270)
params.vaso.foci_dur = 10410e-6;    % FOCI duration (10410e-6)
params.vaso.foci_bw = 150;          % FOCI Bandwidth (150)
params.vaso.f_v_delay = 600e-3;     % FOCI-VASO delay (600e-3)
params.vaso.v_b_delay = 10e-3;      % VASO-BOLD delay (10e-3)
params.vaso.b_f_delay = 5e-3;       % BOLD-FOCI delay (5e-3)

%% Set system limits
if params.gen.field_strength == 7
    % 7T, SC72 gradient
    lims = mr.opts('MaxGrad',65,'GradUnit','mT/m',...
        'MaxSlew',200,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 100e-6,'adcDeadTime', 10e-6, 'B0',6.98);  % To read it in VM I need rfDeadtime = 180e-6, 100e-6 for scanner
elseif params.gen.field_strength == 9
    % 9.4T, AC84 gradient
    lims = mr.opts('MaxGrad',80,'GradUnit','mT/m',...
        'MaxSlew',333,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 150e-6,'adcDeadTime', 10e-6, 'B0',9.38); 
elseif params.gen.field_strength == 7i
    % 7T, impuse gradient
    lims = mr.opts('MaxGrad',198,'GradUnit','mT/m',...
        'MaxSlew',910,'SlewUnit','T/m/s',...
        'rfRingdownTime', 20e-6,'rfDeadtime', 100e-6,'adcDeadTime', 10e-6, 'B0',6.98); 
end
params.gen.lims = lims;

% Set gradient files for PNS check
if params.gen.field_strength == 7
    grad_file = '/home/amonreal/Documents/PhD/IDEA/gradient_files/MP_GPA_K2259_2000V_650A_SC72CD_EGA.asc';
elseif params.gen.field_strength == 7i
    grad_file = '/home/amonreal/Documents/PhD/IDEA/gradient_files/MP_GPA_K2298_2250V_1250A_AC207_Base.asc';
else
    grad_file = '';
end

%% Calculations and restrictions
[params,ro_blocks] = prepare_fix_parameters(params);

%% Create blocks
% TR-FOCI
if params.gen.seq == 1
    [B1_foci,phase,rf_complex,Gz_foci,fa] = tr_foci(params);
    rf_foci = mr.makeArbitraryRf(rf_complex,fa*pi/180, 'system', lims);%, 'Delay',2e-4);
end

% MT pulse
if params.gen.seq == 2
    MT = mr.makeGaussPulse(params.mt.alpha*pi/180,lims,'Duration',params.mt.trf,'FreqOffset',params.mt.delta);
end

% Fat sat
if params.gen.fat_sat
    sat_ppm = -3.45;
    [rf_fs,gx_fs,gy_fs,gz_fs] = prepare_fat_sat(params, sat_ppm);
end

% Kz blips
gz_blips = prepare_gz_blips(params);

% Spoilers
mag_spoil = 40e-3;     % mT
sr_spoil = 180;        % mT/m/s
[gx_spoil,gy_spoil,gz_spoil] = prepare_spoilers(params,mag_spoil,sr_spoil);

%% Preparing readout elements
[rf_phase_offset,adc_phase_offset] = rf_adc_phase(params);
if params.gen.ro_type == 's'
    [spiral_grad_shape,adcSamples,adcDwell,params] = prepare_spirals_rf_grad_adc(params);
        if params.spi.type == 0
            [gx,gy,~,~,adc,params] =  create_spirals_adc_pulseq(spiral_grad_shape,adcSamples,adcDwell,params);
        elseif params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
            [gx,gy,gx_pre,gy_pre,adc,params] =  create_spirals_adc_pulseq(spiral_grad_shape,adcSamples,adcDwell,params);
        end
elseif params.gen.ro_type == 'c'
    [gx,gy_blips,gx_pre,gy_pre,gy_blip_up,gy_blip_down,adc,params]  = create_epi_adc_pulseq(params);
end

% Flip angle, 3.5e-3 (approx rf + rephasing time)
tr_tmp = 3.5e-3+mr.calcDuration(gz_blips(1))+mr.calcDuration(mr.calcDuration(gx)*params.gen.echos)+mr.calcDuration(gx_spoil)+params.gen.te+params.gen.tr_delay; % aprox TR
if and(params.gen.ro_type=='s',or(params.spi.type == 1,params.spi.type == 2)); tr_tmp=tr_tmp+mr.calcDuration(gx_pre(1));end 
params = prepare_flip_angle(tr_tmp,gx, params);

% Create RF and Gz
[rf0,gz,gzReph,params] =  create_rf_gz_pulseq(params,tr_tmp);

%% Prepare Delays and triggers
if params.gen.te > 0; te_delay = mr.makeDelay(round(params.gen.te-((mr.calcDuration(rf0)/2)+mr.calcDuration(gzReph)+mr.calcDuration(gz_blips)),4)); end         % TE delay
if params.gen.tr_delay > 0; tr_delay = mr.makeDelay(params.gen.tr_delay); end                       % TR delay
if params.vaso.f_v_delay > 0; f_v_delay = mr.makeDelay(params.vaso.f_v_delay);          end         % FOCI-VASO delay
if params.vaso.v_b_delay > 0; v_b_delay = mr.makeDelay(params.vaso.v_b_delay);          end         % VASO-BOLD delay
if params.vaso.b_f_delay > 0; b_f_delay = mr.makeDelay(params.vaso.b_f_delay);          end         % BOLD-FOCI delay
% TR delay for Cartesian, to avoid PNS
if params.gen.ro_type == 'c'
    if params.epi.tr_delay > 0
        tr_delay_epi = mr.makeDelay(params.epi.tr_delay);  
    end
end

% Skope delays
if params.gen.skope ~= 0
    % ToDo: Need to confirm the times of this delays
    sk_pre_delay    = mr.makeDelay(0.5);
%     sk_int_delay    = mr.makeDelay(200e-6);  % original
    sk_int_delay = mr.makeDelay(mr.calcDuration(gz_blips(1))-10e-6); % To compensate for missing gz blip in center partition
    sk_post_delay   = mr.makeDelay(3);
    sk_min_tr_delay = gzReph(1).riseTime+gzReph(1).flatTime+gzReph(1).fallTime + ...
                      gz_blips(1).riseTime + gz_blips(1).flatTime + gz_blips(1).fallTime ...
                      + mr.calcDuration(gx(1)); %+ mr.calcDuration(gx_ramp(1));
    sk_min_tr_delay = mr.makeDelay(120e-3-sk_min_tr_delay);                         % Here I take Skope minTR = 110ms
    % Delay to account for the missing RF pulses...
    sk_no_rf_delay = mr.calcDuration(gz_fs)+mr.calcDuration(rf0(1))+mr.calcDuration(gzReph(1));
    sk_no_rf_delay = mr.makeDelay(sk_no_rf_delay);
    if params.gen.field_strength == 7i
        skope_trig = mr.makeDigitalOutputPulse('ext1','duration',10e-6,'system',lims);  % Skope trigger
    else
        skope_trig = mr.makeDigitalOutputPulse('osc0','duration',10e-6,'system',lims);  % Skope trigger
    end   
    
end

% Fieldmap scan delays
if params.gen.seq == 3
    if params.epi.te > 0
        tmp = (mr.calcDuration(rf0(1))/2)+mr.calcDuration(gzReph(1))+mr.calcDuration(gz_blips(1))+mr.calcDuration(gx_pre);
        tmp = tmp + (round(round(params.gen.n(2)/params.epi.ry)/2)*(mr.calcDuration(gx(1)))+mr.calcDuration(gy_blip(1)));
        for i=1:length(params.epi.te)
            fm_te_delay(i) = mr.makeDelay(round(params.epi.te(i)-tmp,4));
        end
    end
end

% Delay for repetitions without MT pulse, only for ABC (seq=2)
if params.gen.seq == 2
    no_mt_delay = mr.makeDelay(mr.calcDuration(MT)+mr.calcDuration(gx_spoil));
end

% No Gz blip delay
no_blip_delay = mr.makeDelay(mr.calcDuration(gz_blips(1)));

% No FatSat delay...
% fs_delay = mr.makeDelay(mr.calcDuration(rf_fs)+mr.calcDuration(gx_fs));

% External trigger for fMRI
ext_trig = mr.makeDigitalOutputPulse('ext1','duration',10e-6,'system',lims);        % External trigger
dummy_delay = mr.makeDelay(10e-3);                                                  % Delay to use as dummy anywhere

% Save partitions and segments in a new variable
last_part = params.gen.n(3);
last_seg = params.gen.seg;

%% Add blocks to Skope seq
if params.gen.skope == 1 || params.gen.skope == 3
    seq_sk=mr.Sequence();          % Create a new sequence object
%     seq_sk.addBlock(sk_pre_delay); % Skope initial delay
    % Skope sync scans
    for i_sk_pre = 1:params.gen.skope_sync
%         seq_sk.addBlock(mr.makeLabel('SET','ONCE',1));
%         seq_sk.addBlock(sk_pre_delay);
        seq_sk.addBlock(skope_trig);
        seq_sk.addBlock(sk_int_delay);     % Gradient free interval
        seq_sk.addBlock(no_blip_delay);
        if params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
            tmp_delay = mr.makeDelay(mr.calcDuration(gx_pre(1)));
            seq_sk.addBlock(tmp_delay);                       
        end
        tmp_delay = mr.makeDelay(mr.calcDuration(gx(1)));
        seq_sk.addBlock(tmp_delay,adc);
        tmp_delay = mr.makeDelay(mr.calcDuration(gz_spoil));
        seq_sk.addBlock(tmp_delay);
        % rewinder if multiple echos
        if params.gen.echos > 1
            tmp_delay = mr.makeDelay(mr.calcDuration(gx_pre(1)));
            seq_sk.addBlock(tmp_delay)
        end
        % Delay after sync scans
        if i_sk_pre == params.gen.skope_sync; seq_sk.addBlock(sk_post_delay); end
%         seq_sk.addBlock(mr.makeLabel('SET','ONCE',0));
    end
    % Partitions loop
    for i=1:last_part
        % Spiral/Cartesian Interleaves/segments loop
        % for m=1:params.gen.n(3)
            for j = 1:last_seg
                adc.phaseOffset = adc_phase_offset(i);
%                 seq_sk.addBlock(gzReph(i));
                if params.gen.te > 0; seq_sk.addBlock(te_delay); end % TE delay
                % EPI navigators
                if params.gen.ro_type == 'c'
                    for i_nav = 1:3
                        gx.amplitude = -gx.amplitude;
                        seq_sk.addBlock(gx,adc); 
                    end
                end
                seq_sk.addBlock(sk_no_rf_delay);      % Delay to account for the missing RF pulses... 
                % We only add the trigger in the center partition, remove if if I want it in all partitions
                if i == floor((params.gen.n(3)/2)+1)
                    sk0 = seq_sk.duration();           % To get Skope triger-ADC delay
                    seq_sk.addBlock(skope_trig);
                    seq_sk.addBlock(sk_int_delay);     % Gradient free interval
                end
                % Gz blip
                if i < floor((params.gen.n(3)/2)+1) % && params.gen.kz_enc == 0
                        tmp_blip = gz_blips(i,j);
                elseif or(i > floor((params.gen.n(3)/2)+1) && params.gen.kz_enc == 0, i > 1 && params.gen.kz_enc == 1)
                        tmp_blip = gz_blips(i-1,j);
                end
                if params.gen.ro_type == 's'
%                     if k==1 || k==params.gen.echos-1
                        if or(i == floor((params.gen.n(3)/2)+1) && params.gen.kz_enc == 0, i == 1 && params.gen.kz_enc == 1)
                            seq_sk.addBlock(no_blip_delay);
                        else
                            seq_sk.addBlock(tmp_blip);
                        end
%                     end
                    % if spiral-in or in-out, we need in-plane pre-phasing
                    if params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
                        seq_sk.addBlock(gx_pre(i,j),gy_pre(i,j));                       
                    end
                    if i == floor((params.gen.n(3)/2)+1); sk1 = seq_sk.duration();    end                % To get Skope triger-ADC delay
                    seq_sk.addBlock(gx(i,j),gy(i,j),adc);
                    seq_sk.addBlock(gz_spoil);
                    % rewinder if multiple echos
                    if params.gen.echos > 1; seq_sk.addBlock(gx_pre(i,j),gy_pre(i,j)); end
                    if params.gen.tr_delay > 0; seq_sk.addBlock(tr_delay); end % TR delay
                elseif params.gen.ro_type == 'c'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq_sk.addBlock(gx_pre,gy_pre)
                    else
                        seq_sk.addBlock(gx_pre,gy_pre,tmp_blip)
                    end
                    for k = 1:round(params.gen.n(2)/params.epi.ry)-((round(params.gen.n(2)/params.epi.ry)-round(params.gen.n(2)/params.epi.ry*params.epi.pf))/2)
                        gx.amplitude = -gx.amplitude;
                        if k == 1
                            seq_sk.addBlock(gx,gy_blip_up,adc);
                        elseif k == round(params.gen.n(2)/params.epi.ry*(params.epi.pf))
                            seq_sk.addBlock(gx,gy_blip_down,adc);
                        else
                            seq_sk.addBlock(gx,gy_blips,adc);
                        end
                    end
                    if gx.amplitude > 0
                        gx.amplitude = -gx.amplitude;
                    end
                    % seq_sk.addBlock(gx_spoil,gy_spoil,gz_spoil);
                    if params.gen.tr_delay > 0; seq_sk.addBlock(tr_delay); end % TR delay
                    if params.epi.tr_delay > 0; seq_sk.addBlock(tr_delay_epi); end
                end
%                 seq_sk.addBlock(sk_min_tr_delay); % Do I need this if I only trigger center partition?
            end
        % end
    end
    params.gen.skope_trig_adc_delay = sk1-sk0;
end

%% Actual scan
seq = mr.Sequence();
if params.gen.skope == 2 || params.gen.skope == 3
    % Skope sync scans
    for i_sk_pre = 1:params.gen.skope_sync
%         seq.addBlock(mr.makeLabel('SET','ONCE',1));
        seq.addBlock(sk_pre_delay);
        seq.addBlock(skope_trig);
        seq.addBlock(sk_int_delay);     % Gradient free interval
        seq.addBlock(no_blip_delay);
        if params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
            tmp_delay = mr.makeDelay(mr.calcDuration(gx_pre));
            seq.addBlock(tmp_delay);                       
        end
        tmp_delay = mr.makeDelay(mr.calcDuration(gx(1)));
        seq.addBlock(tmp_delay,adc);
        tmp_delay = mr.makeDelay(mr.calcDuration(gz_spoil));
        seq.addBlock(tmp_delay);
        % rewinder if multiple echos
        if params.gen.echos > 1
            tmp_delay = mr.makeDelay(mr.calcDuration(gx_pre));
            seq.addBlock(tmp_delay) 
        end
        % Delay after sync scans
        if i_sk_pre == params.gen.skope_sync; seq.addBlock(sk_post_delay); end
%         seq.addBlock(mr.makeLabel('SET','ONCE',0));
    end
end
% MT specific 
if params.gen.seq == 2
    if params.mt.mt_prep > 0 && params.mt.mt == 1
%         seq.addBlock(mr.makeLabel('SET','ONCE',1));
        for i_mt_prep = 1:params.mt.mt_prep
            seq.addBlock(MT);                                      % MT pulse 
            seq.addBlock(gz_spoil);
        end
        seq.addBlock(mr.makeDelay(1));                            % Delay after dummy MT-prep 
%         seq.addBlock(mr.makeLabel('SET','ONCE',0));
    end
    if params.mt.rf_spoil == 0
        rf_phase_offset = rf_phase_offset.*0;
        adc_phase_offset = adc_phase_offset.*0;
    end
end
% seq.addBlock(mr.makeLabel('SET','REP',0));
seq.addBlock(ext_trig);                                      % External trigger
if params.gen.seq == 1
    if params.vaso.foci; seq.addBlock(rf_foci); end              % FOCI
    if params.vaso.f_v_delay > 0; seq.addBlock(f_v_delay); end   % FOCI-VASO delay
end

for i_ro_blocks = 1:ro_blocks
    for i=1:last_part
        % Adding MT pulse every specified time
        if params.gen.seq == 2
            if params.mt.mt == 1 
                    seq.addBlock(MT);                                      % MT pulse 
                    seq.addBlock(gz_spoil);                                 % Original
            else
                    seq.addBlock(no_mt_delay);
            end
        end
        for j=1:last_seg
            if params.gen.fat_sat; seq.addBlock(rf_fs,gx_fs,gy_fs,gz_fs);   end      % fat-sat
%             if params.gen.fat_sat && j>1; seq.addBlock(fs_delay);   end        % fat-sat delayy...
            if i==1 && j==1 && i_ro_blocks==1; tr0 = seq.duration(); end         % save seq dur to calc TR
            rf = rf0(i);
            rf.phaseOffset = rf_phase_offset(i,j);
            adc.phaseOffset = adc_phase_offset(i,j);
            adc_post.phaseOffset = adc_phase_offset(i,j);
            seq.addBlock(rf,gz(i));
            if i==1; te0 = seq.duration(); end                           % save seq dur to calc TE
            seq.addBlock(gzReph(i));
            % EPI navigators
            if params.gen.ro_type == 'c'
                for i_nav = 1:3
                    gx.amplitude = -gx.amplitude;
                    seq.addBlock(gx,adc); 
                end
            end
            % Gz blip
            if i < floor((params.gen.n(3)/2)+1) % && params.gen.kz_enc == 0
                    tmp_blip = gz_blips(i,j);
            elseif or(i > floor((params.gen.n(3)/2)+1) && params.gen.kz_enc == 0, i > 1 && params.gen.kz_enc == 1)
                    tmp_blip = gz_blips(i-1,j);
            end
            % Echos loop
            for k=1:params.gen.echos
                if params.gen.ro_type == 's'
                    if k==1 || k==params.gen.echos-1
                        if or(i == floor((params.gen.n(3)/2)+1) && params.gen.kz_enc == 0, i == 1 && params.gen.kz_enc == 1)
                            if params.gen.skope == 2 || params.gen.skope == 3
                                seq.addBlock(skope_trig);
                                seq.addBlock(sk_int_delay);     % Gradient free interval
                            else
                                seq.addBlock(no_blip_delay);
                                if i_ro_blocks==1; params.vaso.ti1 = seq.duration();end
                                if i_ro_blocks==2; params.vaso.ti2 = seq.duration();end
                            end
                        else
                            seq.addBlock(tmp_blip);
                        end
                    end
                    if and(params.gen.seq == 3, params.epi.te > 0); seq.addBlock(fm_te_delay(i_ro_blocks)); end
                    if params.gen.te > 0; seq.addBlock(te_delay); end       % TE delay
                    % if spiral-in or in-out, we need in-plane pre-phasing
                    if params.spi.type == 1 || params.spi.type == 3 || params.spi.type == 4
                        seq.addBlock(gx_pre(i,j),gy_pre(i,j));                       
                    end
                    if i==1; te1 = seq.duration(); end                           % save seq dur to calc TE
                    seq.addBlock(gx(i,j),gy(i,j),adc);
%                     seq.addBlock(gx_ramp(i,j),gy_ramp(i,j));
                    if params.gen.dork; seq.addBlock(adc_post); end
                    % Spoil
%                     if k==params.gen.echos; seq.addBlock(gz_spoil); end  % Original
                    if k==params.gen.echos; seq.addBlock(gz_spoil,gx_spoil);end
                    % rewinder if multiple echos
                    if params.gen.echos > 1; seq.addBlock(gx_pre(i,j),gy_pre(i,j)); end
                    if params.gen.tr_delay > 0; seq.addBlock(tr_delay); end % TR delay
                elseif params.gen.ro_type == 'c'
                    if i == floor((params.gen.n(3)/2)+1)
                        seq.addBlock(gx_pre,gy_pre)
                    else
                        seq.addBlock(gx_pre,gy_pre,tmp_blip)
                    end
                    for k = 1:params.epi.n_lines
                        gx.amplitude = -gx.amplitude;
                        if k == 1
                            seq.addBlock(gx,gy_blip_up,adc);
                        elseif k == params.epi.n_lines
                            seq.addBlock(gx,gy_blip_down,adc);
                        else
                            seq.addBlock(gx,gy_blips,adc);
                        end
                        if i==1 && k==ceil((params.epi.n_lines/2*params.epi.pf)+1); te1=seq.duration(); end                           % save seq dur to calc TE
                    end
                    if gx.amplitude > 0
                        gx.amplitude = -gx.amplitude;
                    end
                    seq.addBlock(gz_spoil);
                    if params.gen.tr_delay > 0; seq.addBlock(tr_delay); end % TR delay
                    if params.epi.tr_delay > 0; seq.addBlock(tr_delay_epi); end
                end
            end
            if i==1 && j==1 && i_ro_blocks==1; tr1 = seq.duration(); end                           % save seq dur to calc TR
%         seq.addBlock(dummy_delay);
        end
%         seq.addBlock(mr.makeLabel('INC','PAR',1));
        dur1 = seq.duration();
    end
    if params.gen.seq == 1 && i_ro_blocks == 1
        if params.vaso.v_b_delay > 0; seq.addBlock(v_b_delay); end      % VASO-BOLD delay
    elseif params.gen.seq == 1 &&  i_ro_blocks == 2
        if params.vaso.b_f_delay > 0; seq.addBlock(b_f_delay); end   % BOLD-FOCI delay
    end
%     seq.addBlock(mr.makeLabel('INC','REP',1));
%     seq.addBlock(mr.makeLabel('SET','PAR',0));
end

%% Check timing of the sequence
[ok, error_report]=seq.checkTiming;
if (ok)
    fprintf('Timing check passed successfully\n');
else
    fprintf('Timing check failed! Error listing follows:\n');
    fprintf([error_report{:}]);
    fprintf('\n');
end

%% Getting k-space trajectories
ks_traj = create_ks_trajectory(seq,adc,params);

%% Adding some extra parameters to params directory
params = prepare_add_parameters(seq,ks_traj,gx,rf,adc,te0,te1,tr0,tr1,params);

%% Check accoustic resonance frequencies, Des script
check_accoustic_fq_pns(seq,params,grad_file)

%% Set definitions
seq.setDefinition('MaxAdcSegmentLength',params.gen.adc_split);
seq.setDefinition('FOV', params.gen.fov);
seq.setDefinition('Name', seq_name);
seq.setDefinition('baseResolution', 64) ;

%% Saving files
% Check if folder exist, if no it creates it
if exist(sprintf('./data/%s',folder_name)) == 0
    tmp = sprintf('./data/%s',folder_name);
    system(sprintf('mkdir %s',tmp));
    system(sprintf('mkdir %s/acq',tmp));
    system(sprintf('mkdir %s/analysis',tmp));
    system(sprintf('mkdir %s/ismrmd',tmp));
    system(sprintf('mkdir %s/ismrmd/2d',tmp));
    system(sprintf('mkdir %s/ismrmd/3d',tmp));
    system(sprintf('mkdir %s/raw',tmp));
    system(sprintf('mkdir %s/raw/twix',tmp));
    system(sprintf('mkdir %s/raw/dicom',tmp));
    system(sprintf('mkdir %s/raw/nifti',tmp));
    system(sprintf('mkdir %s/recon',tmp));
    system(sprintf('mkdir %s/recon/2d',tmp));
    system(sprintf('mkdir %s/recon/3d',tmp));
    system(sprintf('mkdir %s/tmp',tmp));
    system(sprintf('mkdir %s/sim',tmp));
end

namestr = strcat('./data/',folder_name,'/acq/',seq_name);
seq.write(strcat(namestr,'.seq'));
if params.gen.skope == 1 || params.gen.skope == 3
    seq_sk.write(strcat(namestr,'_sk.seq'));
end
if params.spi.type == 3 && params.gen.ro_type == 's'
    % ToDo : Fix this part so it works for more echos..
    ks_traj_full = ks_traj;
    ks_traj = ks_traj_full.e1;
    save(strcat(namestr,'_e1_ks_traj_nom.mat'),'ks_traj')
    ks_traj = ks_traj_full.e2;
    save(strcat(namestr,'_e2_ks_traj_nom.mat'),'ks_traj')
else
    save(strcat(namestr,'_ks_traj_nom.mat'),'ks_traj')
end
save(strcat(namestr,'_params.mat'),'params')