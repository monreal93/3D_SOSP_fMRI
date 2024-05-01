
function [rf_phase_offset,adc_phase_offset] = rf_adc_phase(params)
    %% 3D RF and ADC phase
    % This part is coming from ./Emily/RF_design/calcRFphase_3D.m
    % I don't know how this works for 3D; I kinda know now
    delta_x = 0;
    Gslice = 9.5571;
    Tpulse = 2560;
    Asym = 0.5;
    InitPhase = 90;
    
    f = floor(42.5756*delta_x*Gslice+0.5);
        RF.freqOffset = f;
        RF.InitPhaseSet = -(f*360/1e6*Tpulse*Asym)+InitPhase;
        RF.InitPhaseNeg = -(f*360/1e6*Tpulse*(1-Asym))-InitPhase;

    % Per partition and per interlave
    % AMM: Original phase...
    tmp = 1;
    for i = 1:params.gen.n(3) %MrProt.private.l_additionalslice+MrProt.sliceGroupList(1).sliceperslab
        for j = 1:params.spi.interl
            rf_phase(i,j) = RF.InitPhaseSet+25*(tmp^2)+175*tmp+300;
            rf_phase_offset(i,j) = mod(rf_phase(i,j),360)*pi/180;
            
            adc_phase_offset(i,j) = rf_phase(i,j)-RF.InitPhaseSet;
            adc_phase_offset(i,j) = mod(adc_phase_offset(i,j),360)*pi/180;
            tmp = tmp+1;
        end
    end
        
%     % AMM: Original phase...
%     for i = 1:params.gen.n(3) %MrProt.private.l_additionalslice+MrProt.sliceGroupList(1).sliceperslab
%         rf_phase(1,i) = RF.InitPhaseSet+25*i*i+175*i+300;
%         rf_phase(2,i) = RF.InitPhaseNeg+25*i*i+175*i+300;
%         rf_phase_offset(i) = mod(rf_phase(1,i),360)*pi/180;
%         
%         adc_phase_offset(i) = rf_phase(1,i)-RF.InitPhaseSet;
%         adc_phase_offset(i) = mod(adc_phase_offset(i),360)*pi/180;     
%     end
    
end