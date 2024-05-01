function [xf, yf, sumres, ysine]=gradFreqPlot(t,y,Fs,axlabel)
% function [xf, yf, sumres, ysine]=gradFreqPlot(t,y,Fs,axlabel)
% take output from getDSV (t, y, Fs) to plot and 
% output the freq samples: xf, 
% fft of y normalised to total input length (numel(t)./Fs): yf
% sums of yf^2 in the acoustic resonance bands: sumres
% reference sine wave at acoustic resonance at 1 mT/m: ysine

%SC72 on 7 T Plus
acousticRes = [1100, 550, 367];
acousticBw = [300, 100, 55];
slewMax = 200; %mT/(m*ms)
gAmpMax = 42; %mT/m comments: 70 only for diffusion gradients

%AC84 on 9.4 T
%acousticRes = [3000];
%acousticBw = [100];
%slewMax = 333; %mT/(m*ms)
%gAmpMax = 50; %mT/m comments: 80 only for diffusion gradients

% allowed max oscillation amp
acousticMax = 4.999999999999; %mT/m

%axis label for plot
if nargin < 4 || isempty(axlabel) || ~ischar(axlabel)
    axlabel = '_u';
else
    axlabel = ['_' axlabel];
end

% bonus max grad and max slew rate check 
if any(abs(y)>gAmpMax)
    warning(['G' axlabel ' (' num2str(max(abs(y)),'%.1f') ' mT/m) exceeded the maximum amplitude of ' num2str(gAmpMax) ' mT/m!']);
end
if any(abs( diff(y)*Fs*1e-3 ) > slewMax)
    warning(['G' axlabel ' (' num2str(max(abs( diff(y)*Fs*1e-3 )),'%.1f') ' mT/(m*ms)) exceeded the maximum slew rate of ' num2str(slewMax) ' mT/(m*ms)!']);
end

%make reference of pure sine wave, worst case scenario 
ysine = sin(2.*pi.*repmat(t,[numel(acousticRes) 1]).* repmat((acousticRes.*1e-6).',[1 numel(t)] ));
%add the max amp set above in kT, and fft
ysinef = fftshift(fft(ysine.*acousticMax.*1e-6,[],2),2)./(numel(t)./Fs);

%convert units of y from mT to kT
y = y.*1e-6;

xf = [-numel(t)/2:(numel(t)/2-1)].*Fs./numel(t);
yf = fftshift(fft(y))./(numel(t)./Fs);

%calculate integral in banned range
resf = abs(xf)>(acousticRes.'-acousticBw.'/2) & abs(xf)<(acousticRes.'+acousticBw.'/2);
resfLorentzian = 1./(1+((xf-acousticRes.')./(acousticBw.'/4) ).^2 ) + 1./(1+((xf+acousticRes.')./(acousticBw.'/4) ).^2 ) ;
resf = resf.*resfLorentzian;
%yfres = resf.* abs(yf).^2;
%ysinefres = resf.*abs(ysinef).^2;
yfres = abs(resf.*yf).^2;
ysinefres = abs(resf.*ysinef).^2;
%resfLorzentian = 1./(1+((xf-acousticRes.')./(acousticBw.'/2) ).^2 ) + 1./(1+((xf+acousticRes.')./(acousticBw.'/2) ).^2 ) ;
%deltaf = xf(2)-xf(1);
%intres = sum(yfres, 2).*deltaf; %has boardening effect with the number of data points
sumres = sum(yfres, 2);
sumsineres = sum(ysinefres, 2);
sumresstr = ['sum G' axlabel '(f)^2: '];
bpass = true;
for ii=1:numel(acousticRes)
    sumresstr = [sumresstr '(' num2str(acousticRes(ii)) 'Hz): ' num2str(sumres(ii),3) ' '];
    if sumres(ii)>=sumsineres(ii) %equivalent of a pure sine wave at resonance with amplitude set above
        warning([num2str(ii) ' Acoustic Res f: ' num2str(acousticRes(ii)) ' Hz, bw: '...
            num2str(acousticBw(ii)) ' Hz. Sum(G' axlabel '(f)^2): '...
            num2str(sumres(ii),3) ' > ' num2str(sumsineres(ii),3) ' (' num2str(acousticMax,3) ' mT/m ref sine).' ]);
        sumresstr = [sumresstr '> ' num2str(sumsineres(ii),3) ' (' num2str(acousticMax,3) ' mT/m ref sine)!!! '];
        bpass = false;
    else
        disp([num2str(ii) ' Acoustic Res f: ' num2str(acousticRes(ii)) ' Hz, bw: '...
            num2str(acousticBw(ii)) ' Hz. Sum(G' axlabel '(f)^2): '...
            num2str(sumres(ii),3) ' < ' num2str(sumsineres(ii),3) ' (' num2str(acousticMax,3) ' mT/m ref sine).']);
        sumresstr = [sumresstr '< ' num2str(sumsineres(ii),3) ' (' num2str(acousticMax,3) ' mT/m ref sine) '];
    end 
    if ii~=numel(acousticRes)
        sumresstr = [sumresstr newline];
    end
end

linlim = 2.*[min(abs(yf).^2) max(abs(yf).^2)];
if linlim(1)==linlim(2)
    linlim(1) = linlim(1)-0.05;
    linlim(2) = linlim(2)+0.05;
end
figure();
for ii=1:2
    subplot(3,1,ii);
    plot(xf,abs(yf).^2, ...
        [1 1].*(acousticRes(1)), linlim, 'r-',...
    [1 1].*(-acousticRes(1)), linlim, 'r-',...
    [1 1].*(acousticRes(2)), linlim, 'r-',...
    [1 1].*(-acousticRes(2)), linlim, 'r-',...
        [1 1].*(acousticBw(1)/2+acousticRes(1)), linlim, 'r--',...
    [1 1].*(acousticBw(1)/2-acousticRes(1)), linlim, 'r--',...
    [1 1].*(acousticBw(2)/2+acousticRes(2)), linlim, 'r--',...
    [1 1].*(acousticBw(2)/2-acousticRes(2)), linlim, 'r--',...
    [-1 -1].*(acousticBw(1)/2+acousticRes(1)), linlim, 'r--',...
    [-1 -1].*(acousticBw(1)/2-acousticRes(1)), linlim, 'r--',...
    [-1 -1].*(acousticBw(2)/2+acousticRes(2)), linlim, 'r--',...
    [-1 -1].*(acousticBw(2)/2-acousticRes(2)), linlim, 'r--'...
        );
    
    ylabel(['G' axlabel '(f)^2/t_t_o_t (a.u.)']);
    ylim( linlim./2.*1.05 );
    xlabel('f (Hz)');
    if ii==2
        xlim([-1 1].*1.5.*max(acousticRes));
        title(sumresstr);
    elseif ii==1
        if bpass
            title(['G' axlabel ' passed resonance test.']);
        else
            title(['G' axlabel ' failed resonance test!!!']);
        end
    end
end %endfor
%plot grad in time
subplot(3,1,3);
plot(t*1e-6,y*1e6);
ylabel(['G' axlabel '(t) (mT/m)']); xlabel('t (s)');
xlim([t(1)-10 (t(end)+10)]*1e-6)
shg;



end %end function
