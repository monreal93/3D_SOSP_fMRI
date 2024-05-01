%GETDSV reads a Siemens IDEA file in DSV format
% tmax = getDSV(f)
% f     filename
% tmax  maximum time in seconds
%
% [t, y, Fs] = getDSV(f, tmax)
% t     time of samples in microseconds
% y     sample values
% Fs    sampling rate in Hz
%
% Author: Christian Labadie, 2008-2012

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [t, y, Fs]=getDSV(f,tmax)

fp=fopen(f,'r');
key='';
while ~strcmp(key,'[VALUES]')
    st=fgetl(fp);
    if ~ischar(st), break; end
    [key,val]=strtok(st,'='); val=strtok(val,'=');
    switch key
        case 'SAMPLES',     ny=str2double(val);
        case 'HORIDELTA',   dt=str2double(val);  % microseconds
        case 'VERTFACTOR',  dy=1/str2double(val);
        case 'STARTTIME',   tbeg=str2double(val) * 1000; % in milli-seconds
    end
end

if nargin == 1
    t = ny * dt; % tmax in file
    fclose(fp);
    return
end
dsv = textscan(fp, '%n','commentStyle', ';');
dsv = dsv{1}(:);
fclose(fp);

[d f]=fileparts(f);
fprintf('%s: ',f);

% DSV encoding :
% The values are signed 32-bit integers. The real value will be obtained
% after division by VERTFACTOR. If the same sample value is repeated
% several times, the value is only listed twice, and the third value
% specifies the number of repetitions. Only the first sample is stored as
% an absolute value. All other values are deltas to the preceeding value.
% (This way it is possible to compress linear slopes.)
%
% [ 1,1,1,1,5,5,8,8,8,10,11] will be encoded as [1,1,2,5,5,0,8,8,1,10,11]

% convert DSV
ndsv=length(dsv);
ny = min(ny, floor(tmax*1e6/dt));
nyd = 0;
t = (1:ny) * dt;
y = t * 0;
i = 1;
n = 1;
y(1)=dsv(1);
prev=dsv(1);
while i<ndsv && n < ny
    i = i + 1;
    y(n+1) = y(n) + dsv(i);
    if dsv(i) ~= prev
        n = n + 1;
        prev = dsv(i);
    else
        nrep = dsv(i+1);
        if nrep == 0        % zero repetition
            i = i + 2;
            n = n + 1;
        else
            n = n + 1;
            y(n+1:n+nrep) = y(n) + (1:nrep) * prev; % to speed up this loop
            n = n + nrep;
            i = i + 2;
        end

        if i <= ndsv
            y(n+1) = y(n) + dsv(i);
            n = n + 1;
            prev = dsv(i);
        end
    end

    if n - nyd >= ny/10, fprintf(1,'%.0f%% ', n/ny*100); nyd = n; end

end
if n ~= ny, y = y(1:ny); end
y = y * dy;
Fs = 1e6/dt;  % in Herz

fprintf(1,'\n')
