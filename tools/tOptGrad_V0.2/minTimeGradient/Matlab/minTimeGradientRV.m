function [C,time,g,s,k, phi, sta, stb] = minTimeGradientRV(C,g0, gfin, gmax, smax,T, ds, show)

% [C,time,g,s,k] = minTimeGradient(C, g0, gfin, gmax, smax,T,ds, show)
% 
% This function computes the rotationally variant solution to the time
% optimal gradient design. 
%  
%   C       -   The Curve in k-space given in any parametrization [1/cm]
%               C should be inputed as an Nx2 or Nx3 trajectory.
%   g0      -   Initial gradient amplitude (leave empty for g0 = 0)
%   gfin    -   Gradient value at the end of the trajectory. If not possible, 
%               the result would be the largest possible ampltude.
%               (Leave empty if you don't care to get maximum gradient.)
%   gmax    -   Maximum gradient [G/cm] (3.9 default)
%   smax    -   Maximum slew [G/Cm/ms]  (14.5 default)
%   T       -   Sampling time interval [ms] (4e-3 default)
%   ds      -   step size for ODE integration, leave empty to use default value
%   show    -   Show plots while optimizing (Warning: This will make the
%               process considerably slower!)
%   
% return values:
%   C       - reparametrized curve, sampled at T[ms]
%   time    - total time to get to the end
%   g       - gradiet waveform [G/cm]
%   s       - slew rate [G/cm/ms]
%   k       - exact k-space corresponding to gradient g (This function reparametrizes
%             C, then takes a derivative. Numerical errors in the derivative can lead to 
%             deviation.  
%   phi     - Geometry constraints on the amplitude vs. arclength
%   sta     - Solution for the forward ODE
%   stb     - Solution for the backward ODE
if nargin<1
    error('You gotta give me a curve!');
end

if nargin<2
    g0 = [];
end

if nargin<3
    gfin = [];
end

if nargin<4
    gmax = 4;
end

if nargin<5
    smax = 15;
end

if nargin<6
    T = 4e-3;
end

if nargin < 7
    ds = [];
end

if nargin<8
    show=0;
end

dt = T;

gamma = 4.257;

disp('Const arc-length parametrization');

% represent the curve using spline with parametrization p
Lp = length(C);
p = [0:Lp-1].';
PPx = spline(p,C(:,1));
PPy = spline(p,C(:,2));
PPz = spline(p,C(:,3));

% interpolate curve for gradient accuracy
dp = 1e-1;
CCx = ppval(PPx, [0:dp:Lp-1]');
CCy = ppval(PPy, [0:dp:Lp-1]');
CCz = ppval(PPz, [0:dp:Lp-1]');

Cpx = (CCx([2:end,end]) - CCx)/dp;
Cpx(end) = Cpx(end-1);
Cpx(1) = (CCx(2)-CCx(1))/dp;

Cpy = (CCy([2:end,end]) - CCy)/dp;
Cpy(end) = Cpy(end-1);
Cpy(1) = (CCy(2)-CCy(1))/dp;

Cpz = (CCz([2:end,end]) - CCz)/dp;
Cpz(end) = Cpz(end-1);
Cpz(1) = (CCz(2)-CCz(1))/dp;

Cp = (Cpx.^2 + Cpy.^2 + Cpz.^2).^0.5;
% find length of curve

s_of_p = cumtrapz(abs(Cp))*dp;
L = s_of_p(end);


% decide ds and compute st for the first point
stt0 = (gamma*smax) ; % always assumes first point is max slew
st0 = stt0*dt/2; % start at half the gradient for accuracy close to g=0
s0 = st0*dt;
if isempty(ds)
    ds = s0/1.5; % smaller step size for numerical accuracy
end
s = [0 : ds : L].';
s_half = [0:ds/2:L].';
sta = s*0;

if isempty(g0)
    g0 = 0;
end

sta(1) = min(max(g0*gamma+st0), gamma*gmax);
p_of_s_half = interp1(s_of_p, [0:dp:Lp-1]', s_half, 'spline'); 
p_of_s = p_of_s_half(1:2:end);


disp('Compute geometry dependent constraints');
% compute constraints (forbidden line curve)
[phi,kx, ky, kz, xp, yp, zp] = sdotMax(PPx, PPy, PPz, p_of_s_half,s_half, gmax);
k = ((kx).^2 + (ky).^2 + (kz).^2).^(0.5);
k = k([1:end,end,end]); % extend for the Runge-Kutte method
kx = kx([1:end,end,end]);
ky = ky([1:end,end,end]);
kz = kz([1:end,end,end]);

xp = xp([1:end,end,end]);
yp = yp([1:end,end,end]);
zp = zp([1:end,end,end]);

disp('Solve ODE forward');
% solve ODE forward
for n=2:length(s)
    dstds = RungeKutte(ds,sta(n-1),kx([(n-1)*2-1:(n-1)*2+1]), ky([(n-1)*2-1:(n-1)*2+1]), kz([(n-1)*2-1:(n-1)*2+1]),...
                                   xp([(n-1)*2-1:(n-1)*2+1]), yp([(n-1)*2-1:(n-1)*2+1]), zp([(n-1)*2-1:(n-1)*2+1]), smax);
    tmpst = sta(n-1) + dstds;
    
    sta(n) = min(tmpst, phi(n*2-1));
    if mod(n,1000)==1 & show
        s(n)/L*100
        figure(100), plot(s,sta,linspace(s(1),s(end),length(phi)),phi); , axis([0,s(end),0,gmax*gamma*sqrt(2)]);, drawnow
    end
end


stb = 0*s;

if isempty(gfin)
    stb(end) = sta(end);
else
    stb(end) = min(max(gfin*gamma,st0), gamma*gmax);
end

% solve ODE backwards
disp('Solve ODE backwards');
for n=length(s)-1:-1:1
    dstds = RungeKutte(ds, stb(n+1), kx([(n+1)*2-1:-1:(n)*2-1]), ky([(n+1)*2-1:-1:(n)*2-1]), kz([(n+1)*2-1:-1:(n)*2-1]), ...
                                     xp([(n+1)*2-1:-1:(n)*2-1]), yp([(n+1)*2-1:-1:(n)*2-1]), zp([(n+1)*2-1:-1:(n)*2-1]), smax);
    tmpst = stb(n+1) + dstds;
    if isreal(tmpst)
        stb(n) = min(tmpst,phi(n*2-1));
    else
        stb(n) = phi(n*2-1);
    end
    if mod(n,1000)==1 & show
        s(n)/L*100
        figure(100), plot(s,stb,linspace(s(1),s(end), length(phi)),phi,s,sta);, axis([0,s(end),0,gmax*gamma*sqrt(2)]);,drawnow
    end
end

disp('Final Interpolations')
% take the minimum of the curves
st_of_s = min([sta,stb],[],2);

% compute time
t_of_s = cumtrapz(1./st_of_s*ds);
t = 0:dt:t_of_s(end);

s_of_t = interp1(t_of_s, s, t,'spline');
p_of_t = interp1(s, p_of_s, s_of_t,'spline');

Cx = ppval(PPx, p_of_t)';
Cy = ppval(PPy, p_of_t)';
Cz = ppval(PPz, p_of_t)';
C = [Cx Cy Cz];

gx = diff(Cx)/gamma/dt; gx=[gx; gx(end) + gx(end)-gx(end-1)]; 
gy = diff(Cy)/gamma/dt; gy=[gy; gy(end) + gy(end)-gy(end-1)]; 
gz = diff(Cz)/gamma/dt; gz=[gz; gz(end) + gz(end)-gz(end-1)]; 
g = [gx gy gz];

kx = cumtrapz(gx)*dt*gamma;
ky = cumtrapz(gy)*dt*gamma;
kz = cumtrapz(gz)*dt*gamma;
k = [kx ky kz];

sx = diff(gx)/dt;
sy = diff(gy)/dt;
sz = diff(gz)/dt;
s = [sx sy sz];

time = t(end);

disp('Done');

function sdotdot = sdotdot(xpp, ypp, zpp, xp, yp, zp, st, smax)
gamma = 4.257;
sxx = (-xpp*st*st + gamma * smax) / xp;
syy = (-ypp*st*st + gamma * smax) / yp;
szz = (-zpp*st*st + gamma * smax) / zp;
sdotdot = min([sxx, syy, szz]);

function [res] = RungeKutte(ds,st,xpp, ypp, zpp, xp, yp, zp, smax)
% Solve ODE using Runge-Kutte method.

k1 = ds * (1/st)* sdotdot(xpp(1), ypp(1), zpp(1), xp(1), yp(1), zp(1), st, smax);
k2 = ds * (1/(st + k1/2)) * sdotdot(xpp(2), ypp(2), zpp(2), xp(2), yp(2), zp(2), st+k1/2, smax);
k3 = ds * (1/(st + k2/2)) * sdotdot(xpp(2), ypp(2), zpp(2), xp(2), yp(2), zp(2), st+k2/2, smax);
k4 = ds * (1/(st + k3/3)) * sdotdot(xpp(3), ypp(3), zpp(3), xp(3), yp(3), zp(3), st+k3/2, smax);

res = k1/6 + k2/3 + k3/3 + k4/6;

function [sdot, kx, ky, kz, Csx, Csy, Csz] = sdotMax(PPx, PPy, PPz, p_of_s, s, gmax)

gamma = 4.257;

s = s(:);
dp_p = p_of_s([2:end,end]) - p_of_s; , dp_p(end) = dp_p(end-1);
dp_m = p_of_s - p_of_s([1,1:end-1]);, dp_m(1) = dp_m(2);
ds_p = s([2:end,end]) - s; , ds_p(end) = ds_p(end-1);
ds_m = s - s([1,1:end-1]);, ds_m(1) = ds_m(2);


Cs_px = (ppval(PPx,p_of_s + dp_p) - ppval(PPx, p_of_s))./ds_p;
Cs_mx = (ppval(PPx,p_of_s) - ppval(PPx, p_of_s-dp_m))./ds_m;
Csx = Cs_px/2 + Cs_mx/2;
Cssx = (Cs_px - Cs_mx)./(ds_m/2+ds_p/2);

Cs_py = (ppval(PPy,p_of_s + dp_p) - ppval(PPy, p_of_s))./ds_p;
Cs_my = (ppval(PPy,p_of_s) - ppval(PPy, p_of_s-dp_m))./ds_m;
Csy = Cs_py/2 + Cs_my/2;
Cssy = (Cs_py - Cs_my)./(ds_m/2+ds_p/2);

Cs_pz = (ppval(PPz,p_of_s + dp_p) - ppval(PPz, p_of_s))./ds_p;
Cs_mz = (ppval(PPz,p_of_s) - ppval(PPz, p_of_s-dp_m))./ds_m;
Csz = Cs_pz/2 + Cs_mz/2;
Cssz = (Cs_pz - Cs_mz)./(ds_m/2+ds_p/2);

kx = abs(Cssx);   ky = abs(Cssy);   kz = abs(Cssz);
Csx = abs(Csx);
Csy = abs(Csy);
Csz = abs(Csz);
% fix edge numerical problems
kx(end) = kx(end-1);  ky(end) = ky(end-1);  kz(end) = kz(end-1);
kx(1) = kx(2);    ky(1) = ky(2);    kz(1) = kz(2);

% curve curvature dependent constraint
sdot1 = gamma*gmax./(Csx+eps);
sdot2 = gamma*gmax./(Csy+eps);
sdot3 = gamma*gmax./(Csz+eps);

% calc total constraint
sdot = min([sdot1, sdot2, sdot3],[],2);