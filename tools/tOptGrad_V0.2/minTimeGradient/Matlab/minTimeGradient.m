function [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C,rv, g0, gfin, gmax, smax,T, ds,show)

%   [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C, RIV/RV, g0, gfin, gmax, smax,T, ds, show)
%   [C,time,g,s,k, phi, sta, stb] = minTimeGradient(C, [], [], [], gmax, smax, T, [], [])
%  
%   Given a k-space trajectory C = [x y z], gradient and slew constraints,
%   This function will return a new parameterization that will meet these constraints
%   while getting from one point to the other in minimum time.
%   The constraints can be either magnitude (|g| < Gmax and |s| < Smax) or as a set of
%   individual for each direction (|gx|, |gy|, |gz| < Gmax and |sx|, |sy|, |sz| < Smax).
%   It will thus either call the function minTimeGradientRIV.m or
%   minTimeGradientRV.m based on the specified solution type. 
%  
%   Input Values :
%  
%   C       -   The Curve in k-space given in any parametrization [1/cm]
%               C should be inputed as an Nx2 or Nx3 trajectory.
%   RIV/RV   -  0 for rotationally invariant solution (magnitude constraints),
%                1 for rotationally variant solution (individual gradient/slew constraints).
%                Leave empty for default (0 for the rotationally invariant solution).
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

if nargin < 2
    rv = [];
end

if nargin<3
    g0 = [];
end

if nargin<4
    gfin = [];
end

if nargin<5
    gmax = 4;
end

if nargin<6
    smax = 15;
end

if nargin<7
    T = 4e-3;
end

if nargin < 8
    ds = [];
end

if nargin<9
    show=0;
end

if isempty(rv)
    rv = 0;
end

if rv == 1
    disp('Rotationally Variant Solution')
    [C, time, g, s, k, phi, sta, stb] = minTimeGradientRV(C, g0, gfin, gmax, smax, T,ds, show);
else
     disp('Rotationally Invariant Solution')
    [C, time, g, s, k, phi, sta, stb] = minTimeGradientRIV(C, g0, gfin, gmax, smax, T,ds, show);
end
