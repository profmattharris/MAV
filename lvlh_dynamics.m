%- Description
%
%   LVLH_DYNAMICS defines the nonlinear differential equations of motion
%   for the MAV written in a Local Horizontal Local Vertical (LVLH) frame.
%   x and y are the in-plane curvilinear distance and altitude, u and v are
%   the in-plane horizontal and vertical velocity components, z is the out-
%   of-plane curvilinear distance and w is the out-of-plane velocity. tau
%   is the thrust acceleration.
%
%-Inputs
%   t               current time (s)
%
%   s               initial state vector      
%
%   R               radius of Mars (m)
%
%   MU              gravitational parameter (m^3/s^2)
%
%   Ve              exhaust velocity (m/s)
%
%   a               ratio (alpha = m0/b) (s)
%
%   th              thrust pitch angle (rad)
%
%   psi             thrust yaw angle (rad)
%
%   tstart          iginition time (s)
%
%-Output
%
%   sdot           the derivative of the state vector
%
%-Reference
%
%   Hull, D. G. and Harris, M. W., "Optimal Solutions for Quasi-Planar 
%   Ascent Over a Spherical Moon," Journal of Guidance, Control, and 
%   Dynamics, Vol. 35, July- August 2012, pp. 1218-1223.
%-&




function sdot = lvlh_dynamics(t,s,R,MU,Ve,a,th,psi,tstart)

x = s(1); y = s(2); z = s(3);

u = s(4); v = s(5); w = s(6);

r = R + y;

g = MU/r^2;

tau = -Ve/(t-tstart-a);

xdot = R*u/r;

ydot = v;

zdot = R*w/r;

udot = tau*cos(th)*cos(psi) + u*w/r*tan(z/R) - u*v/r;

vdot = tau*sin(th) - g + u^2/r + w^2/r;

wdot = tau*cos(th)*sin(psi) - u^2/r*tan(z/R) - v*w/r;

sdot = [xdot; ydot; zdot; udot; vdot; wdot];

end