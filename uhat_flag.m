%-Description
%
%   UHAT_FLAG implements the uhat-trigger.
%
%-Inputs
%   
%   s0      current state 
%
%   sf      target state
%
%   tf      SRM2 burn time (s)
%
%   Ve      SRM2 exhaust velocity (m/s)
%
%   a       ratio (m0/b) (s)
% 
%   R       radius of Mars (m)
%
%   MU      gravitaional parameter (m^3/s^2)
%
%-Outputs
%
%   flag            scalar that stores either 0 or 1. 
%                   flag = 0 and flag = 1 corresponds to trigger being 
%                   inactive and active, respectively. 
%
%-Assumption 
%
%   Quasiplanar flight, small thrust angles
%
%-Reference
%
%   Hull, D. G. and Harris, M. W., "Optimal Solutions for Quasi-Planar 
%   Ascent Over a Spherical Moon," Journal of Guidance, Control, and 
%   Dynamics, Vol. 35, July- August 2012, pp. 1218-1223.
%-&


function flag = uhat_flag(s0,sf,tf,Ve,a,R,MU)

 
x0 = s0(1); y0 = s0(2); z0 = s0(3); u0 = s0(4); v0 = s0(5); w0 = s0(6);

xf = sf(1); yf = sf(2); zf = sf(3); uf = sf(4); vf = sf(5); wf = sf(6);


[r,v] = lvlh2rv(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),R);

[~,~,~,~,~,nu] = rv2orbel(r,v,MU);

nu = nu *180/pi;

gm = MU/( R+y0 )^2 - u0^2/R;

 

[F,G] = FG_int(u0,tf,Ve,a,R);

Vy = vf-v0 + gm*tf - F;

Vz = wf-w0;

Y  = yf-y0-v0*tf + gm*tf^2/2 - G;

Z  = zf-z0-w0*tf;

 

L = -Ve*log(1-tf/a);

S = (-a+tf)*L + Ve*tf;

J = -S + L*tf;

Q = -Ve*tf^2/2 + a*S;
 
L2 = (Vy*S-Y*L)/(L*Q-J*S);

C2 = (Vy*Q-Y*J)/(L*Q-J*S);

L3 = (Vz*S-Z*L)/(L*Q-J*S);

C3 = (Vz*Q-Z*J)/(L*Q-J*S);

 
uhat = u_est(u0,L2,C2,L3,C3,tf,Ve,a);


if uhat >= uf || nu >= 180

    flag = 1;

else

    flag = 0;

end


end

% Support functions

function [r,v] = lvlh2rv(x,y,z,u,v,w,R)

lon = x/R;

lat = z/R;

TT = LVLH2N(lon,lat);

r = TT*[R+y;0;0];

v = TT*[v;u;w];

end


function TT = LVLH2N(lon,lat)

TT = [cos(lon)*cos(lat), -sin(lon), -cos(lon)*sin(lat);
     sin(lon)*cos(lat),  cos(lon), -sin(lon)*sin(lat);
     sin(lat),           0,         cos(lat)];

end


function [F,G] = FG_int(u0,tf,Ve,alpha,R)

% Compute the centrifugal acceleration integrals F and G

xx    = 1 - tf/alpha;

xx2   = xx^2;

lnxx  = log(xx);

lnxx2 = lnxx^2;

 

I1 = -alpha*(xx*lnxx-xx+1);

I2 = -alpha*(xx*lnxx2-2*xx*lnxx+2*xx-2);

I5 = -(alpha/2)*(xx2*lnxx-xx2/2+1/2);

I6 = -(alpha/2)*(xx2-1);

I3 = -alpha*(I5-I6+tf);

I7 = -(alpha/2)*(xx2*lnxx2-xx2*lnxx+xx2/2-1/2);

I4 = -alpha*(I7-2*I5+2*I6-2*tf);

 

F = (u0^2*tf-2*u0*Ve*I1+Ve^2*I2)/R;

G = (u0^2*tf^2/2-2*u0*Ve*I3+Ve^2*I4)/R;

end


function uhat = u_est(u0,L2,C2,L3,C3,tf,Ve,a)

% Compute uhat using quadrature

t = linspace(0,tf);

L5 = -L2*t + C2;

L6 = -L3*t + C3;

denpsi  = sqrt(1+L6.^2);

sinpsi  = L6./denpsi;

cospsi  = 1./denpsi;

psi     = atan2(sinpsi,cospsi);

denth = sqrt(1+L5.^2+L6.^2);

sinth = L5./denth;

costh = sqrt(1+L6.^2)./denth;

th    = atan2(sinth,costh);

uhat = u0 + trapz(t,-Ve./(t-a).*cos(th).*cos(psi));

end
