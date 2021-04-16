%-Description
%
%   Q_GUIDANCE computes the thrust pitch and yaw angle using Q-guidance
%   scheme. 
%
%-Inputs
%   
%   s0      current state 
%
%   sd      target state
%
%   t       remaining burn time (tdone-tnow) (s)
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
%   th      thrust pitch angle (rad)
%
%   psi     thrust yaw angle (rad)
%
%-Assumption 
%
%   Circular final orbit.
%
%-Reference
%
%   Battin, R. H., An Introduction to the Mathematics and Methods of 
%   Astrodynamics, Revised Edition, American Institute of Aeronautics and Astronautics, 1999.
%
%-&

function [th,psi] = Q_guidance(s0,sf,t,Ve,a,R,MU)

% Convert current and desired states to the inertial frame

[r0,v0] = lvlh2rv(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),R);

[rf,vf] = lvlh2rv(sf(1),sf(2),sf(3),sf(4),sf(5),sf(6),R);

% Calculate unit vectors

ir = r0/norm(r0);

in = cross(rf,vf)/norm( cross(rf,vf) );

tau = -Ve/(t-a);

S = [0 -in(3) in(2); in(3) 0 -in(1); -in(2) in(1) 0];

vr = S*r0*sqrt(MU/norm(r0)^3);

vg = vr - v0;

Q = sqrt(MU/norm(r0)^3)*S*( eye(3)-3/2*ir*ir.' );

p = -Q*vg;

ig = vg/norm(vg);

q = sqrt( tau^2 - norm(p)^2 + dot(p,ig)^2 );


aT = p + (q-dot(p,ig))*ig;

rhat = r0 / norm(r0);

lon = atan2( rhat(2), rhat(1) );

lat = asin( rhat(3) );

T = LVLH2N(lon,lat);

dum = inv(T)*aT;

ay = dum(1);

ax = dum(2);

az = dum(3);

aT = [ax; ay; az];

aTu = aT/norm(aT);

psi = atan2(aTu(3),aTu(1));

th  = asin(aTu(2));

end


% Support Functions

function [r,v] = lvlh2rv(x,y,z,u,v,w,R)

lon = x/R;

lat = z/R;

T = LVLH2N(lon,lat);

r = T*[R+y;0;0];

v = T*[v;u;w];

end

function T = LVLH2N(lon,lat)

T = [cos(lon)*cos(lat), -sin(lon), -cos(lon)*sin(lat);
     sin(lon)*cos(lat),  cos(lon), -sin(lon)*sin(lat);
     sin(lat),           0,         cos(lat)];

end