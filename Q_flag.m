%-Description
%
%   Q_FLAG implements the Q-trigger as described in the paper by Everett.
%    
%-Inputs
%   
%   s0          current state 
%
%   sf          target state
% 
%   R           radius of Mars (m)
%
%   MU          gravitaional parameter (m^3/s^2)
%
%   dV_SRM2     delta-v capability of SRM2 (m/s)
%
%-Outputs
%
%   flag        scalar that stores either 0 or 1. 
%               flag = 0 and flag = 1 corresponds to trigger being 
%               inactive and active, respectively. 
%
%-Assumption 
%
%   Circular final orbit.
%
%-Reference
%
%   Everett, J. M., "A Generalized Guidance Approach to In-Space Solid-Propellant Vehicle
%   Maneuvers," AAS Guidance, Navigation, and Control Conference, 2020.
%
%-&


function flag = Q_flag(s0,sf,R,MU,dV_SRM2)

% Convert current and desired states to the inertial frame

[r0,v0] = lvlh2rv(s0(1),s0(2),s0(3),s0(4),s0(5),s0(6),R);

[rf,vf] = lvlh2rv(sf(1),sf(2),sf(3),sf(4),sf(5),sf(6),R);

% Calculate current eccentric anomaly and mean anomaly

[sma0,ecc0,inc0,W0,w0,nu0] = rv2orbel(r0,v0,MU);

E0 = sqrt((1-ecc0)/(1+ecc0))*tan(nu0/2);

E0 = 2*atan2(E0,1);

M0 = E0-ecc0*sin(E0);

nu0 = nu0*180/pi;

n = sqrt(MU/sma0^3);

if nu0 >= 177 && nu0 < 180
    
    % Calculate desired semi-major axis
    
    [smaf,~,~,~,~,~] = rv2orbel(rf,vf,MU);

    % Calculate unit vectors
    
    ix = r0/norm(r0);
    
    iy = -cross(rf,vf)/norm( cross(rf,vf) );
    
    iz = cross(ix,iy)/norm(cross(ix,iy));

    tig = 0;
    
    error = 1.0e12;
    
    tol = 1.0e-5;
    
    tol_Kep = 1.0e-8;

    while error > tol

        M = M0 + n*tig;
        
        E = KepEqE(M,ecc0,tol_Kep);
        
        nu = 2*atan2(sqrt(1+ecc0)/sqrt(1-ecc0)*tan(E/2),1);
        
        [r,v] = orbel2rv(sma0,ecc0,inc0,W0,w0,nu,MU); 

        vd = sqrt(MU*(2/norm(r)-1/smaf))*iz;
        
        vgo = vd - v;

        dum = -MU/(norm(vd)*norm(r)^2)*((r/norm(r))'*v)*iz ...
            - norm(vd)*((v'*iz)/(r'*ix))*ix + MU/(norm(r)^3)*r;
        
        dvgo_dt = (vgo/norm(vgo))'*dum;

        tig_new = tig - (norm(vgo) - dV_SRM2)/dvgo_dt;
        
        error = abs(tig_new - tig);

        tig = tig_new;

    end

    if tig <= 0 
        
        flag = 1;
        
    else
        
        flag = 0;
        
    end

elseif nu0 >= 180   
    
    flag = 1;   
    
else    
    
    flag = 0;
    
end

end

% Support Functions

function E = KepEqE(M,e,tol)

% Solves Kepler's equation

    if (-pi < M < 0) || (M > pi)
        
         Eold = M - e;
         
    else
        
         Eold = M + e;
         
    end

    err = 1.0e12;

    while err > tol

    Enew = Eold + (M-Eold + e*sin(Eold))/(1-e*cos(Eold));
    
    err = abs(Enew-Eold);
    
    Eold = Enew;

    end

    E = Enew;
    
end


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

