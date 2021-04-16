%-Description
%
%   BOX_FLAG implements the box-flag trigger. The trigger is activated the
%   first time the state enters the semi-major axis - eccentricity box.
%
%-Inputs
%   
%   s0              current state 
%
%   sf              target state
%
%   tnow            current time (s)
%
%   Ve              SRM2 exhaust velocity (m/s)
%
%   a               ratio (m0/b) (s)
% 
%   R               radius of Mars (m)
%
%   MU              gravitaional parameter (m^3/s^2)
%
%   tburn_SRM2      SRM2 burn time (s)
%
%   TN2P            Inertial to in-plane transformation matrix
%
%   Guidance_type   the scalar string name of the guidance type
%
%-Outputs
%
%   flag            scalar that stores either 0 or 1. 
%                   flag = 0 and flag = 1 corresponds to trigger being 
%                   inactive and active, respectively. 
%-&


function flag = box_flag(s0,sf,tnow,Ve,a,R,MU,tburn_SRM2,TN2P, Guidance_type)

global sFINAL

sFINAL = zeros(1,6);

options = odeset('RelTol',1e-10);

tstart = tnow;

tdone = tnow + tburn_SRM2;

tf = tdone - tnow;

dt = 1;


x0 = s0(1); y0 = s0(2); z0 = s0(3); 

u0 = s0(4); v0 = s0(5); w0 = s0(6);

[r0,v0] = lvlh2rv(x0,y0,z0,u0,v0,w0,R);

[~,~,~,~,~,nu] = rv2orbel(r0,v0,MU);

nu = nu*180/pi;

if nu >= 177 && nu < 180

    while tf > 0
        
        if strcmp(Guidance_type, 'Q-guidance')
                
           [th,psi] = Q_guidance(s0,sf,tf,Ve,a,R,MU);
                
        elseif strcmp(Guidance_type, 'lvlh') || strcmp(Guidance_type, 'socp')
                
           [th,psi] = lvlh_guidance(s0,sf,tf,Ve,a,R,MU);
                
        else
                
           error('Error: Please input "Q-guidance" or "lvlh" or "socp"')
                
        end
        
        tnext = tnow + min( [dt,tdone-tnow] );
        
        [~,s] = ode45(@lvlh_dynamics,[tnow,tnext],s0,options,R,MU,Ve,a,th,psi,tstart);
        
        s0 = s(end,:);
        
        sFINAL = s0;
        
        tnow = tnext;
        
        tf = tdone - tnow;
        
    end
    
    x = sFINAL(1); y = sFINAL(2); z = sFINAL(3);
    
    u = sFINAL(4); v = sFINAL(5); w = sFINAL(6);
    
    [rf(:,1),vf(:,1)] = lvlh2rv(x,y,z,u,v,w,R);

    % Convert to inertial frame
    
    RF(:,1) = inv(TN2P)*rf(:,1);
    
    VF(:,1) = inv(TN2P)*vf(:,1);

    % Convert to orbital elements
    
    [aF,eF,~,~,~,~] = rv2orbel(RF(:,1),VF(:,1),MU);
    
    SMA = aF - R;
    
    if (SMA <= 3.43e5) && (eF <= 0.0005)
        
        flag = 1;
        
    else
        
        flag = 0;
        
    end
    
elseif nu >= 180
    
    flag = 1;
    
else
    
    flag = 0;
    
end

end


% Support Functions

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