%-Abstract
%   
%   MAV_ASCENT_TO_CIRCULARIZE computes the thrust pitch and yaw angles
%   corresponding to the second burn and also returns the final state. 
%   The MAV is launched from the Jezero Crater using an open-loop guidance 
%   scheme (the thrust pitch and yaw angles are assumed constant for this 
%   launch phase). This is followed by a coasting or pre-burn phase where 
%   closed-loop guidance calculations targeting a desired final point is  
%   performed. During this phase, an appropriate Guidance_flag is used to 
%   ignite the second stage. Once ignition occurs (burn phase), the thrust 
%   pitch and yaw angles are computed based on a chosen guidance scheme to 
%   perform the continuous thrust manuever.
% 
%-Inputs
%
%   Guidance_type           the scalar string name of the guidance type.
%
%   Guidance_flag           the scalar string name of the guidance flag.
%
%   N                       a double precision scalar representing the 
%                           number of runs. 
%
%   SIGMA_m                 a double precision scalar representing the
%                           3-sigma limit for the mass at launch.
%
%   SIGMA_Ve_SRM2           a double precision scalar representing the
%                           3-sigma limit for the stage 2 delta-v.
%
%-Outputs
%
%   T                       a double precision array representing time
%                           post first burn (s)  
%
%   THETA                   a double precision array respresenting the 
%                           thrust pitch angle corresponding to T (rad)
%                           
%   PSI                     a double precision array representing the
%                           thrust yaw angle corresponding to T (rad)
%                           
%   HA                      a double precision scalar or array representing
%                           the final apoapse altitude (km)
%
%   HP                      a double precision scalar or array representing
%                           the final periapse altitude (km)
%
%   SMA                     a double precision scalar or array representing
%                           the final semi-major axis (km)
%
%   ECC                     a double precision scalar or array representing
%                           the final eccentricity
%
%   INC                     a double precision scalar or array representing
%                           the final inclination (deg)
%
%   RAAN                    a double precision scalar or array representing
%                           the final right ascension of the ascending node
%                           (deg)
%-&


function [T, THETA, PSI, HA, HP, SMA, ECC, INC, RAAN] = MAV_ascent_to_circularize(Guidance_type, Guidance_flag, N, SIGMA_m, SIGMA_Ve_SRM2)


global R MU W lon0 lat0 grav slope

global m0_SRM1 tburn_SRM1 b_SRM1 Ve_SRM1

global m0_SRM2 tburn_SRM2 b_SRM2 a_SRM2 dV_SRM2 Ve_SRM2

grav = 0; slope = 0; 


for k = 1:N
    
    clearvars -except HA HP SMA INC RAAN ECC R MU W lon0 lat0 ...
                k S Guidance_flag Guidance_type SIGMA_m SIGMA_Ve_SRM2 ...
                m0_SRM1 tburn_SRM1 b_SRM1 Ve_SRM1 m0_SRM2 tburn_SRM2 ...
                b_SRM2 a_SRM2 dV_SRM2 Ve_SRM2

    % SRM1 DISPERSIONS
    
    delta_m     = SIGMA_m*1/3*randn(1);
    
    a_act_SRM1  = (m0_SRM1 + delta_m)/b_SRM1;
 
    % SRM2 DISPERSIONS
    
    delta_Ve_SRM2   = SIGMA_Ve_SRM2*1/3*randn(1)/100;
    
    a_act_SRM2      = (m0_SRM2 + delta_m)/b_SRM2;
    
    dV_act_SRM2     = dV_SRM2*(1 + delta_Ve_SRM2);
    
    Ve_act_SRM2     = -dV_act_SRM2/log(1-tburn_SRM2/a_act_SRM2);

    %---------------------------% -----------------------------------------
    % LAUNCH INITIAL CONDITIONS %
    %---------------------------%

    % Initial Conditions in Inertial Frame
    
    ROT = LVLH2N(lon0,lat0);
    
    r0_ascent = ROT*[R;0;0];
    
    v0_ascent = cross( [0;0;W] , r0_ascent );

    % Initial Conditions in LVLH Frame
    
    [x0,y0,z0,u0,v0,w0] = rv2lvlh(r0_ascent,v0_ascent,R);

    %-------------------%--------------------------------------------------
    % ASCENT SIMULATION %
    %-------------------%
    
    s0 = [x0; y0; z0; u0; v0; w0];
    
    T_ascent  = [];
    
    S_ascent  = [];

    for i = 0:1:tburn_SRM1-1
        
        [th,psi] = launch_guidance(s0,s0);
        
        [t,s]    = ode45(@lvlh_dynamics,[i,i+1],s0,[],R,MU,Ve_SRM1,a_act_SRM1,th,psi,0);
        
        T_ascent = [T_ascent; t];
        
        S_ascent = [S_ascent; s];
        
        s0       = s(end,:);
        
    end

    % Final state after stage 1
    
    xf_ascent = s0(1); yf_ascent = s0(2); zf_ascent = s0(3);
    
    uf_ascent = s0(4); vf_ascent = s0(5); wf_ascent = s0(6);

    % Convert the final state to the Inertial Frame
    
    [rf,vf] = lvlh2rv(xf_ascent, yf_ascent, zf_ascent, uf_ascent, vf_ascent, wf_ascent,R);

    % Convert the final state to orbital elements
    
    [~,~,~,Wf,wf,~] = rv2orbel(rf, vf, MU);
    
    % Target state in Inertial Frame
    
    [Rd,Vd] = orbel2rv(R+343e3, 0, 25*pi/180, Wf, wf, pi, MU);

    % Inertial to In-plane transformation
    
    TN2P = N2P(Rd, Vd);

    % In-plane initial and target states
    
    r0 = TN2P*rf;
    
    v0 = TN2P*vf;
    
    rd = TN2P*Rd;
    
    vd = TN2P*Vd;

    % LVLH initial and target states
    
    [x0,y0,z0,u0,v0,w0] = rv2lvlh(r0,v0,R);
    
    [xd,yd,zd,ud,vd,wd] = rv2lvlh(rd,vd,R);

    %-----------------------------------------% ---------------------------
    % COASTING AND CIRCULARIZATION SIMULATION %
    %-----------------------------------------%
    
    dt     = 1;
    
    tnow   = 0-dt;
    
    tnext  = tnow+dt;
    
    tstart = inf;
    
    tdone  = inf;
    
    tfinal = 600;

    s0 = [x0; y0; z0; u0; v0; w0]';
    
    sd = [xd; yd; zd; ud; vd; wd]';
    
    T  = []; 
    
    S  = []; 
    
    THETA = []; PSI = [];

    options = odeset('RelTol',1e-10);


    while tnext < tfinal

        tnow  = tnext;

        % Navigation Model
        
        s0hat = nav_model(s0);
        
        % Pre-burn Phase

        if tnow < tstart 
            
            if strcmp(Guidance_flag, 'Q-flag')
                
                flag_guid = Q_flag(s0,sd,R,MU,dV_SRM2);
                
            elseif strcmp(Guidance_flag, 'box-flag')
                
                flag_guid = box_flag(s0,sd,tnow,Ve_SRM2,a_SRM2,R,MU,tburn_SRM2,TN2P, Guidance_type);
                
            elseif strcmp(Guidance_flag, 'uhat-flag')
                
                flag_guid = uhat_flag(s0,sd,tburn_SRM2,Ve_SRM2,a_SRM2,R,MU);
                
            else
                error('Error: Please input "Q-flag" or "box-flag" or "uhat-flag"')
            end 

            VE = 0; th = 0; psi = 0; tnext = tnow+dt;
            
            if flag_guid == 1

                tstart = tnow; 

                tdone  = tnow + tburn_SRM2;

            end
            
            j = 1;
            
        end
        
        % Burn Phase

        if tnow>=tstart && tnow<tdone 

            if strcmp(Guidance_type, 'Q-guidance')
                
                [th_g, psi_g] = Q_guidance(s0hat,sd,tdone-tnow,Ve_SRM2,a_SRM2,R,MU);
                
                j = 1;
                
            elseif strcmp(Guidance_type, 'lvlh')
                
                [th_g, psi_g] = lvlh_guidance(s0hat,sd,tdone-tnow,Ve_SRM2,a_SRM2,R,MU);
                
                j = 1;
                
            elseif strcmp(Guidance_type, 'socp')
                
                num = 3;        % defines the desired number of calls to socp_guidance
                
                g_step = fix(tburn_SRM2/num);
                
                tgo_int = fix(tdone-tnow);
                
                if (mod(tgo_int,g_step) == 0) && (tgo_int > 0)
                
                    if (tgo_int == g_step)
                        g_step = g_step + 1;
                    end

                    [th_g, psi_g] = socp_guidance(s0hat,sd,tnow,tdone,tstart,Ve_SRM2,a_SRM2,R,MU,g_step);

                    j = 1;
                end
                
            else
                
                error('Error: Please input "Q-guidance" or "lvlh" or "socp"')
                
            end 

            % Control Model
            
            [th,psi,VE] = ctrl_model(th_g(j), psi_g(j), Ve_act_SRM2);
            
            
            tnext = tnow + min( [dt,tdone-tnow] );
            
        end
        
        % Post-burn Phase

        if tnow >= tdone 
            
            VE = 0; th = 0; psi = 0; tnext = tfinal;
            
        end


        [t,s] = ode45(@lvlh_dynamics,[tnow,tnext],s0,options,R,MU,VE,a_act_SRM2,th,psi,tstart);
        
        s0 = s(end,:);

        THETA = [THETA; repmat(th,length(t),1)];
        
        PSI   = [PSI; repmat(psi,length(t),1)];
        
        T   = [T; t];
        
        S   = [S; s];
        
        j   = j + 1;

    end


    %-----------------% ---------------------------------------------------
    % POST-PROCESSING %
    %-----------------%

    % Convert final point to in-plane 
    
    for i = 1:length(T)
        
        x(i) = S(i,1); y(i) = S(i,2); z(i) = S(i,3);
        
        u(i) = S(i,4); v(i) = S(i,5); w(i) = S(i,6);
        
        [rf(:,i),vf(:,i)] = lvlh2rv(x(i),y(i),z(i),u(i),v(i),w(i),R);

        % Convert to inertial frame
        
        RF(:,i) = inv(TN2P)*rf(:,i);
        
        VF(:,i) = inv(TN2P)*vf(:,i);

        % Convert to orbital elements
        
        [aF(i),eF(i),iF(i),WF(i),wF(i),nuF(i)] = rv2orbel(RF(:,i),VF(:,i),MU);
        
    end

    nuF = nuF*180/pi; nuF(nuF<100) = nuF(nuF<100)+180;

    HA(k) = (aF(end)*(1+eF(end))-R)/1000;
    
    HP(k) = (aF(end)*(1-eF(end))-R)/1000;

    SMA(k) = (aF(end))/1000;
    
    ECC(k) = eF(end);
    
    INC(k) = iF(end)*180/pi;
    
    RAAN(k) = WF(end)*180/pi;

end

end



%-------------------% -----------------------------------------------------
% SUPPORT FUNCTIONS %
%-------------------%

function [th,psi] = launch_guidance(~,~)

th  = (42.0123)*pi/180;

psi = (19.4388)*pi/180;

end

function s = nav_model(s)

s = s + 0*randn(1,6);

end

function [th,psi,VE] = ctrl_model(th,psi,Ve)

th  = th;

psi = psi;

VE  = Ve;

end


function [r,v] = lvlh2rv(x,y,z,u,v,w,R)

lon = x/R;

lat = z/R;

TT = LVLH2N(lon,lat);

r = TT*[R+y;0;0];

v = TT*[v;u;w];

end


function [x,y,z,u,v,w] = rv2lvlh(r0,v0,R)

rhat = r0 / norm(r0);

lon = atan2( rhat(2), rhat(1) );

lat = asin( rhat(3) );

x = R*lon;

y = norm(r0) - R;

z = R*lat;

TT = LVLH2N(lon,lat);

dum = inv(TT)*v0;

v = dum(1);

u = dum(2);

w = dum(3);

end


function TT = LVLH2N(lon,lat)

TT = [cos(lon)*cos(lat), -sin(lon), -cos(lon)*sin(lat);
     sin(lon)*cos(lat),  cos(lon), -sin(lon)*sin(lat);
     sin(lat),           0,         cos(lat)];

end


function TN2P = N2P(r,v)

% Inertial directions

ix = [1;0;0];

iy = [0;1;0];

iz = [0;0;1];

% In-plane directions

i1 = r/norm(r);

i3 = cross(r,v) / norm( cross(r,v) );

i2 = cross(i3,i1);

% Rotation matrix from Inertial to In-plane

TN2P = [dot(ix,i1), dot(iy,i1), dot(iz,i1);
        dot(ix,i2), dot(iy,i2), dot(iz,i2);
        dot(ix,i3), dot(iy,i3), dot(iz,i3)];
    
end
