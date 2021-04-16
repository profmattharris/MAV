%-Description
%
%   SOCP_GUIDANCE computes the thrust pitch and yaw angle using SOCP-guidance
%   scheme. 
%
%-Inputs
%   
%   s0          current state 
%
%   sf          target state
%
%   tnow        current time (s)
%
%   tdone       cut-off time (s)
%
%   tstart      ignition time (s)
%
%   Ve          SRM2 exhaust velocity (m/s)
%
%   a           ratio (m0/b) (s)
% 
%   R           radius of Mars (m)
%
%   MU          gravitaional parameter (m^3/s^2)
%
%   g_step      guidance step, defines the size of control variables "th" and
%               "psi".
%
%-Outputs
%
%   th      thrust pitch angle (rad)
%
%   psi     thrust yaw angle (rad)
%
%-Assumption 
%
%   the gravity term is approximated with a linearly time-varying estimate "gnow".
%
%-&

function [th,psi] = socp_guidance(s0,sf,tnow,tdone,tstart,Ve,a,R,MU,g_step)

global grav; global slope;

x0 = s0(1); y0 = s0(2); z0 = s0(3); u0 = s0(4); v0 = s0(5); w0 = s0(6);

xf = sf(1); yf = sf(2); zf = sf(3); uf = sf(4); vf = sf(5); wf = sf(6);


tf = tdone - tnow;

tf_i = round(tf);

tf_f = tf - tf_i;

dt_i = 1;

dt_f = tf_f;

if (tf_i == 24)
    
    grav = -MU/( R+y0 )^2 + (u0^2)/(R+y0);                  
    
    slope = (-MU/( R+yf )^2 + (uf^2)/(R+yf) - grav)/tf;
    
    gnow = grav;
    
else
    
    gnow = grav + slope*(tnow - tstart);
    
end

c = [0; gnow; 0];



if tf_i == 0
    
    N = 2;
    
else
    
    N = tf_i/dt_i + 2;
    
end


s0 = s0';

sf = sf';

A = [zeros(3,3), eye(3,3); zeros(3,6)];

B = [zeros(3,3); eye(3,3)];

% Discretize the systems

sysc = ss(A,B,eye(6),0);

sysd = c2d(sysc,dt_i);

AD_i = sysd.A;

BD_i = sysd.B;

CD_i = sysd.C;

DD_i = sysd.D;

sysd = c2d(sysc,dt_f);

AD_f = sysd.A;

BD_f = sysd.B;

CD_f = sysd.C;

DD_f = sysd.D;


% YALMIP setup
opts            = sdpsettings;

opts.solver     = 'gurobi';

opts.verbose    = 1;

yalmip('clear')

s           = sdpvar(6,N);

u           = sdpvar(3,N-1);

g           = sdpvar(1,N-1);


     LMI         = [s(:,1) == s0];
  
for k = 1:N-1
    
    t   = tnow + (k-1)*dt_i;

    if (k < N-1)
        
        AD  = AD_i; 
        
        BD  = BD_i;
        
    else
        
        AD  = AD_f; 
        
        BD  = BD_f;
        
    end

    tau         = -Ve/(t-tstart-a);
    
    LMI         = [LMI, s(:,k+1) == AD*s(:,k) + BD*u(:,k)+ BD*c];
    
    LMI         = [LMI, norm( u(:,k) ) <= tau];
    
end

    obj = (s(2,N)/yf-1)^2 + (s(3,N)-zf)^2 + (s(4,N)/uf-1)^2 + (s(5,N)-vf)^2+ (s(6,N)-wf)^2; 
    
    sol = solvesdp(LMI,obj,opts);
    
 
if (N == 2)
    
    g_step = 1;
    
end

for i = 1:g_step    
    
    tau1 = double(u(1,i));
    
    tau2 = double(u(2,i));
    
    tau3 = double(u(3,i));
    
    tau_mag = sqrt(tau1^2 + tau2^2+ tau3^2);
    
    psi(i) = atan2(tau3,tau1);
    
    th(i)  = asin(tau2/tau_mag); 
    
end


end






