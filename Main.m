%-Abstract
%   
%   MAIN is the driver file that computes the thrust pitch/yaw
%   angle profiles, and the final state for the Mars Ascent Vehicle (MAV).   
%   First, it declares the variables associated with the data for  
%   Mars, the launch site, the first stage solid rocket motor (SRM1), and    
%   the second stage solid rocket motor (SRM2). Second, it provides the 
%   user the options of choosing a guidance type and an appropriate 
%   guidance flag. Third, the user is given the flexibility to perform
%   a single-run analysis or a Monte-Carlo (multiple-run) analysis. For a
%   Monte-Carlo run, two primary variables - (1) mass at launch and (2) 
%   stage 2 delta-v can be dispersed. Finally, MAV_ASCENT_TO_CIRCULARIZE 
%   computes the thrust pitch/yaw profiles and the final state(s).
%   PLOT_FIGURES plots the results.
%
%-Inputs
%
%   Guidance_type           the scalar string name of the guidance type.
%                           Options include:
%                           (1) Q-guidance
%                           (2) lvlh
%                           (3) socp
%
%   Guidance_flag           the scalar string name of the guidance flag.
%                           Options include:
%                           (1) Q-flag
%                           (2) box-flag
%                           (3) uhat-flag
%
%   run_count               a double precision scalar representing the 
%                           number of runs. Set to 1 for a single run or
%                           a specific number for Monte-Carlo analysis.
%                           e.g. run_count = 100.
%
%   SIGMA_m                 a double precision scalar representing the
%                           3-sigma limit for the mass at launch (kg).
%
%   SIGMA_Ve_SRM2           a double precision scalar representing the
%                           3-sigma limit for the stage 2 delta-v (percent).
%
%-Output
%
%   For a single run, the pitch and yaw angle thrust profiles are 
%   plotted. For a multiple-run, figures showing variations of (1) final 
%   periapsis altitude (HP) with final apoapsis altitude (HA), and (2) 
%   final semi-major axis altitude (SMA) with final eccentricity (ECC) are
%   plotted.
%
%-&

clear all; close all; clc

% All units are m, s, rad, kg

global R MU W lon0 lat0 grav slope

global m0_SRM1 tburn_SRM1 b_SRM1 Ve_SRM1

global m0_SRM2 tburn_SRM2 b_SRM2 a_SRM2 dV_SRM2 Ve_SRM2

grav = 0; slope = 0; 


% Mars Data

R  = 3389.5e3;                          % Radius of Mars (m)

MU = 4.2828372e13;                      % Gravitational Parameter (m^3/s^2)

W  = ( 868.22e3 / 3600 ) / R;           % Angular Rate (rad/s)


% Launch Site Data

lon0 = 77.50 * pi/180;

lat0 = 18.43 * pi/180;


% SRM1 Data

m0_SRM1     = 393.79;               % First stage Initial mass (kg)   

tburn_SRM1  = 76;                   % SRM1 burn time (s)

b_SRM1      = (215.619/tburn_SRM1); % Propellant mass flow rate (kg/s)

a_SRM1      = m0_SRM1/b_SRM1;       % alpha = m0/b (s)

T_SRM1      = 7925.77;              % SRM1 avg. thrust (N)

Ve_SRM1     = T_SRM1/b_SRM1;        % Constant exhaust velocity (m/s)


% SRM2 Data

m0_SRM2         = 16+1.208+31.838+14.296+54.467;        % Second stage Initial mass (kg)

tburn_SRM2      = 24.45;                                % SRM2 burn time (s)

b_SRM2          = 54.467/tburn_SRM2;                    % Propellant mass flow rate (kg/s)

a_SRM2          = (m0_SRM2)/b_SRM2;                     % alpha = m0/b (s)

dV_SRM2         = 1696;                                 % SRM2 avg. thrust (m/s)

Ve_SRM2         = -dV_SRM2/log(1-tburn_SRM2/a_SRM2);    % Constant exhaust velocity (m/s)




%Guidance_type = 'Q-guidance'; Guidance_flag = 'Q-flag'; 

% Guidance_type = 'lvlh'; Guidance_flag = 'uhat-flag'; 

 Guidance_type = 'socp'; Guidance_flag = 'box-flag'; 

run_count = 5;

SIGMA_m = 1/3;                    % units in kg

SIGMA_Ve_SRM2 = 1/3;              % percent 



[T, THETA, PSI, HA, HP, SMA, ECC, INC, RAAN] = MAV_ascent_to_circularize(Guidance_type, Guidance_flag, run_count, SIGMA_m, SIGMA_Ve_SRM2);

res = plot_figures(T, THETA, PSI, HA, HP, SMA, ECC, run_count);


