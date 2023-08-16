clc; close all;



%sampling_n = 10;

%% System parameters :

V_bat = 380;             % battery supply voltage [V]
R_i = 22e-3;              % wire resistance [Ohm]

R_i_Bat = R_i;

R_3AC_Load = 100;        % load at 3AC output
R_load_1 = 100;            % load resistance [Ohm]
R_load_2 = 100;

MPPT_Voc = 360;         % MPPT open circuit voltage
MPPT_delta_D = 0.0005;
MPPT_epsilon = 0.01;
MPPTPERIOD = 5.0e-04; 
Kp_Ipv = 1;
Ki_Ipv = 10;
MPPT_reset = 0;
Kp_MPPT_I = 18.75;
Ki_MPPT_I = 165;

PvShortCircuitCurrent = 3; % [A]



V_Bat_C_init = 0;
V_Load_C_init = 0;
V_PV_C_init = 0;
V_AC_a_C_init = 0;
V_AC_b_C_init = 0;
V_AC_c_C_init = 0;


V_DC_init = max(V_bat,MPPT_Voc); % like in real experiment the C will be charged through the forward diodes

V_DC_ref_init = 650;
V_load_1_ref_init = 230; 
V_load_2_ref_init = 96;

%% PARAMETERS

Lb = 3.3e-3;            % Inductance used for the boost [H]
Rb = 22e-3;             % Resistance of the inductor used for the boost [Ohm]

C = 500e-6;             % Capacitance of the DC bus [F]
C_Bat = C;

R_IGBT = 1e-3;          % Resistor of IGBT in on state

L_Load_1 = 2.5e-3;
L_Load_2 = 2.5e-3;

K_Load_1 = 200;
K_Load_1_alpha = 10000;

K_Load_2 = 200;
K_Load_2_alpha = 10000;

K_Bat_alpha = 25;
K_Bat = 30;
K_Bat_dot = 0.1;

Ki_Bat = 2500;
Ki_Bat_alpha = 100;

L_Bat = 2.5e-3;

%% Time parameters:

fs = 20e3;              % Interrupt frequency [Hz]
Ts = 1/fs;              % Interrupt period [s]
fsw = fs;               % Switching frequency [Hz]
Tsw = 1/fsw;            % Switching period [s]
Td = 1*Ts + 0.5*Tsw;    % Control delay [s]

%% Boost currents regulators

K1 = 1/Rb;
T1 = Lb/Rb;
Tn = T1;
Ti = 2*K1*Td;
%Kp = Tn/Ti;
%Ki = 1/Ti;


%% init Voltage Capacitors
%V_C10_init = 20;
%V_C11_init = 10;
%V_C12_init = 10;


C_Load_1 = 10e-3;
C_Load_2 = 10e-3;

C_DC = 0.1*C; 

%% ---------------------------------------------------------------------------- Test
% This file is used to compute the values of Kp and Ki for the current regulator.
% Place this file in the same folder than the Simulink file.


% PARAMETERS
fs = 20e3;              % Sampling, interrupt and switching frequency [Hz]
Ts = 1/fs;              % Sampling period [s]

Lb = 2.5e-3;            % Inductance used for the boost [H]
Rb = 22e-3;             % Resistance of the inductor used for the boost [Ohm]

Vb = 100;               % Battery voltage [V]
Rload = 79;             % Load resistance [Ohm]

Vdc0 = V_DC_ref_init;             % Inital DC bus voltage [V]
Ib0 = 0;                % Inital boost current [A]


%% BOOST CURRENT REGULATOR -------------------------------------------------

K1 = 1/Rb;
T1 = Lb/Rb;

Td1 = (1/2 + 1/2 + 1/3)*Ts; % 1/2Ts of sampling, 1/2Ts of computation, 1/3 for the triangle carrier

Tn1 = T1;
Ti1 = 2*K1*Td1;

Kp1 = Tn1/Ti1;
Ki1 = 1/Ti1;

% new PI Version:

% Batterie Boost cascaded PI:
%Kp_Bat_V = 0.9375;
%Ki_Bat_V = 439.45;
%Kp_Bat_I = 1.875;
%Ki_Bat_I = 1652;

% Load Buck cascaded PI:
%Kp_Load_V = 0.01;
%Ki_Load_V = 10;
%Kp_Load_I = 0.01;
%Ki_Load_I = 10;

% Batterie Boost cascaded PI:
Kp_Bat_V = 5; %0.4875;
Ki_Bat_V = 300; %228.5156;
Kp_Bat_I = 3; %18.75;
Ki_Bat_I = 300; %165;

% Load Buck cascaded PI:
Kp_Load_V = 0.1; %0.01;
Ki_Load_V = 2; %10;
Kp_Load_I = 0.001 %0   .01;
Ki_Load_I = 1; %10;

Kp_Load_2_V = Kp_Load_V; %0.01;
Ki_Load_2_V = Ki_Load_V; %10;
Kp_Load_2_I = Kp_Load_I %0   .01;
Ki_Load_2_I = Ki_Load_I; %10;

% BOOST CASCADED VOLTAGE REGULATOR ----------------------------------------
T2 = C;
Tdeq = 2*Td1;

a = 4;
Tn2 = a^2*Tdeq;
Ti2 = a^3*Tdeq^2/T2;

Kp2 = Tn2/Ti2;
Ki2 = 1/Ti2;


%% PV

PV_cell_series_n = 8;           % each of the cells has Voc=10V
PV_cell_parallel_n = 10;  
V_oc_per_cell = 10;

%%  3AC

SwitchingFreq = 20e3;           % Switching frequency [Hz]
Deadtime = 1e-6;              % Deadtime [s]

ControlFreq = SwitchingFreq;    % Control frequency [Hz]
ControlPeriod = 1/ControlFreq;  % Control period [s]


GridFreq = 50;                  % Grid frequency [Hz]
%Vg = 12;                       % RMS Grid phase voltage [V]
Lg = 2.5e-3;                    % Grid-side inductor value [H]
Rg = 22e-3;                     % Series equivalent resistance of Lg [ohm]


Tdtot = ControlPeriod;          % Total control loop delay
                                % See PN142 for more details

%Kp_dq = Lg/(2*Tdtot);              % Proportional gain of the controller
%Ki_dq = Rg/(2*Tdtot);              % Integral gain of the controller

Kp_dq = 10;              % Proportional gain of the controller
Ki_dq = 500;              % Integral gain of the controller


%%% Values for nonlinear

R_l = 0.1;
L1 = 2.5e-3;
%Supercap
C1= 500e-6;
C2= 500e-6;
R1=R_i;
R2=0.1;
R3=R_i;

R01=1e-3;
R03=R01;

L3 = 2.5e-3;

K3=3000;
K31=50000;

c2= 500;
c3= 10000;

%Battery
C4=0.01;
C5=0.01;
R4=0.1;
R5=0.1;
R04=0.01;
L6=0.0033;

K6=5/L6;
K61=10/L6;

%PV
C7=0.01;
C8=0.01;
R7=0.1;
R8=0.1;
R07=0.01;
R08=0.01;
L9=0.0033;

K9=10/L9;
K91=10/L9;

%Load
C11=0.01;
C12=0.01;
R11=0.1;
R12=0.1;
R011=0.01;
R012=0.01;
L13=0.0033;

K11=10/R11;
K111=0;

K13=10/L13;
K131=50/L13;

K10 = 1;
K101 = 10;

%Supercap
C1=0.01;
C2=0.01;
R1=0.1;
R2=0.1;
R01=0.01;
L3=0.0033;

K3=3000;
K31=50000;

c2= 500;
c3= 10000;

%Battery
C4=0.01;
C5=0.01;
R4=0.1;
R5=0.1;
R04=0.01;
L6=0.0033;

K6=5/L6;
K61=10/L6;

%PV
C7=0.01;
C8=0.01;
R7=0.1;
R8=0.1;
R07=0.01;
R08=0.01;
L9=0.0033;

K9=10/L9;
K91=10/L9;

%Load
C11=0.01;
C12=0.01;
R11=0.1;
R12=0.1;
R011=0.01;
R012=0.01;
L13=0.0033;

K11=10/R11;
K111=0;

K13=10/L13;
%K131=0;

%Train
C14=0.01;
C15=0.01;
R14=0.1;
R15=0.1;
R014=0.01;
R015=0.01;
L16=0.0033;

K14=1050;
K141=50000;

K16=1200;
K161=200000;

%ACgrid
Pnom=100e3;
Vnom=400;

R17=0.1;
C17=500e-3;

R=0.002;
L=500e-6;

Kd=5000/L;
Kd1=50000/L;

Kq=5000/L;
Kq1=50000/L;

C10=0.0001;
%fs=20000;


%Omer model
C = 500e-6;
R_i = 0.1;
R_l = 0.1;
L = 2.5e-3;
R_on = 0.01;

%battery
C2=C;
R1=R_i;
R2=R_l;
R01=R_on;
L3=L;

K3=3000;
K31=50000;
c2= 500;
c3= 10000;

%PV
C8=C;
R7=R_i;
R8=R_l;
R07=R_on;
R08=R_on;
L9=L;

K9=10/0.0033;
K91=10/0.0033;

%Load
C11=0.01; %Inserted in the Buck converter of the Load
C12=20*C; %This one is different from Omer model
R11=R_i;
R12=R_l;
R011=R_on;
R012=R_on;
L13=L;

K11=1/0.1;
K111=2/0.1;
K13=1/0.0033;
K131=2/0.0033;

C10=C/2; %capacitor inserted in the DC bus
fs=20000; %frequency switching


% PI calculated Version:
%Bat
%Kp_Bat_PI_1 = 18.75;
%Ki_Bat_PI_1 = 165;
%Kp_Bat_PI_2 = 0.9375;
%Ki_Bat_PI_2 = 439.4531;

%Kp_Bat_PI_1 = 0.9375;
%Ki_Bat_PI_1 = 439.4531;
%Kp_Bat_PI_2 = 18.75;
%Ki_Bat_PI_2 = 165;

Kp_Bat_PI_1 = 0.4;
Ki_Bat_PI_1 = 20;
Kp_Bat_PI_2 = 2.5;
Ki_Bat_PI_2 = 625;

Kp_Bat2_PI_1 = 0.4;
Ki_Bat2_PI_1 = 20;
Kp_Bat2_PI_2 = 2.5;
Ki_Bat2_PI_2 = 625;

%L1
Kp_Load_PI_1 = 0.1;
Ki_Load_PI_1 = 5;
Kp_Load_PI_2 = 0.001;
Ki_Load_PI_2 = 0.1;
Kp_Load2_PI_1 = Kp_Load_PI_1;
Ki_Load2_PI_1 = Ki_Load_PI_1;
Kp_Load2_PI_2 = Kp_Load_PI_2;
Ki_Load2_PI_2 = Ki_Load_PI_2;


%adaptive virtual inertia Filipe Parameter:
Dp = 50;
Kw = 20;
M = 2;
fn = 50;
Vg = 230;
