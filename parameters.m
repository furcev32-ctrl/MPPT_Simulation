% This file is used to compute the values of Kp and Ki for the current and
% voltage regulators. Place it in the same folder as the Simulink file.


% SYSTEM PARAMETERS -------------------------------------------------------
CTRLPERIOD = 1/20e3;        % Control is at 20kHz
MPPTPERIOD = 1/2e3;         % MPP tracking is at 2 kHz

Ts = 1.5*CTRLPERIOD;        % Control delay

Lb = 2.5e-3;                % Inductor for the boost
Rb = 22e-3;                 % Resistance of the inductor used for the boost
Lg = 2.5e-3;                % Grid's inductors
Rg = 22e-3;                 % Resistance of grid's inductors


Isc = 50;                   % Short-circuit current (PV panel)
Voc = 300;                  % Open-circuit voltage (PV panel)

% INITIALIZATION VALUES ---------------------------------------------------
Vdc0 = 650;                 % Initial DC bus voltage


% BOOST CURRENT REGULATOR -------------------------------------------------
K1 = 1/Rb;
T1 = Lb/Rb;

Tn1 = T1;
Ti1 = 2*K1*Ts;

Kp1 = Tn1/Ti1;
Ki1 = 1/Ti1;

% BOOST CASCADED VOLTAGE REGULATOR ----------------------------------------
T2 = Cdc;
Ts2 = 13*CTRLPERIOD;

Tn2 = 4*Ts2;
Ti2 = 8*Ts2^2/T2;

Kp2 = Tn2/Ti2;
Ki2 = 1/Ti2;

% INVERTER CURRENT REGULATOR ----------------------------------------------
Ts3 = Ts;

K3 = 1/Rg;
T3 = Lg/Rg;

Tn3 = T3;
Ti3 = 2*K3*Ts3;

Kp3 = Tn3/Ti3;
Ki3 = 1/Ti3;