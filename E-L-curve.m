%project of MGTN duty cycle mac protocol : energy conservation
% /*********** Alejandro Montero ***********/
% /***********   Chetan KC    ***********/
clear;
clc;
close all;

% Get the time when we start computations:
start_time = clock;
% problem constants
P     = 32;            % Payload [byte]
R     = 31.25;         % CC2420 Radio Rate [kbyte/s]
D     = 8;             % number of levels
C     = 5;             % neighbors size (connectivity)
N     = C*D^5;         % number of nodes
Lmax  = 1000;          % Maximal allowed Delay (ms)
Emax  = 5;             % Maximal Energy Budjet (J)
Fs   = 1/(60*30*1000); % Min traffic rate 1 pkt/half_hour = 1/(60*30*1000) pk/ms

% Parameter Bounds
Tw_max  = 400;         % Maximum Duration of Tw in ms
Tw_min  = 100;         % Minimum Duration of Tw in ms
L_pbl = 4;             % preamble length [byte]
L_hdr = 9 + L_pbl;     % header length [byte]
L_ack = 9 + L_pbl;     % ACK length [byte]
L_ps  = 5 + L_pbl;     % preamble strobe length [byte]
Tal  = 0.95;           % ack listen period [ms]
Thdr = L_hdr/R;        % header transmission duration [ms]
Tack = L_ack/R;        % ACK transmission duration [ms]
Tps  = L_ps/R;         % preamble strobe transmission duration [ms]
Tcw  = 15*0.62;        % Contention window size [ms]
Tcs  = 2.60;           % Time [ms] to turn the radio into TX and probe the channel (carrier sense)
Tdata = Thdr + P/R + Tack; % data packet transmission duration [ms]

Tw = 100:400;
%fprintf('value of Tw = % d\n', Tw)
%Fs = 0.005:5; % to get different value of E_xmac and L_xmac change Fs and
%comment the above Fs and fix the Tw to 200.
% network energy consumption
d = 1;
I_d = ((2*d)+1)/((2*d)-1);
F_I_d = Fs * (((D^2)-(d^2))/(2*d)-1);
F_d_out = Fs * ((((D^2)-(d.^2)) + ((2.*d)-1)) / (2.*d)-1);
F_B_d = (C - I_d)*F_d_out;
a1 = norm(Tcs + Tal + ((3/2) * Tps) * (((Tps+Tal)/2) + Tack + Tdata) * F_B_d);
a2 = norm(F_d_out/2);
a3 = norm(((((Tps+Tal)/2) + Tcs + Tal + Tack + Tdata) * F_d_out)+(((((3/2)*Tps)+Tack+Tdata)...
     * F_I_d) + ((3/4)*Tps*F_B_d)));
E_xmac = ((a1./Tw) + (a2.*Tw) + a3);
%End to End delay (e2e) 

syms i
B1 = symsum((1/2), i, 1, D);   
B2 = symsum(((Tcw/2) + Tdata), i, 1, D);
L_xmac= ((B1 * Tw) + B2);

%printing
fprintf('Energy_Comsumptio E_xmac = %f \n End_2_End_Delay L_xmac = %f\n\n ',...
     E_xmac,L_xmac)
%plot
figure,
plot(E_xmac, L_xmac);
title( ['E-L curve for different Tw']);
xlabel('energy comsumption');
ylabel('End_2_End Delay');