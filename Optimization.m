%project of MGTN duty cycle mac protocol : energy conservation

% /*********** Alejandro Montero ***********/
% /***********   Chetan KC    ***********/
clear;
clc;
close all;
% problem constants
P     = 32;            % Payload [byte]
R     = 31.25;         % CC2420 Radio Rate [kbyte/s]
D     = 8;             % number of levels
C     = 5;             % neighbors size (connectivity)
N     = C*D^2;         % number of nodes
Lmax  = 0.1:5;          % Maximal allowed Delay (ms)
Emax  = 0.5:5;             % Maximal Energy Budjet (J)
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
Tw = 200;
%Fs = 0.005:5;
%network energy consumption
d = 1;
% calculated for question 1.
I_d = ((2*d)+1)/((2*d)-1);
F_I_d = Fs * (((D^2)-(d^2))/(2*d)-1);
F_d_out = Fs * ((((D^2)-(d.^2)) + ((2.*d)-1)) / (2.*d)-1);
F_B_d = (C - I_d)*F_d_out;
a1 = (Tcs + Tal + ((3/2) * Tps) * (((Tps+Tal)/2) + Tack + Tdata) * F_B_d);
a2 = norm(F_d_out/2);
a3 = ((((Tps+Tal)/2) + Tcs + Tal + Tack + Tdata) * F_d_out)+(((((3/2)*Tps)...
    +Tack+Tdata)* F_I_d) + ((3/4)*Tps*F_B_d));

E_xmac = max(((a1/Tw) + (a2*Tw) + a3),N);

%End to End delay (e2e) 
syms i
b1 = symsum((1/2), i, 1, D);   
b2 = symsum(((Tcw/2) + Tdata), i, 1, D);

L_xmac = max(((b1*Tw) + b2),N);

%calculations formulas
T_tx = ((Tw/(Tps + Tal))*((Tps + Tal)/2)) + Tack + Tdata;
E_tx = (Tcs+Tal+T_tx) * F_d_out;

%fprintf('B1 is %f \n',B1)
%minimization of energy
%n = size(100,400);
cvx_begin
  variable E_xmac(Tw);
  variable L_xmac(Tw);
  minimize(norm(E_xmac(Tw)));
  subject to 
  L_xmac(Tw) <= Lmax;
  Tw >= Tw_min;
  norm(C)* E_tx <= 1/4;
cvx_end

%minimization of delay
cvx_begin
  variable E_xmac(Tw);
  variable L_xmac(Tw);
  minimize(norm(L_xmac(Tw)));
  subject to 
  E_xmac(Tw) <= Emax;
  Tw >= Tw_min;
  norm(C)* E_tx <= 1/4;
cvx_end

%printing
%fprintf('E_xmac = %f \n L_xmac = %f \n\n ', E_xmac,L_xmac)

%plot
figure,
plot(norm(E_xmac), norm(L_xmac));
%title( ['E-L curve for different Tw']);
%xlabel('energy comsumption');
%ylabel('End_2_End Delay');

