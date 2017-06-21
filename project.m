%project of MGTN duty cycle mac protocol : energy conservation

% /*********** Alejandro Montero ***********/
% /***********   Chetan KC    ***********/

clear all; 
% Get the time when we start computations:
start_time = clock;
% problem constants
P     = 32;            % Payload [byte]
R     = 31.25;         % CC2420 Radio Rate [kbyte/s]
D     = 8;             % number of levels
C     = 5;             % neighbors size (connectivity)
N     = C*D^2;         % number of nodes
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
%%%%%%%%%%%%%%%%%
d = 0:D;
Tw = 100:400;
% number of node in ring d
if d == 0 
    N_d = 1;
else
    N_d = ((2.*d)-1)*C;
end
% Average number of input links
if d == D
    I_d = 0;
elseif d == 0
        I_d = C;
else
    I_d = ((2.*d)+1)/((2.*d)-1);
end
%input frequency (number of packets that enter a node)
if d == 0
    F_I_d = Fs*D^2*C;
else 
    F_I_d = Fs * (((D^2)-(d.^2))/ (2.*d)-1);
end
% output frequecy (number of packers that leaves the packet)
if d == D
    F_d_out = Fs;
else
    F_d_out = Fs * ((((D^2)-(d.^2)) + ((2.*d)-1)) / (2.*d)-1);
end
% background node's traffic frequency
F_B_d = C -((I_d)*F_d_out);

% energy of n node
Ttx = ((Tw/(Tps + Tal))*((Tps + Tal)/2)) + Tack + Tdata;
E_cs = ((Tcs));
E_tx = (Tcs+Tal+Ttx)*F_d_out;
E_rx = (((3/2)* Tps)+ Tack+Tdata)*F_I_d;
E_ovr = (((3/2)*(Ttx/Tw)) * Tps)* F_B_d;
E_n = E_cs + E_tx + E_rx + E_ovr;

% delay of n node
n = 1:N;
%p = (d.^n);
Tcw  = 15*0.62; 

L_n = ((Tw/2)+ (Tcw/2) + Tdata);

fprintf('energy per node = %f \nlatency per node = %f\n\n',...
    E_n, L_n)
% network energy consumption
a1 = Tcs + Tal + ((3/2) * Tps) * (((Tps+Tal)/2) + Tack + Tdata) * F_B_d;
a2 = F_d_out/2;
a3 = ((((Tps+Tal)/2) + Tcs + Tal + Tack + Tdata) * F_d_out)+(((((3/2)*Tps)+Tack+Tdata) * F_I_d) + ((3/4)*Tps*F_B_d));
E_xmac = ((a1./Tw) + (a2.*Tw) + a3);
%[n,E_xmac] = fminunc(@(n)-E_xmac,N)

% End to End delay (e2e) 
B1 = symsum((1/2), i, 1, d.^n);

B2 = (((Tcw/2) + Tdata));
L_xmac = ((B1 * Tw) + B2);

fprintf('Energy_Comsumptio of network = %f \n End_2_End_Delay = %f\n\n ',...
    E_xmac, L_xmac)
% test plot
figure,
plot(E_xmac, L_xmac);
figure,
plot (E_n , L_n);





