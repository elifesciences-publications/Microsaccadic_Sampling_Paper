% This function calculate the light induced current inside each microvillus
% in voltage clamp condition.

% It includes the changes of reversal potential of TRP, but it does not
% include the membrane voltage feedback, neither the Mg2+ blockage of TRP
% channel.

% The input Ntrp is the channel opennings simulated by Gillespie algorithm,
% TRP_Rev_pre is the TRP reversal potential calculated at the last
% simulation step

% The outputs are I_m ICa INa IMg IK Rev, where I_m is the total current
% influx to the microvillus, in unit of pA. ICa is the Calcium influx, INa
% is the sodium current influx, IMg is the meganesium current influx, IK is
% the potassiunm current influx, Rev is the TRP reversal potential at this
% simulation step

% Cai is the intracellular free calcium concentration simulated by
% Gillespie algorithm.

% GHK current equation
% The GHK flux equation for an ion S (Hille 2001):
% 
%     \Phi_{S} = P_{S}z_{S}^2\frac{V_{m}F^{2}}{RT}\frac{[\mbox{S}]_{i} - [\mbox{S}]_{o}\exp(-z_{S}V_{m}F/RT)}{1 - \exp(-z_{S}V_{m}F/RT)} 
% 
% where
% 
%     * ΦS is the flux across the membrane carried by ion S, measured in amperes per centimeter squared (A = C·s?1cm?2)
%     * PS is the permeability of ion S measured in m·s?1
%     * zS is the valence of ion S
%     * Vm is the transmembrane potential in volts 
%     * F is the Faraday constant, equal to 96,485 C·mol?1 or J·V?1·mol?1
%     * R is the gas constant, equal to 8.314 J·K?1·mol?1
%     * T is the absolute temperature, measured in kelvins (= degrees Celsius + 273.15)
%     * [S]i is the intracellular concentration of ion S, measured in mol·m?3 or mmol·l?1
%     * [S]o is the extracellular concentration of ion S, measured in mol·m?3

% The calculation of TRP reversal potential based on GHK equations follow
% the loop of Si,So,Vm,I_influx----->fq----->P1----->Pq=wq*P1----->Iq-----
%              |  |       |<----------------------------------------------|    

% Zhuoyi Song, updated in 06/2017


function [I_m ICa INa IMg IK Rev] = TRP_Rev_GHK_online(Ntrp,Cai,TRP_Rev_pre)
% Initialisation 
g_TRP = 8;

% Parameters 
Vrev = 0.07; % membrane reversal potential 70mV
Cai0 = 0.00016; %mM
Nao = 120; Nai = 8; %mM 
Mgo =4; Mgi = 3;%mM
Ko = 5; Ki = 140; % mM Ionic concentrations from Postma 1999

Vm = -70*0.001; %mV*0.001 = volts
Lm = 1.5; % length of microvilli um
Dm = 0.06; % outer diameter of a microvillus, um
dm = 0.05; % inner diameter of a microvillus, um
A_villi =  pi * Lm* Dm * 10^(-8); % cm^2
V_villi = pi*(dm^2)*Lm/4; % um^3

Zca = 2; Zna = 1; Zmg = 2; Zk =1;
F = 96485; % C*mol^-1
R = 8.314; %J*K^-1*mol^-1
T = 293;   %20oC
Cao = 1.5; %mM

Wca = 0.85; Wmg = 0.11; Wna = 0.02; Wk = 0.02; % the permeability ratios (Hardie and Minke, 1992)

% Ca
beta = F/R/T; 
fCa =  Zca.*beta.*Vm.*((Cai-Cao.*exp(-Zca.*Vm.*beta))./(1-exp(-Zca.*Vm.*beta)));
fNa =  Zna.*beta.*Vm.*((Nai-Nao.*exp(-Zna.*Vm.*beta))./(1-exp(-Zna.*Vm.*beta)));
fMg = Zmg.*beta.*Vm.*((Mgi- Mgo.*exp(-Zmg.*Vm.*beta))./(1-exp(-Zmg.*Vm.*beta)));
fK = Zk.*beta.*Vm.*((Ki- Ko.*exp(-Zk.*Vm.*beta))./(1-exp(-Zk.*Vm.*beta)));

Episilon = A_villi.*F.*(Zca.*Wca.*fCa + Zmg.*Wmg.*fMg + Zna.*Wna.*fNa +Zk.*Wk.*fK).*10^12;% in pA
I_in = (Ntrp+0.0001)*g_TRP*(TRP_Rev_pre + Vrev);%pA
P1 = I_in/Episilon;
Pca = 0.85*P1;
Pmg = 0.11*P1;
Pna = 0.02*P1;
Pk = 0.02*P1;
Rev = R*T*log((Pna*Nao + Pca*Cao + Pmg*Mgo + Pk*Ko)/(Pna*Nai + Pca*Cai + Pmg*Mgi + Pk*Ki))/F; % at rest 6.3mV

ICa = Pca.*fCa.*F.*Zca.*A_villi*10^12; %pA  outward current is defined here as positive
INa = Pna.*fNa.*F.*Zna.* A_villi * 10^12; %pA
IMg = Pmg.*fMg.*F.*Zmg.* A_villi * 10^12; %pA
IK = Pk.*fK.*F.*Zk.* A_villi * 10^12;%pA
I_m = ICa + INa + IMg + IK; 

