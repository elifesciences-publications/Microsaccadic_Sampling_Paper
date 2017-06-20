% This function calculated the net current of NaCa exchanger and Cacium pump in 
% the cell body. 

% The input Vol is the cell membrane voltage in volts, Cai is the free
% calcium concentration inside the cell body, in unit of mM, Nai is the free
% socium concentration inside the cell body, in unit of mM.
% The output Inaca is the net current of Exchanger, in unit of pA, inward
% current as positive, Ica is the net current of calcium pump.

% Zhuoyi Song, updated in 06/2017

%Reference: Juusola et al 2017 ELife
function [Inaca Ica] = NaCaPump_Body(Vol,Cai,Nai)
F=96485;
R=8.314472;
T=293;
nao=120;
nai=8;
cao=1.5;
gama=0.65;
nnaca = 3; %here if 0.004 is used, F/2/R/T is used, the feedback is 1/10 of cai, or if the cai is amplified 4 times, the feedback is 1/40
Sb = 1.57*10^-5; % cell body membrane area;
Kmcai = 0.5;

vo = Vol*0.001;
beta=0.0036;
knaca_pasi=0.001;

dnaca_pasi = 0.001;
Inaca_nu = exp(gama.*(nnaca-2).*vo.*F./2./R./T).*(Nai.^3).*cao-exp(-(1-gama).*(nnaca-2).*vo.*F./2./R./T).*(nao.^3).*Cai;
Inaca_de = (87.5^3+nao^3)*(1.38+cao)*(1+dnaca_pasi*exp((gama-1)*(nnaca-2).*vo.*F./2./R./T));%the dominator have a phenomenon of slow tail
Inaca = -1200*Inaca_nu./Inaca_de;

Aca = 2000*4*Sb; %nA
Ica = Aca .*(1./((Kmcai./(cao-Cai)).^1+1));
