% This function calculated the net current of NaCa exchanger in each
% microvillus. It is assumed here that NaCa exchanger is located inside
% each microvillus and colocalized with TRP.

% The input Vol is the cell membrane voltage in volts, Cai is the free
% calcium concentration inside the microvillus, in unit of mM.
% The output Inaca is the net current of Exchanger, in unit of pA, inward
% current as positive

% Zhuoyi Song, updated in 06/2017

%Reference: Juusola et al 2017 ELife
function Inaca = NaCaPump(Vol,Cai)

F=96485;
R=8.314472;
T=293;
nao=120;
nai=8;
cao=1.5;
gama=0.65;
nnaca = 3;

beta=0.0036;
knaca_pasi=0.001;
dnaca_pasi = 0.1;
Inaca_nu = exp(gama.*(nnaca-2).*Vol.*F./2./R./T).*(nai.^3).*cao-exp(-(1-gama).*(nnaca-2).*Vol.*F./2./R./T).*(nao.^3).*Cai;
Inaca_de = (87.5^3+nao^3)*(1.38+cao)*(1+dnaca_pasi*exp((gama-1)*(nnaca-2).*Vol.*F./2./R./T));%the dominator have a phenomenon of slow tail
Inaca = 0.000004*Inaca_nu./Inaca_de; % in unit of nA, outward current as positive, this is an inward current, so the current is negative
Inaca = -1000*Inaca; % in unit of pA, since in my simulations, inward current is defined as positive, hence there is a minus
