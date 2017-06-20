% This function calculated the net current of NaK pump in the whole cell
% NaK pump is mainly located in the cell body.

% The input Vol is the cell membrane voltage in volts, Nai is the free
% sodium concentration inside the cell body, in unit of mM.
% The output Inak is the net current of NaK pump, in unit of pA, inward
% current as positive

% Zhuoyi Song, updated in 06/2017

%Reference: Juusola et al 2017 ELife
function Inak = NaKPump(Vol,Nai)
vo = Vol*0.001;
F=96485;
R=8.314472;
T=293;
nao=120;
Ko = 5;
Kp = 33; % Kp  is the intracellular Na concentration for half-maximal Na+ current of the pump
A = -3.7/3/5.1; % A is the maximal Na+ current for the pump
Sb = 1.57*10^-5; % cell body membrane area;
Inak = A .* Nai./(Nai+Kp);




