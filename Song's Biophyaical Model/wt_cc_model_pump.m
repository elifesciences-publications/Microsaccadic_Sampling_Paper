% This is the programme for simulating voltage response of Drosophila
% photoreceptor offline. It includes the shab, shaker and new K
% conductances.

% Modified and streamlined function based in Mikko Vahasoyrinki's model
% Zhuoyi Song, updated in 06/2017

function  [y Inaca Inak Ileak Inet Ica Icleak Icl Ishaker Ishab Isti]= wt_cc_model_pump(I,param,samprate)

% +++++++++++ Parameters as used by this function are:
% param(1)= V0 (resting potentail)
% param(2)= Cm (membrane capacitance)
% param(3)= Vl (Chloride leak rev. potential)
% param(4)= gl (Chloride leak conductance)
% param(5)= VK (Potassium rev. potential)
% param(6)= gl2 (Potassium leak conductance)
% param(7)= Vlic (LIC rev. potential)
% param(8)= glic (Light leak conductance)
% param(9)= Vsyn (SYN rev. potential)
% param(10)= gsyn (Synaptic leak conductance)
% param(11)= gKsm (max Shab conductance)
% param(12)= gKAm (max Shaker conductance)
% param(13)= gnewm (max novel K+ conductance)
% param(14)= gsynm (max synaptic conductance)
% Hodgkin-Huxley equations solved using ODE23S to have equally
% spaced points for analysis.

% -------------fix the time intervals by the stimulus -------------
time_length=length(I)/samprate*1000;            % duration of complete experiment in ms
time_resol=1000/samprate;                       % distance between points in a ms
time_stimulus_resol=1000/samprate;              % time resolution 0.1 ms
tspan = 0:time_stimulus_resol:time_length-time_stimulus_resol;
% tspan, time points for solutions
%
% ---------------initial values for the voltage and the conductances
V0=param(1); % typically -66 mV, resting potential
h0=1/(1+exp((-25.7-V0)/-6.4));     % start point at Shab inact. curve
n0=(1/(1+exp((-1-V0)/9.1)))^(1/2); % start point at Shab act. curve
mKA0=(1/(1+exp((-23.7-V0)/12.8)))^(1/3); % start point at Shaker act. curve
hKA0=(0.8/(1+exp((-55.3-V0)/-3.9))+0.2/(1+exp((-74.8-V0)/-10.7))); % start point at Shaker inact. curve
gnew0=1/(1+exp((-14-V0)/10.6)); % start point at new K+ conduct. act. curve
nai0 = 8; % Nai at rest is 8mM
kai0 = 140;
cai0 = 0.00016;
gKsm=param(11);     % S/cm^2, Shab max conductance
gKAm=param(12);     % S/cm^2, Shaker max conductance
gnewm=param(13);    % S/cm^2, new-potassium max conductance
Sb = 1.57*10^-5; % cell body membrane area;
Sm = 1.57*10^(-5)+ 30000*3.77*10^(-9); % the membrane area of the cell
if param(12)==0, gKAmf=0; % S/,resol,channcm^2, partial failure in Shaker inactivation is 0 if no Shaker conducntance
else gKAmf=0.087*10^-3;   % otherwise  gKAmf=0.087*10^-3;
end % S/cm^2, partial failure in Shaker inactivation
%
%----- initial setting for ODE; solving stiff differential equations, low order method.
options=odeset('RelTol',1e-4); % Error tolerance of calculations
% --------------------- integrates and stores V- and g-values ------------
for n=1:size(I,2),  % loops through the number of stimulus traces
    Im = I(:,n);  % selects the stimulus trace
    yinit=[V0;h0;n0;mKA0;hKA0;gnew0;nai0;kai0;cai0]; % initial values for diff. equations
    % [t,y]=ode23s(@cur_clamp,tspan,yinit,options,Im,time_resol,param);% ODE solver
    [t,y]=ode45(@cur_clamp,tspan,yinit,options,Im,time_resol,param);%
    
    [Inaca Ica] = NaCaPump_Body(y(:,1),y(:,9),y(:,7));
    Inak = NaKPump(y(:,1),y(:,7)); % in unit of nA
    % ------------------ Cl- channel dynamics----------------------------------
    gclm = 0.585e-4*5;
    v05 = - 0.0986;
    A = 0.0045; % RT/z'/F calculated according to G.Ugarte 2005. z' = 5.6
    gcl = gclm/(1+exp((y(1)*0.001-v05)/A));
    Icleak = param(4).*(param(3)-y(:,1)).*Sb.*1e6;
    Icl = gcl.*(param(3)-y(:,1)).*Sb.*1e6;
    
    gKs=(y(:,3).^2).*y(:,2).*gKsm; % Shab (or delayed rectifier)
    gKA=((y(:,4).^3).*y(:,5).*gKAm)+((y(:,4).^3).*gKAmf); % Shaker
    gnew=y(:,6).*gnewm; % Novel K+ conductance
    Ileak = (param(5)-y(:,1)).*param(6).*Sb.*1e6; % in unit of nA
    Ishaker = gKA.*(param(5)-y(:,1))*Sb*1e6;
    Ishab = gKs.*(param(5)-y(:,1))*Sb*1e6;
    Inew = gnew.*(param(5)-y(:,1))*Sb*1e6;
    Inet = Im + Ileak + Inaca + Inak + Ica +Icleak +Icl;
    Isti = Inet + Ishaker + Ishab + Inew;
end

% --------------------------------------------------------------------------


%%%%%%%%%%%%%%%%%%%%%%%%- WT photoreceptor model -%%%%%%%%%%%%%%%%%%%%%%%%%
function dy=cur_clamp(t,y,Im,time_resol,param)
% Current injection protocol
time_point=floor(t/time_resol+1);   % Current stimulus at time t
if time_point>length(Im), I=0; else I=Im(time_point); end
% -------------- Fixed parameters needed for the calculatios -------------
Cm=param(2);        % 4(*10^-6 F/cm^2), membrane capacitance
Vl=param(3);        % -30 mV, Chloride leak rev. potential
gl=param(4);        % 0.0585*10^-3 S/cm^2, Chloride leak conductance
VK=param(5);        % -85 mV, Potassium rev. potential
gl2=param(6);       % 0.0855*10^-3 S/cm^2, Potassium leak conductance
Vlic=param(7);      % 10 mV, LIC rev. potential
glic=param(8);      % 0 S/cm^2, Light leak conductance
Vsyn=param(9);      % -5 mV, synaptic rev. potential
gsyn=param(10);     % 0 S/cm^2, synaptic leak conductance
gKsm=param(11);     % 3*10^-3 S/cm^2, Shab max conductance
gKAm=param(12);     % 0.8*10^-3 S/cm^2, Shaker max conductance
gnewm=param(13);    % 0.11*10^-3 S/cm^2, Novel K+ channel max conductance
if param(12)==0, gKAmf=0; % S/cm^2, partial failure in Shaker inactivation
    % is 0 if no Shaker conducntance
else gKAmf=0.087*10^-3;   % otherwise  gKAmf=0.087*10^-3;
end % S/cm^2, partial failure in Shaker inactivation

% calculation of other currents

% ----------------- Activation-inactivation parameters --------------------
% Steady-state activation and inactivation curves for voltage-gated channs
% are fitted with Boltzman functions: [g/gmax]=1/(1+exp((V50-Vm)/s), where
% V50 is the voltage producing a steady-state conductance of 50% max and
% s is the slope factor of the function. Fits to experimental data (Hardie
% 1991; Hevers & Hardie, 1995).
%------------- HH activation and inactivation rate constants --------------
% the Hodgkin-Huxley rate constants are calculates as:
% alpha = {[1/(1+exp((V50-Vm)/s)]^1/n_gate}/tau/temp_Q10
% beta = {1-[1/(1+exp((V50-Vm)/s)]^1/n_gate}/tau/temp_Q10
% where n_gate is the number of gating particles, tau is the time contant
% for activation/inactivation and temp_Q10 is the temperature correction of
% 1.35 (here from 20C to 25C).
%------------- Estimation of the time constants, tau ----------------------
% Time constants were derived from experimental data. Experimental data was
% fitted with the bell shaped function:
% tau = 1/{[p1*exp((p2-Vm)/p3)] + [(p4*(p5-Vm))/(exp((p5-Vm)/p6)-1)]}
% where pi are the free parameters for fitting
%--------------------------------------------------------------------------

% ---------------------- Shab K+ channel dynamics -------------------------
% ------- Inactivation: V50=-25.7, s=-6.4; tau=1200 ms
ah=(1/(1+exp((-25.7-(y(1)))/-6.4)))/(1200/1.35);
bh=(1-(1/(1+exp((-25.7-(y(1)))/-6.4))))/(1200/1.35);
% ------- Activation: V50=-1; s=9.1;
an=((1/(1+exp((-1-(y(1)))/9.1)))^(1/2))/(1/(1.35*(0.116258*exp((-y(1)...
    -25.6551)/32.1933)+0.00659219*(-y(1)-23.8032)/(exp((-y(1)-23.8032)/...
    1.34548)-1))));
bn=((1-(1/(1+exp((-1-(y(1)))/9.1)))^(1/2)))/(1/(1.35*(0.116258*exp((-y(1)...
    -25.6551)/32.1933)+0.00659219*(-y(1)-23.8032)/(exp((-y(1)-23.8032)/...
    1.34548)-1))));

% --------------------- Shaker K+ channel dynamics ------------------------
% ------- Activation: V50=-23.7, s=12.8
amKA=((1/(1+exp((-23.7-y(1))/12.8)))^(1/3))/(1/(1.35*...
    (0.008174*exp((-y(1)+1.61882)/24.6538)+0.058139*...
    (-y(1)-59.639)/(exp((-y(1)-59.639)/4.50122)-1))));
bmKA=(1-((1/(1+exp((-23.7-y(1))/12.8))))^(1/3))/(1/(1.35*...
    (0.008174*exp((-y(1)+1.61882)/24.6538)+0.058139*(-y(1)-59.639)/...
    (exp((-y(1)-59.639)/4.50122)-1))));
% ------- Inactivation: V50(1)=-55.3, s(1)=-3.9; V50(2)=-74.8, s(2)=-10.7
% ------- two components: V50(1) contibutes 80% and V50(2) 20%
ahKA=(0.8/(1+exp((-55.3-y(1))/-3.9))+0.2/(1+exp((-74.8-y(1))/...
    -10.7)))/(1/(1.35*(0.230299*exp((-y(1)-192.973)/31.31961)+...
    0.0437316*(-y(1)+13.4859)/(exp((-y(1)+13.4859)/11.11)-1))));
bhKA=(1-(0.8/(1+exp((-55.3-y(1))/-3.9))+0.2/(1+exp((-74.8-y(1))/...
    -10.7))))/(1/(1.35*(0.230299*exp((-y(1)-192.973)/31.31961)+...
    0.0437316*(-y(1)+13.4859)/(exp((-y(1)+13.4859)/11.11)-1))));

% -------------------- Novel K+ channel dynamics --------------------------
if gnewm==0,anew=0; bnew=0; % so that gnew, y(6), would not be calculated
else    % unnecessarily, otherwise...
    anew=(1/(1+exp((-14-y(1))/10.6)))/((13+(6232/(30*sqrt(pi/2)))*...
        exp(-2*((y(1)+19.4)/30).^2))/1.35);
    bnew=(1-(1/(1+exp((-14-y(1))/10.6))))/((13+(6232/(30*sqrt(pi/2)))*...
        exp(-2*((y(1)+19.4)/30).^2))/1.35);
end

Sm = 1.57*10^(-5)+ 30000*3.77*10^(-9);
Sb = 1.57*10^(-5);% cell body membrane area , cm^2;
% ------------------ Cl- channel dynamics----------------------------------
gclm = 5*0.585e-4;
v05 = - 0.0986;
A = 0.0045; % RT/z'/F calculated according to G.Ugarte 2005. z' = 5.6
gcl = gclm/(1+exp((y(1)*0.001-v05)/A));

% ----------- Conductances
gKs=(y(3)^2)*y(2)*gKsm; % Shab (or delayed rectifier)
gKA=((y(4)^3)*y(5)*gKAm)+((y(4)^3)*gKAmf); % Shaker
gnew=y(6)*gnewm; % Novel K+ conductance

[Inaca Ica]= NaCaPump_Body(y(1),y(9),y(7)); % in unit of nA
Inak = NaKPump(y(1),y(7));

Inew = gnew*(y(1)-param(5))*Sb*1e6;
Ishaker = gKA*(y(1)-param(5))*Sb*1e6;
Ishab = gKs*(y(1)-param(5))*Sb*1e6;
Ileak = gl2*(y(1)-param(5))*Sb*1e6;

% ----------- Differential equations
% y(1) is the voltage, y(2) h, y(3) n, y(4) mKA, y(5) hKA, y(6) gnew

dy=[-y(1)*((gKA+gKs+gl+glic+gl2+gnew+gsyn+gcl)/(0.001*Cm))+(((gnew*VK)...
    +(gKA*VK)+(gKs*VK)+(gl*Vl)+(gcl*Vl)+(glic*Vlic)+(gsyn*Vsyn)+(gl2*VK))/(0.001*Cm))+((I+Inaca+Inak+Ica)/(1000*Cm*Sb))
    -y(2)*(ah+bh)+ah
    -y(3)*(an+bn)+an
    -y(4)*(amKA+bmKA)+amKA
    -y(5)*(ahKA+bhKA)+ahKA
    -y(6)*(anew+bnew)+anew
    (I*0.2054 + 3*Inaca + 3*Inak)*1000/2.92/96485 % sodium concentration
    (I*0.2401 - 2*Inak -Inew -Ishaker -Ishab -Ileak)*1000/2.92/96485 % potassium concentration
    (I*0.41 - 2*Inaca + Ica)*1000/2.92/96485/2]; % calcium concentration cell body volume = 2.826*10^-3nl, cell volume 2.92*10^-3 nl ,y(9) cai, y(8) ki y(7)nai
