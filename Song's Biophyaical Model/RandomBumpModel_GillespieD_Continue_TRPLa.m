% Simulation of phototransduction cascade inside a single microvillus

% The input nstep is the simulation steps, as the time interval of
% simulation step is 0.1ms, nstep = 10*length(LL). LL is the
% photon inpurt series, the time interval between two points of LL is 1ms.
% y is the initial state of the transduction state,
% para is the parameter of the transduction model

% The output tstep stores the real time of simulation, yy is the output
% matrix, with each row store the transduction sate at each time point.

% Zhuoyi Song, updated in 06/2017

function [tstep, yy] = RandomBumpModel_GillespieD_Continue_TRPLa(nstep,y,LL,para,fnstr,la)
%%%%LL is the light input
%% Ca dynamics is fast and deterministic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% to execute run
% [t,yy]=modelF_st_2(nsteps,[B* Ca G* Rh* PLC* A* C* G_t 0 0],C_ex);
%
% ** nstep : number of Gillespie steps (nstep > 10*length(LL) is appropriate).
%
% ** [B* Ca G* Rh* PLC* A* C* G_t 0 0 ] initial condition :
%   B* = 0
%   Ca = 1
%   G* = 0
%   Rh* = 1 : one Rh activated at t = 0
%   PLC* = 0
%   A* = 0
%   C* = 0
%   G_t = 50 (number of available G-proteins).
%
% ** C_ex : external Calcium contentration.
% In our units, physiological (1.5mM) corresponds to 5000.

% WT parameters used here :
%--------------------------

% % feedback parameters :
%
kBp =para(1);% Kp
kBn = para(2);
mp=para(3);
mn=para(4);

% Rhodopsin deactivation :

gamma_Rh=para(5);
gRh =para(6);

% G-dynamics :
kap_G = para(7);
gamma_Gi = para(8);

% nb : G_t specified by the initial condition.
% for PLC* dynamics :

gamma_P = para(9);
gplc = para(10);
kap_P = para(11);
P_t = para(12);

% for A* dynamics :

nd= para(13); % cooperativity of the A* to activate B*.
kap_A = para(14);
gamma_A = para(15);
gAn = para(16);
gAp = para(17);

% B* dynamics :

kap_B=para(18); % kap_B
K_A = para(19); % K_A ** arbitrary decomposition of the ratio kap_B/K_A^3.

gBn = para(20);
gBp = para(21);
gamma_B = para(22);

% C* dynamics :

K_U = para(23);
K_R = para(24);


% Ca dynamics :

K_Ca = para(25);
T_t = para(26);
gamma_GAP = para(27);

G_t = y(8);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%DYNAMICAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% y(1) -- B*
% y(2) -- Ca or [cai]i
% y(3) -- G*
% y(4) -- Rh*
% y(5) -- PLC*
% y(6) -- A*
% y(7) -- C*
% y(8) -- G
% y(9)=fbn;
% y(10)=fbp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%% LOOP %%%%%%%%%%%%%%%%%%%%%%%%%
sam = 10; % This is the number of samples in each ms.
ii=1/sam;
tstep = zeros(100000,1);

yy = zeros(nstep/sam,size(y,2)); % to store the data on ms basis

yy(1,:) = y;
k=1;
T_rh = find(LL~=0); % The matrix when light input are not zeros
T_rhin = 1;
t_rho = zeros(length(T_rh)+1,1);%t_rho is in ms base, index of T_rh
starh = 0; % to see if Rh is updated according to photostimuli, if y(4) is incremented, starh plus 1
TRP_Rev_pre = 0.009; %Volts
while ii<nstep/sam;  %ii is the count of ms
    %%%%%%%%%%%   LOOP %%%%%%%%%%%%%%%%%
    ii;
    k=k+1;
    save_in_now = ceil(ii); % this is the index of ii, when it is counted in 1ms base
    %% Feedback functions
    fbp=(y(2)/kBp)^mp/(1+(y(2)/kBp)^mp);
    fbn=fnstr*(y(7)/kBn)^mn/(1+(y(7)/kBn)^mn);
    
    z(1,1) = kap_B*((y(6)/K_A)^nd)*(1+gBp*fbp)*(T_t-y(1));
    z(1,2)=y(1)*gamma_B*(1+gBn*fbn);
    
    z(2,1)=0;z(2,2)=0;
    %% G-protein
    z(3,1)=kap_G*y(8)*y(4);
    % removal of G* is linked to PLC* activation
    z(3,2)=gamma_GAP*y(3)*y(5);
    %removal of G is linked to G* activation
    z(8,2)=0;
    %recovery of G from a refractory pool
    z(8,1)=gamma_Gi*(G_t-y(5)-y(3)-y(8)); % a5 = gamma_Gi*(G_t-PLC*-G*-G);
    %% Rh* inactivation
    if rem(ii,1)==0 && t_rho(T_rhin)==0
        rh_in = ii;
        y4_pre = y(4);
        y(4) = LL(rh_in)+y(4);
        t_rho(T_rhin)=1;
        if y(4)>y4_pre
            starh = starh+1;% to see if Rh is updated according to photostimuli, if y(4) is incremented, starh plus 1
            Rh_Up = 1; % To indicate that Rhodopsin has just been updated
        end
    end
    
    z(4,1)=0;
    z(4,2)=y(4)*(1+gRh*fbn)*gamma_Rh;
    z(5,1)=kap_P*(P_t-y(5))*y(3);
    z(5,2)=gamma_P*y(5)*(1+gplc*fbn);
    
    %% A* :
    z(6,1)=kap_A*y(5)*(1+gAp*fbp);
    z(6,2)=gamma_A*y(6)*(1+gAn*fbn);
    
    %%  Ca buffer which drives negative feedback (C*)
    Oc = y(7)/0.5/6.022/3/100;
    z(7,1)=K_U*y(2)*(1-Oc); %the number of bined calmodulin
    z(7,2)=K_R*y(7);
    Icam = 1.16*(z(7,1)-z(7,2))/1806.6;
    %% Ca
    
    % The calculations of iteration of TRP_Rev is wrong, the result is up to
    % 1 volts, which is wrong
    Ntrp = y(1);
    Cai = y(2)/6.022/3/100; % 1/6.022/2/100 is the scaling factor to change mM to numbers in 3*10*(-9) nl volumn of a microvillus
    [I_m ICa INa IMg IK TRP_Rev] = TRP_Rev_GHK_online(Ntrp,Cai,TRP_Rev_pre);
    TRP_Rev_pre = TRP_Rev;
    Inaca = NaCaPump(-0.07,Cai);% in unit of pA
    
    y(2) = 6.022*3*100*((ICa)*1e9/1.002/2/96485/3/1000 + 2*K_R*Oc + 2e-4)/(K_Ca+2*K_U*(1-Oc)+7.24); %if Knaca = 4*10^-3, the parameter for NaCa exchanger to have pA
    % formula for steady state Calcium \frac{Ica,net/(2*V*F)+n*Bi*K_R*Oc}{K_Ca+n*Bi*K_U*(1-Oc)}
    y(9)=fbn;
    y(10)=fbp;
    y(11) = I_m;
    y(12) = ICa;
    y(13) = TRP_Rev;
    y(14) = Inaca;
    y(15) = I_m +Inaca-Icam;
    clear I_m ICa TRP_Rev Inaca Icam
    %      y(14) = INa;
    %      y(15) = IMg;
    %        y(15) = Icam;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% stochastic update
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    up=z(:,1);
    norm_1=sum(up);
    if norm_1>0
        zup=cumsum(up)/norm_1;      % Establish up-wheel
    end
    down=z(:,2);
    norm_2=sum(down);
    if norm_2>0
        zdown=cumsum(down)/norm_2;  % Establish down wheel
    end
    
    if norm_1+norm_2<0
        break;
    end
    
    if norm_1+norm_2==0 %in this case, simulation is continued without the variables updated
        yy(save_in_now,:)=y;
        ii = ii+1/sam;
        continue;
    elseif norm_1+norm_2>0
        r=rand;
        q=norm_1/(norm_1+norm_2);   % Probability: there is a up change
        if r<q                      % if up happens
            r=r/q;                  % normalizing: 100% up happened
            tst=find(zup>r); nn=tst(1);  % which reaction happens
            y(nn)=round(y(nn))+1;            %  Increase by one
            if nn==5
                y(3)=y(3)-1;        % G* down
            end
            if nn==3
                y(8)=round(y(8))-1;        % one G protein down
            end
        else                     % if down happens
            r=(r-q)/(1-q);         % normalizing: 100% down happened
            tst=find(zdown>r);nn=tst(1); % find which reaction down
            y(nn)=round(y(nn))-1;               % reaction down by one
            if nn==4 && Rh_Up==1
                y(nn)=round(y(nn))+1;
                Rh_Up = 0; % to change back Rh_Up to zero.
            end
        end
    end
    
    r2=0.001+(0.999)*rand;
    tt = ii+log(1/r2)/abs(la+(norm_1+norm_2));
    t_update = round(tt*sam)/sam; % this is in ms basis, but with 0.1ms precision
    
    %update the next simulation timestep and store information
    if T_rhin<length(T_rh) || T_rhin==length(T_rh)
        if   t_update < T_rh(T_rhin)|| t_update == T_rh(T_rhin)
            ii = t_update;  % only update if there is no photon stimuli in the next dt
        elseif  t_update > T_rh(T_rhin)
            ii = T_rh(T_rhin);
            T_rhin = T_rhin+1; % the rolling index of nonzero light input
        end
    elseif T_rhin > length(T_rh) && t_update < nstep/sam
        ii = t_update;  % only update if there is no photon stimuli in the next dt
    elseif T_rhin > length(T_rh) && t_update >= nstep/sam
        ii = nstep/sam;
    end
    tstep(k) = ii; %t_step saves the update states
    save_in_future= ceil(ii);
    yy(save_in_now:save_in_future,:) = repmat(y,save_in_future-save_in_now+1,1);
    clear z
end
