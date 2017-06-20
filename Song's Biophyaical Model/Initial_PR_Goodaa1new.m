% This files set the initial parameters for the biophysical model 

% Zhuoyi Song, updated in 06/2017

%Reference: Juusola et al 2017 ELife
% Parameters regulated to produce the same dynamics as Postma 1999 paper
yini =[0    1       0        0        0        0         0   50       1         0];
yini_online = [yini,zeros(1,5)];
para = zeros(29,1);
% feedback parameters :
aa = 1;

% parameters for when calcium is 1-1.5mM
para(1) =0.3*6.022*3*100;% kBp =para(1);% Kp %different
para(2) = 0.18*6.022*3*100; % kBn = para(2);

para(3) = 2; % mp=para(3);
para(4) = 3; % % mn=para(4);

% % Rhodopsin deactivation :

para(5) = 0.0037*aa;% gamma_Rh=para(5);
para(6) = 40*1;% gRh =para(6);

% % G-dynamics :
para(7) = 0.0047*0.3*5;% kap_G = para(7);
para(8) = 0.0035*aa;% gamma_Gi = para(8);% may be worth to include the calcium positive feedback

% % nb : G_t specified by the initial condition.
% % for PLC* dynamics :
para(9) = 0.048*aa*3; % gamma_P = para(9);
para(10) = 11.1;% gplc = para(10);
para(11) = 0.0039*aa*4;% kap_P = para(11);
para(12) = 100;% P_t = para(12);

% % for A* dynamics :
para(13) = 2; % nd = para(13), the cooperativity of A* to B*
para(14) = 1.3*aa;% kap_A = para(14);
para(15) = 0.004*aa;% gamma_A = para(15);
para(16) = 37.8; % gAn = para(16);
para(17) = 0;% gAp = para(17);

% % B* dynamics :
para(18) = 0.15*aa; % kap_B=para(18); % kap_B
para(19) = 100; % K_A = para(19); % K_A ** arbitrary decomposition of the ratio kap_B/K_A^3.

%
para(20) = 10;% gBn = para(20);
para(21) = 11.5;% gBp = para(21);
para(22) = 0.025*aa; % gamma_B = para(22);

% % C* dynamics :

para(23) = 0.03; % K_U = para(23);
para(24) = 0.0055; % K_R = para(24);
%
% % Ca dynamics :

para(25) = 1; % K_Ca = para(25);
para(26) =27; % T_t = para(26);

para(27) = 0.003*aa; % gamma_GAP = para(27);
