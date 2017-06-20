% This is the main program that simulates light input for the biophysical
% model, taking into account the receptive field movement effects.

% Different with 'SimulateAn', the photoreceptor movement programe uses 
% 'An_simulate_movement_SkewRecepJouniTube'. This model considers the effect of
% cone cell movement, which effectively make the center of the receptive
% field misaligned with the center of the optical axis. 

% Zhuoyi Song, updated in 06/2017

%Reference: Juusola et al 2017 ELife
kk = 1;
%Sp = [600,650,700,800,1000,1200]; % speed of moving light dots, in degree/s
Sp = [10,20,45,60,100,200,300,400,600]; % speed of moving light dots, in degree/s
b = 1.6; % movement range, in degree
di = 2;  % two dot distance, in degree

for sp = 1:length(Sp)
    %for di=1:4
    speed = Sp(sp);
    hw     = 8.0603;  % half width of receptive field
    a     = 9/5;  % tigger zone/hw = 7/5, trigger zone is the width of dot enters the field when the cell starts to respond
    c     = 0;    % alpha, direction of receptive field movement
    FixedParams.hw = hw;
    FixedParams.a  = a;
    FixedParams.c  = c;
    FixedParams.b(kk) = b;
    %FixedParams.di  = di;
    VaryParams.speed(kk) = speed;
    VaryParams.di(kk) = di;
    
    [y,ys,yskew,yskew_inter,C_x1,C_y1,t] = An_simulate_movement_SkewRecepJouniTube(hw,speed,a,b,c,di);
    %[y,ys,yskew,yskew_inter,C_x1,C_y1,t] = Singledot_movement_SkewRecep(hw,speed,a,b,c);

    %Response.IntensityMove(kk,:)     = y;  % light stimulus with gaussian receptive field movement
    Response.IntensitySkew(kk,:)      = yskew; % light stimulus with skewed receptive field movement
    %Response.IntensitySkewInter(kk,:) = yskew_inter; % interpolated light stimulus
    %Response.IntensityRefer(kk,:)     = ys; % reference light stimulus without photoreceptor movement
    Response.xpos(kk,:)               = C_x1;
    Response.ypos(kk,:)               = C_y1;
    Response.t(kk,:)                  =t;
    kk = kk+1;
    %end
end

save LightInput_di20_c0_JouniExp FixedParams VaryParams Response

