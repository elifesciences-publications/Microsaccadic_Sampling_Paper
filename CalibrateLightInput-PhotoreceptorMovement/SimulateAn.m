% This is a very simple program to investigate how the photoreceptor
% movement will impact on the shape of two dot light stimuli.

% This programe uses an earlier version of the the photoreceptor movement model (An_simulate_r16_movement),
% which only considers symmetrical gaussian receptive field, and receptive
% field center moves together with the physical photorecepotor optical
% axis.

% Zhuoyi Song, updated in 06/2017

%Reference: Juusola et al 2017 ELife
clear all
kk = 1;
%Sp = 100:20:600; % speed of dot movement, in degrees/s
Sp = [400,500,600];
for sp = 1:length(Sp)
    for b = 1.6
        speed = Sp(sp);
        kk
        hw    = 5;   % half width of receptive field
        %speed = 50;  % speed of moving light dots, in degree/s
        a     = 7/5; % tigger zone/hw = 7/5, trigger zone is the width of dot enters the field when the cell starts to respond
        %b     = 5;   % movement range, in degree
        c     = 0;   % alpha, direction of receptive field movement
        FixedParams.hw = hw;
        FixedParams.a  = a;
        FixedParams.c  = c;
        VaryParams.b(kk) = b;
        VaryParams.speed(kk) = speed;
        [y,ys, C_x1,C_y1,t] = An_simulate_r16_movement(hw,speed,a,b,c);
        Response.Intensity(kk,:)          = y;  % light intensity profile with photoreceptor movement
        Response.IntensityRefer(kk,:)     = ys; % reference light intensity profile, i.e. the light intensity profile if no photoreceptor movement
        Response.xpos(kk,:)               = C_x1;
        Response.ypos(kk,:)               = C_y1;
        Response.t(kk,:)                  =t;
        
        kk = kk+1;
    end
end

%save Intensity_c150_sim FixedParams VaryParams Response