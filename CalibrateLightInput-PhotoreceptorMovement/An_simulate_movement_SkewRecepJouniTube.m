% function to simulate photoreceptor light input with 2 moving dots 

% movement of optical axis has two phases:
% phase 1: move b degree in 50ms in the direction of c degree
% c = 0 means opposite of the moving dots and c = 180 means same direction
% with the dots
% phase 2: slowly return backs to original position in 375ms. Phase 2 could
% be interupted by the second dot getting close to centre of receptive
% field

% In "An_simulate_movement_SkewRecepJouniTube", all receptive field parameters for before/after movement 
% were fitted to Jouni's tube experiment (Jouni lens tube.opj).
% The corresponding parameter changes provide the scaling factor to
% the physical movement of the rhabdomere (MoveA). 

% Zhuoyi Song, updated in 06/2017
%Reference: Juusola et al 2017 ELife
function [y,ys,yskew,yskew_inter,C_x1,C_y1,t] = An_simulate_movement_SkewRecepJouniTube(hw,speed,a,b,c,di)

fs = 1000;    % sampling rate
A = 8e5/fs;   % 8e5 photons per sec, number of photons per sample
t = 0:1/fs:3;
t2p = 0.1;    % receiptive movement time to peak time 100 ms;
td  = 0.008;  % delay time before the receptive field start to move, 8ms
tph2 = 0.5; % time period for phase 2.

% motion of the two dots
D1_x = - 20 + speed*t; % light dot moving along x direction
%D1_x =  20 - speed*t; % light dot moving along x direction
D1_y = 0*t;
D2_y = 0*t;
D2_x = D1_x + di;

% centre of receptive field start at (x,y) = (0,0)
C_x = 0*t;
C_y = 0*t;


di1 = sqrt((D1_x-C_x).^2+(D1_y-C_y).^2); % distance between d1 and centre of receptive field

id = find(di1<=a*hw,1,'first'); % the time dot 1 enters "trigger zone"

if id>1 % if dot 1 enters "trigger zone"
    id1 = id + floor(td*fs); % 8ms delay until movement starts
    fprintf('Dot 1 triggered photoreceptor movement at t = %s \n',num2str(t(id)));
    
    % computing movement of receptive field
    C_x1(1:id1) = C_x(1:id1);
    C_y1(1:id1) = C_y(1:id1);
    
    % phase 1 is undisturbable
    id2 = id1+floor(t2p*fs); % end timeing point of phase 1. move in the next 100ms linearly for b along direction c
    C_x1(id1+1 : id2) = b/t2p*cos(c*pi/180)*(1/fs:1/fs:t2p);
    C_y1(id1+1 : id2) = b/t2p*sin(c*pi/180)*(1/fs:1/fs:t2p);
    
    % Phase 2 can be disturbed by dot 2, move linearly back for b along
    % direction of c, returning time perid is in 375 ms
    id22 = id2 + floor(tph2*fs);  % end timing point of phase 2 if dot 2 do not interrupt
    C_x2(id2+1 : id22) = C_x1(id2) - b/tph2*cos(c*pi/180)*(1/fs:1/fs:tph2);
    C_y2(id2+1 : id22) = C_y1(id2) - b/tph2*sin(c*pi/180)*(1/fs:1/fs:tph2);
    
    % distance between dot2 and centre of receptive field during phase1
    di2ph1 = sqrt((D2_x(id1+1:id2)-C_x1(id1+1:id2)).^2+(D2_y(id1+1:id2)-C_y1(id1+1:id2)).^2);
    
    % distance between dot2 and centre of receptive field during phase2
    di2 = sqrt((D2_x(id2+1:id22)-C_x2(id2+1:id22)).^2+(D2_y(id2+1:id22)-C_y2(id2+1:id22)).^2);
    
    id3ph1 = find(di2ph1<=a*hw,1,'first');
    
    id3 = find(di2<=a*hw,1,'first'); % dot 2 enters "trigger zone", here the assumption is id3>1,
    % which means that phase 1 movement is fast enough that dot2 only
    % enters trigger zone after phase 1. The assumption is that nothing would happen if dot2 enter
    % trigger zone during phase 1.
    
    if id3>1 %if dot 2 triggers movement
        
        id4 = id3 + floor(td*fs) + id2; % 8ms delay until movement starts
        fprintf('Dot 2 triggered photoreceptor movement at t = %s \n',num2str(t(id4)));
        C_x1(id2+1:id4) = C_x2(id2+1:id4);
        C_y1(id2+1:id4) = C_y2(id2+1:id4);
        
        id5 = id4+floor(t2p*fs); % second phase 1 movement
        % Phase1
        C_x1(id4+1:id5) = C_x2(id4) + b/t2p*cos(c*pi/180)*(1/fs:1/fs:t2p);
        C_y1(id4+1:id5) = C_y2(id4) + b/t2p*sin(c*pi/180)*(1/fs:1/fs:t2p);
        
        % Phase 2
        tb = C_x1(id5)/(b/tph2*cos(c*pi/180)); % time needed to return to original position
        tb_id = min(round(tb*fs),(length(t)-id5));
        endpoint = min(id5+round(tb*fs), length(t));
        
        C_x1(id5+1:endpoint) = C_x1(id5) - b/tph2*cos(c*pi/180)*(1/fs:1/fs:tb_id/fs);
        C_y1(id5+1:endpoint) = C_y1(id5) - b/tph2*sin(c*pi/180)*(1/fs:1/fs:tb_id/fs);
        
        C_x1(endpoint:length(t)) = 0;
        C_y1(endpoint:length(t)) = 0;
        
        length(t);size(C_x1); size(D1_x);
        
    else % dot 2 does not trigger movement
        
        C_x1(id2+1 : id22) = C_x2(id2+1 : id22);
        C_y1(id2+1 : id22) = C_y2(id2+1 : id22);
        C_x1(id22+1:length(t)) = 0;
        C_y1(id22+1:length(t)) = 0;
    end
else % if dot1 does not trigger movement
    C_x1 = C_x;
    C_y1 = C_y;
end

% % -------------When receptive field is assumed to be Gaussian--------------
% % computing actual distance between dots and centre of receptive field
% di1_r = sqrt((D1_x-C_x1).^2+(D1_y-C_y1).^2); % distance between dots and receptive field
% di2_r = sqrt((D2_x-C_x1).^2+(D2_y-C_y1).^2);
%
% di1_s = sqrt((D1_x).^2+(D1_y).^2); % distance of the dots
% di2_s = sqrt((D2_x).^2+(D2_y).^2);
%
% % light input
% % the relationship between half-width of a gaussian to it variance is
% % 2*sqrt(2*log(2))*sigma = halfwidth

% hww = hw/(2*sqrt(2*log(2)));  % variance for gausian distribution of the receptive field
% y1 = A*exp(-di1_r.^2/(2*hww^2))./hww;  % light delivered by each point object
% y2 = A*exp(-di2_r.^2/(2*hww^2))./hww;  % 800 photon per mili second
% y = y1 + y2;   % total light input.
%
%
% y1s = A*exp(-di1_s.^2/(2*hww^2))./hww;
% y2s = A*exp(-di2_s.^2/(2*hww^2))./hww;
% ys  = y1s + y2s;
%
% figure(1);hold all;
% plot(y1,'b')
% plot(y,'g');
% plot(ys,'r')


% % -------------When receptive field is assumed to be Skewed Gaussian--------------
MoveA =  max(abs(C_x1)); % maximum photoreceptor movement

% % parameter in mm
% Pgaussian = [0,17.8,51.44];
% Pgau      = [7.9419, -12.8111, 11.1723, 28.3057];

% parameter in degree
% Pgaussian = [0,23.2755,58.831];
% Pgau      = [7.9419, -16.7599, 14.616, 32.3755];

% %%parameter in degree, fitted with Jouni's data/5
% Pgaussian = [0,4.6551,100.5554];
% Pgau      = [7.9419, -3.352, 2.92, 55.3370];

% % parameter in degree, fitted with Jouni's data/5.75
% Pgaussian = [0,4.0479,107.8336];
% Pgau      = [7.9419, -2.9148, 2.5419, 59.3424];

% %% parameter for 1/6.8 of jouni's experimental data: x    = SideExpFit(:,1)*16/34*2.78/6.8
% % parameter in degree, fitted with Jouni's data/6.8
% Pgaussian = [0,3.4229,117.2667];
% Pgau      = [7.9419, -2.4647, 2.1494, 64.5335];

%% parameter for 1/8 of jouni's experimental data: x    = SideExpFit(:,1)*16/34*2.78/6.8
% parameter in degree, fitted with Jouni's data/8
Pgaussian = [0,2.9094,127.1936];
Pgau      = [7.9419, -2.095, 1.827, 69.996];

% Diameter of experimental lens for receptive field movement is 34mm
% Diameter of ommatidia is 16 um
% The degree movement is calculated as x*16/34*2.78 degree, every um
% corresponds to 2.78 degree of movement

DegreeMove = 2*16/34*2.78; % 2.6165 degree, 1um corresponds to 2.6125degree

AlphaChange = (Pgau(1)-0)/DegreeMove; % parameter(degree) change this amount every degree of center movement
MeanChange  = (Pgau(2) - Pgaussian(1))/DegreeMove;
GainSigma = (Pgau(3))/DegreeMove/Pgaussian(2);
% caldulate the relationship between parameters of skewed gaussian to
% DegreeMove

Alpha = MoveA*AlphaChange;
% paramMean labels the off-axis amount of between the center of receptive
% field and the center of photoreceptor optical axis
paramMean = MoveA*MeanChange; 

size(D1_x)
size(C_x1)
di1_r = sqrt((D1_x-C_x1).^2+(D1_y-C_y1).^2);

di2_r = sqrt((D2_x-C_x1).^2+(D2_y-C_y1).^2);
di1_s = sqrt((D1_x).^2+(D1_y).^2); % distance of the dots
di2_s = sqrt((D2_x).^2+(D2_y).^2);

% HalfWidth
% the relationship between half-width of a gaussian to it variance is
% 2*sqrt(2*log(2))*sigma = halfwidth

hww  = hw/(2*sqrt(2*log(2)));  % variance for gausian distribution of the receptive field
hwwM = hww.*(1-0.5*sigmf(MoveA,[100,0.6])); % for bigger the shift, smaller the receptive field.
hwwMH = hwwM*2*sqrt(2*log(2))  % half-width

pulse = ones(4,1)/4;

y1_Gauss = (1/sqrt((2*pi)./hwwM).*exp(-(di1_r-paramMean).^2./2./hwwM.^2));
y1_SkewedGauss = 2*A*y1_Gauss.*normcdf(Alpha*((di1_r-paramMean)./hwwM))/hwwM;


y2_Gauss = (1/sqrt((2*pi)./hwwM).*exp(-(di2_r-paramMean).^2./2./hwwM^2));
y2_SkewedGauss = 2*A*y2_Gauss.*normcdf(Alpha*((di2_r-paramMean)./hwwM))/hwwM;

yskew = y1_SkewedGauss + y2_SkewedGauss;
yskew =  conv(pulse,yskew);
%yskew = yskew + 0.03*max(yskew);

y1 = A*exp(-di1_r.^2/(2*hww^2))*(1/sqrt((2*pi)./hww));  % light delivered by each point object
y2 = A*exp(-di2_r.^2/(2*hww^2))*(1/sqrt((2*pi)./hww));  % 800 photon per mili second
y = y1 + y2;   % total light input.
y =  conv(pulse,y);
%y = y + 0.03*max(y);

%figure(2);hold all; plot(di1_r);plot(di2_r)

% figure(3);hold all;
% plot(y1_Gauss,'b');
% plot(y1_SkewedGauss/A/2,'r')


y1s = A*exp(-di1_s.^2/(2*hww^2)).*(1/sqrt((2*pi)./hww));
y2s = A*exp(-di2_s.^2/(2*hww^2)).*(1/sqrt((2*pi)./hww));
ys = y1s + y2s;
ys =  conv(pulse,ys);
%ys = ys + 0.03*max(ys);

% figure(3);hold all;
% plot(ys,'r');
% plot(y,'k')
% plot(yskew);
% plot(y1_SkewedGauss,'b');plot(y2_SkewedGauss,'k')


%% Interplated receptive field function. 
%% The previous receptive field is a Gaussian
%% fit to the experimental recording. But actual receptive field may not be exactly gaussian, reflected by
%% the longer tail in the receptive field. It was in the hope at one time that
%% by interplolating the real measured data, the biophysical model response will approxiate better
%% real cell recordings, especially the early rising phase of the data. But actually, it does not matter. 

load r2real
AnRecepField = r2real;

%load AnRecepField

xdata = AnRecepField(:,1);
ydata = AnRecepField(:,2) - min(AnRecepField(:,2));
ydata = 2*3*A*ydata/sum(ydata);

y1s_inter = fixpt_interp1(xdata,ydata,di1_s,float('single'),2^-3,float('single'),...
    2^-14);
y2s_inter = fixpt_interp1(xdata,ydata,di2_s,float('single'),2^-3,float('single'),...
    2^-14);
ys_inter = (y1s_inter + y2s_inter);
ys_inter = ys_inter - ys_inter(1) ;
%ys_inter = ys_inter*max(ys)/max(ys_inter);

y1_inter = fixpt_interp1(xdata,ydata,di1_r,float('single'),2^-3,float('single'),...
    2^-14);
y2_inter = fixpt_interp1(xdata,ydata,di2_r,float('single'),2^-3,float('single'),...
    2^-14);
y_inter = y1_inter + y2_inter;
y_inter = y_inter - y_inter(1) ;
%y_inter = y_inter*max(ys)/max(y_inter);

ySkewdata = ydata*2*hwwM/hww;
y1skew_inter = fixpt_interp1(xdata,ySkewdata,(di1_r- paramMean)*hww/hwwM,float('single'),2^-3,float('single'),...
    2^-14);
y2skew_inter = fixpt_interp1(xdata,ySkewdata,(di2_r- paramMean)*hww/hwwM,float('single'),2^-3,float('single'),...
    2^-14);
yskew_inter = (y1skew_inter + y2skew_inter);
yskew_inter = conv(pulse,yskew_inter);
yskew_inter  = yskew_inter - yskew_inter(1) ;

figure(1);hold all;
%plot(ys,'b')
plot(yskew_inter,'r');

plot(yskew,'g')
end






