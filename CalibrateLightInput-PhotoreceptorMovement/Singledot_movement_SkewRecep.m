% function to simulate photoreceptor light input for a single moving dot
% with optical axis movement

% hw: the half-width of receptive field
% a is the "trigger zone": optical axis starts moving when one object is
% within a*hw from the receptive field centre

% movement of optical axis has two phases:
% phase 1: move b degree in 50ms in the direction of c degree
% c = 0 means opposite of the moving dots and c = 180 means same direction
% with the dots

% phase 2: slowly return backs to original position in 375ms. Phase 2 could
% be interupted by the second dot getting close to centre of receptive
% field
% Zhuoyi Song, updated in 06/2017
%Reference: Juusola et al 2017 ELife
function [y,ys,yskew,yskew_inter,C_x1,C_y1,t] = Singledot_movement_SkewRecep(hw,speed,a,b,c)

fs = 1000;    % sampling rate
A = 8e5/fs;   % 8e5 photons per sec, number of photons per sample
t = 0:1/fs:3;
t2p = 0.1;    % receiptive movement time to peak time 100 ms;
td  = 0.008;  % delay time before the receptive field start to move, 8ms
tph2 = 0.5; % time period for phase 2.

% motion of the dot
D1_x = - 20 + speed*t; % light dot moving toward object along x direction
D1_y = 0*t;

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
    
    C_x1(id2+1 : id22) = C_x2(id2+1 : id22);
    C_y1(id2+1 : id22) = C_y2(id2+1 : id22);
    C_x1(id22+1:length(t)) = 0;
    C_y1(id22+1:length(t)) = 0;
    
else % if dot1 does not trigger movement
    C_x1 = C_x;
    C_y1 = C_y;
end



%% light input when receptive field has moved, the shape of receptive feild is skewed (skewed gaussian)
%% The parameters corresponds to the maximum movemnt position

MoveA =  max(abs(C_x1));

% % parameter in mm
% Pgaussian = [0,17.8,51.44];
% Pgau      = [7.9419, -12.8111, 11.1723, 28.3057];

% parameter in degree
% Pgaussian = [0,23.2755,58.831];
% Pgau      = [7.9419, -16.7599, 14.616, 32.3755];

% parameter in degree, fitted with Jouni's data/5
%Pgaussian = [0,4.6551,100.5554];
%Pgau      = [7.9419, -3.352, 2.92, 55.3370];

% % parameter in degree, fitted with Jouni's data/5.75
% Pgaussian = [0,4.0479,107.8336];
% Pgau      = [7.9419, -2.9148, 2.5419, 59.3424];

%% parameter for 1/6.8 of jouni's experimental data: x    = SideExpFit(:,1)*16/34*2.78/6.8
% parameter in degree, fitted with Jouni's data/6.8
Pgaussian = [0,3.4229,117.2667];
Pgau      = [7.9419, -2.4647, 2.1494, 64.5335];

% Diameter of experimental lens for receptive field movement is 34mm
% Diameter of ommatidia is 16 um
% The degree movement is calculated as x*16/34*2.78 degree, every um
% corresponds to 2.78 degree of movement

DegreeMove = 2*16/34*2.78; % 2.6165 degree, 1um corresponds to 2.6125degree

AlphaChange = (Pgau(1)-0)/DegreeMove; % parameter in degree change this amount every degree of center movement
MeanChange  = (Pgau(2) - Pgaussian(1))/DegreeMove;
GainSigma = (Pgau(3))/DegreeMove/Pgaussian(2);

Alpha = MoveA*AlphaChange;
paramMean = MoveA*MeanChange;

size(D1_x)
size(C_x1)
di1_r = sqrt((D1_x-C_x1).^2+(D1_y-C_y1).^2);

di1_s = sqrt((D1_x).^2+(D1_y).^2); % distance of the dot


% HalfWidth
% the relationship between half-width of a gaussian to it variance is 2*sqrt(2*log(2))*sigma = halfwidth
hww  = hw/(2*sqrt(2*log(2)));  % variance for gausian distribution of the receptive field
hwwM = hww*MoveA*GainSigma;
hwwMH = hwwM*2*sqrt(2*log(2));

y1_Gauss = (1/sqrt((2*pi)./hwwM).*exp(-(di1_r-paramMean).^2./2./hwwM.^2));
y1_SkewedGauss = 2*A*y1_Gauss.*normcdf(Alpha*((di1_r-paramMean).^2./hwwM.^2));

yskew = y1_SkewedGauss + 0.03*max(y1_SkewedGauss);

y1 = A*exp(-di1_r.^2/(2*hww^2))*(1/sqrt((2*pi)./hww));  % light delivered by each point object
y = y1 + 0.03*max(y1);   % total light input.

y1s = A*exp(-di1_s.^2/(2*hww^2)).*(1/sqrt((2*pi)./hww));
ys = y1s + 0.03*max(y1s);

%% Interplated the receptive field function. 
%% The previous receptive field is a Gaussian fit to the experimental recording. 
%% But actual receptive field may not be exactly gaussian, reflected by
%% the longer tail in the receptive field. It was in the hope that
%% by interplolating the real measured data, the biophysical model response will approxiate better
%% real cell recordings, especially the early rising phase of the data. It is actually not the case

load AnRecepField
xdata = AnRecepField(:,1);
ydata = AnRecepField(:,2) - min(AnRecepField(:,2));
ydata = 4*A*ydata/sum(ydata);

ys_inter = fixpt_interp1(xdata,ydata,di1_s,float('single'),2^-3,float('single'),...
    2^-14);
ys_inter = ys_inter - ys_inter(1) + 20;

y_inter = fixpt_interp1(xdata,ydata,di1_r,float('single'),2^-3,float('single'),...
    2^-14);
y_inter = y_inter - y_inter(1) + 20;

ySkewdata = ydata*2*hwwM/hww;

yskew_inter = fixpt_interp1(xdata,ySkewdata,(di1_r- paramMean)*hww/hwwM,float('single'),2^-3,float('single'),...
    2^-14);
yskew_inter  = yskew_inter - yskew_inter(1) + 20;

end






