% function to simulate photoreceptor light input with 2 moving dots as
% stimulus and *** optical axis movement ***
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

%Reference: Juusola et al 2017 ELife

function [y,C_x1,C_y1,t] = An_simulate_r16_movement2(hw,speed,a,b,c)

di = 8.525; % di: distance between two dots: (8.525 degree as standard)
fs = 1000; % sampling rate
A = 8e5/fs; % 8e5 photons per sec
t = 0:1/fs:1;
riseb = 0.05; % phase 1: move riseb degree in 100 ms
risec = 0.375; % phase 2: move risec degree in 800 ms

% motion of the two dots
D1_x = 30 - speed*t;
D1_y = 0*t;
D2_y = 0*t;
D2_x = D1_x + di;

% centre of receptive field start at (x,y) = (0,0)
C_x = 0*t;
C_y = 0*t;


di1 = sqrt((D1_x-C_x).^2+(D1_y-C_y).^2); % distance between d1 and centre of receptive field
id = find(di1<=a*hw,1,'first'); % dot 1 enters "trigger zone"
if id>1 % if dot 1 enters "trigger zone"
    id1 = id + floor(0.01*fs); % 10ms delay until movement starts
    fprintf('Dot 1 triggered photoreceptor movement at t = %s \n',num2str(t(id)));
    % computing movement of receptive field
    C_x1(1:id1) = C_x(1:id1);
    C_y1(1:id1) = C_y(1:id1);
    % phase 1 is undisturbable
    id2 = id1+floor(riseb*fs);
    C_x1(id1+1 : id2) = b/riseb*cos(c*pi/180)*(1/fs:1/fs:riseb);
    C_y1(id1+1 : id2) = b/riseb*sin(c*pi/180)*(1/fs:1/fs:riseb);
    
    % Phase 2 can be disturbed by dot 2
    id22 = id2 + floor(risec*fs);
    C_x2(id2+1 : id22) = C_x1(id2) - b/risec*cos(c*pi/180)*(1/fs:1/fs:risec);
    C_y2(id2+1 : id22) = C_y1(id2) - b/risec*sin(c*pi/180)*(1/fs:1/fs:risec);
    
    % distance between d1 and centre of receptive field during phase2
    di2 = sqrt((D2_x(id2+1:id22)-C_x2(id2+1:id22)).^2+(D2_y(id2+1:id22)-C_y2(id2+1:id22)).^2);
    
    id3 = find(di2<=a*hw,1,'first'); % dot 2 enters "trigger zone"
    
    if id3>1 %if dot 2 triggers movement
        id4 = id3 + floor(0.01*fs) + id2; % 10ms delay until movement starts
        fprintf('Dot 2 triggered photoreceptor movement at t = %s \n',num2str(t(id3)));
        C_x1(id2+1:id4) = C_x2(id2+1:id4);
        C_y1(id2+1:id4) = C_y2(id2+1:id4);
        id5 = id4+floor(0.05*fs);
        % Phase1
        C_x1(id4+1:id5) = C_x2(id4) + b/riseb*cos(c*pi/180)*(1/fs:1/fs:riseb);
        C_y1(id4+1:id5) = C_y2(id4) + b/riseb*sin(c*pi/180)*(1/fs:1/fs:riseb);
        
        % Phase 2
        tb = C_x1(id5)/(b/risec*cos(c*pi/180)); % time needed to return to original position 
        tb_id = floor(tb*fs);
        
        C_x1(id5+1:id5+tb_id) = C_x1(id5) - b/risec*cos(c*pi/180)*(1/fs:1/fs:tb);
        C_y1(id5+1:id5+tb_id) = C_y1(id5) - b/risec*sin(c*pi/180)*(1/fs:1/fs:tb);
        
        C_x1(id5+tb_id:length(t)) = 0;
        C_y1(id5+tb_id:length(t)) = 0;
    
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


% computing actual distance between dots and centre of receptive field
di1_r = sqrt((D1_x-C_x1).^2+(D1_y-C_y1).^2);
di2_r = sqrt((D2_x-C_x1).^2+(D2_y-C_y1).^2);

% light input
hww = hw/(2*(2*log(2))^0.5);
y1 = A*exp(-di1_r.^2/(2*hww^2));  % light delivered by each point object
y2 = A*exp(-di2_r.^2/(2*hww^2));  % 800 photon per mili second
y = y1 + y2;   % total light input.

% plotting
s = sum(y);
figure()
head = sprintf('Light distribution over time, total %s photons',num2str(floor(s)));
plot(t,y), title(head);
figure()
plot(t,C_x1),hold all, plot(t,C_y1);
title('X- and Y-locations of receptive field over time');legend('x-locs','y-locs')


end






