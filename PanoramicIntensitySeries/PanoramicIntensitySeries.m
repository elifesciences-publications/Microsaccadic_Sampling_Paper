%Script for calculating intensity series from panoramic images based on yaw
%of fly
%Jouni Takalo 2017 
%Article ref
%This script takes two input:
%rawimage: gray scale panoramic image
%yaw: the yaw value of fly
%Output:
%data: 15 intensity series
totsize =size(rawimage);
quatersize =round(totsize(1)/4);
yawwrap = wrapTo360(yaw);
pickpoint = round(yawwrap/360*totsize(2))+1;

k =1;
for i =quatersize: round(quatersize*2/15):quatersize*3;
lineim=squeeze(rawimage(i,:));
data(:,k) =(double(lineim(pickpoint))/255).^2.5;
k=k+1;
end
