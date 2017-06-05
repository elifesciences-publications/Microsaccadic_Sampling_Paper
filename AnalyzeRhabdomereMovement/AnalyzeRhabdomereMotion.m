%Calculation of avarage movement in position of rhabdomeres
%Input is a tiff stack, where every frame is an individual tiff.
%Jouni Takalo 2017 
%Article ref

dirname = 'C:\Users\Mikko\Desktop\musca\may16\1\tiff'; %%folder of tiff stack
filename ='pu10_0_780'; %% base file name of stack
filename =[filename '*'];
%Loading tiff stack
HH =dir(fullfile(dirname,filename)); 
Toff =length(HH);
baseFileName = HH(1).name;
tiffinfo =imfinfo(fullfile(dirname, baseFileName));
LL =tiffinfo.Height;
WW =tiffinfo.Width;
video =zeros(LL,WW,Toff);
for k = 1:Toff
    baseFileName = HH(k).name;
     fullFilename = fullfile(dirname, baseFileName);
      video(:,:,k) = imread(fullFilename, 'Info', tiffinfo);
     video(:,:,k) =video(:,:,k)-median(median(video(:,:,k)));% removing median of the image
end
video(video<0) =0; %%removing median backround noise based

refframe =149; % reference image frame
row =zeros(1,Toff); % Initiliaze movement vector in up/down direction
col =zeros(1,Toff); % Initiliaze movement vector in left/right direction
%Loop thought the tiff stack
for i =1:Toff
    %Calculate 2D correlation map between currren frame and reference frame
    A = xcorr2(video(y,x,refframe)-mean(mean(mean(video(y,x,refframe)))),video(y,x,i)-mean(mean(mean(video(y,x,i))))); 
    [rowi,coli] = find(A>=max(max(A)*0.95)); %Take points of 2D correlation map, which are 95% or higher compared to max of the correlation map
    %Take values of A, which have 95% or higher
    B =zeros(length(rowi),1);
    for j =1:length(rowi)
        B(j) =A(rowi(j),coli(j));
    end
    %Calculate mean position of peak in 2D correlation map and using  2D
    %correlation map as weights
    col(i) = sum(coli.*B)/sum(B);
    row(i) = sum(rowi.*B)/sum(B);
end
%Delete the point of reference frame
row = row -row(refframe);
col = col -col(refframe);
%row gives movement vector in up/down direction
%col gives movement vector in left/right direction