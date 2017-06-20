% This program is to produce how many photons have been absorbed by different microvilli at each time 
% instant for continuous light stimulation. 

% The process of light absorption is modelled as a multinomial distribution

% Inputs:
% N_micro: the number of microvilli
% u_photon: the light time series, calibrated to be number of photons
% SaveDataFile: the file name to save the 30,000 photon series

% Outputs:
% uu: the output light pattern, summation of all photons absorbed
% Ph1: the photons absorbed by the first microvillus
% Ph2: the photons absorbed by the second microvilus
% Ph3: the photons absorbed by the third microvillus

% Zhuoyi Song, updated in 06/2017

function [uu Ph1 Ph2 Ph3]= LightAbsorptionMultinomail(N_micro,u_photon,SaveDataFile)

save(SaveDataFile,'u_photon');
tend = length(u_photon); 

% Because microvillus population is huge, if the absorbed photons are to be saved in a single matrix
% the matrix would be of a large dimention: tend*N_micro, 'outof memory'
% problem can occur. To solve this problem, this proram choped the time series to 100ms chunk
% alternatively, the matrix can be saved as a sparse matrix

re = rem(tend,100); %reminder of tend to 100,
n = fix(tend./100);

% if tend is not integer times of 100, make the
%reminder to be []; 
if re~=0
    u_photon(n*100+1:end)=[];
end

   
for j = 1:n
    
    eval(['L_ph' num2str(j) '= zeros(100,N_micro);']);
    
    for tt = (j-1)*100+1:j*100
        N_photon = u_photon(tt);
        
        lam_m = N_photon/N_micro; %average number of absorbed photons by each microvillus
        
        Pro = 1/N_micro*ones(N_micro,1);
        Absorb = mnrnd(N_photon,Pro,1);
        Absorb = Absorb';
        ti_in = rem(tt,100)+100*fix(tt./j/100);
        eval(['L_ph' num2str(j) '(ti_in,:)=Absorb;']);
    end
    
    
    save(SaveDataFile,'L_ph*','-append');
    clear L_ph*
end

    
load(SaveDataFile,'L_ph1') 

uu = []; aint = 0;Ph1 = [];Ph2=[];Ph3=[];
for  jj = 1:n
    load(SaveDataFile,['L_ph' num2str(jj)]);
    eval(['aint = L_ph' num2str(jj) ';']);
    uu = [uu;sum(aint,2)];
    Ph1 = [Ph1;aint(:,1)];
    Ph2 = [Ph2;aint(:,2)];
    Ph3 = [Ph3;aint(:,3)];
    clear aint
    clear L_ph*
end

save(SaveDataFile,'uu','Ph1','Ph2','Ph3','-append');


