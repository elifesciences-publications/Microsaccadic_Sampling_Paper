% This programme simulates the whole buntch of 30,000 microvilli and sum
% the LIC together.

% Zhuoyi Song, updated in 06/2017

tim = 1; % tim controls simulation trials

% set stochastic simulator
rndseed = RandStream.create('mt19937ar','seed',sum(10000*tim+sum(clock)));
RandStream.setGlobalStream(rndseed);

% load photon stimuli file
Photonfile = 'LightSkew_b2p3_speed205';
SS = who('-file',Photonfile);
load(Photonfile,'L_ph1');

fnstr = 50;   % global feedback parameter to control bump shape
la    = 0.2;  % global parameter to control latency distribution

SaveDataFile = ['LightSkewRes_c180_b2p3_speed205_Fn50La02_', num2str(tim)];
save(SaveDataFile,'Photonfile','rndseed','fnstr');
s_SS = 0;
LightPhoton = 0;
Vol_clamp = -0.07; % clamped voltage in volts
for s = 1:size(SS,1)
    if size(SS{s},2)>3
        if sum(SS{s}(1:4)=='L_ph')==4
            s_SS = s_SS+1;
        end
    end
end

% To solve the problem of 'out of memory', the microvilli is sectioned to
% several big chunks, each chunk with a group of villi (300)
size_cc = 300; %size of one group of microvilli
ph = zeros(100,size_cc);
PH = [];
n_sec = fix(size(L_ph1,2)./size_cc);
n_rem_col = rem(size(L_ph1,2),300);
%this section of programme is to set the input of reminder microvilli to
%void, so the simulation would only simulate the integrate times of size_cc
%(300)microvilli
if n_rem_col~=0
    for i = 1:s_SS
        eval(['load ' Photonfile ' L_ph' num2str(i) ';']);
        eval(['L_ph' num2str(i) '(:,n_sec*size_cc+1:end)=[];']);
    end
end

%for the purpose of saving memory, the number of microvilli are
%sectioned, n_sec is the number of sections.
for kk = 1:n_sec
    kk
    
    %connect the photon input data along time axis together,
    for i = 1:s_SS
        eval(['load ' Photonfile ' L_ph' num2str(i) ';']);
        eval(['ph = L_ph' num2str(i) '(:,(kk-1)*size_cc+1:kk*size_cc);']);
        PH1 = [PH;ph];
        PH = PH1;
        clear PH1;
        eval(['clear L_ph' num2str(i) ';']);
    end
    
    nstep = size(PH,1)*10; %nstep is the count of 0.1ms
    yini =[0    1        0         0         0         0         0   50.0000         1         0];
    yini_online = [yini,zeros(1,5)];
    C_ex = 5000;
    sii = zeros(1,size(PH,2)); %to note down the size of the response
    CCIN_REV=0;
    CCIN_REV1 = 0; % CCIN_REV and CCIN_REV1 are for store the current when TRP reversal potential is changing
    mm = max(sii);
    eval(['BumpS_state' num2str(kk) '=cell(size_cc,1);']);
    %this is to store the state information of jth section of microvilli
    
    for jj = 1:size_cc
        L1 = PH(:,jj);
        if sum(L1)~=0 % to see whether the photon input is zero, if so, do not simulate
            
            clear para
            Initial_PR_Goodaa1new;
            [t1,yy1] = RandomBumpModel_GillespieD_Continue_TRPLa(nstep,yini_online,L1,para,fnstr,la);
            I_m = yy1(:,15); % the macroscopic current including NaCa exchanger current
            
            eval(['BumpS_state' num2str(kk) '{jj,1}(:,1)=yy1(:,1);']);
            eval(['BumpS_state' num2str(kk) '{jj,1}(:,2)=yy1(:,13);']);
            
            sii(jj)=size(yy1,1);
            clear t1 yy1
            
            if length(CCIN_REV)<sii(jj)
                CCIN_REV1 = [CCIN_REV;zeros(sii(jj)-length(CCIN_REV),1)];
                CCIN_REV = CCIN_REV1;
                clear CCIN_REV1
                CCIN_REV = CCIN_REV + I_m;
            elseif length(CCIN_REV)>sii(jj)
                Im1 = [I_m; zeros(length(CCIN_REV)-sii(jj),1)];
                CCIN_REV = CCIN_REV + Im1;
                clear yyc Im1
            elseif length(CCIN_REV)==sii(jj)
                CCIN_REV = CCIN_REV + I_m ;
            end
            
        else
            eval(['BumpS_state' num2str(kk) '{jj,1}=0;']);
        end % end if sum(L1)~=0
        
    end % end for
    
    eval(['MacroSecRev' num2str(kk) '=CCIN_REV;']);
    clear ccin CCIN_REV
    
    save(SaveDataFile,'BumpS_state*', '-append');
    fclose('all');
    clear *state*;
    save(SaveDataFile,'Macro*', '-append');
    clear Macro*
    fclose('all');
    PH=[];
end
Macro_C = Sum_MacroSecFun(SaveDataFile,100); % sum the LIC from 100 sections of microvilli
eval(['Macro_C' num2str(tim) '=Macro_C;']);
save(SaveDataFile,'para','Macro_C*','-append');
Vol_FeedbackCluster % calculate membrane voltage response

