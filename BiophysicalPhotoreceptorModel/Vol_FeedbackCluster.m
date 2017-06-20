% This program calculates the voltage response by including the voltage
% feedback mechanism, which regulates the TRP driving force.

% Zhuoyi Song, updated 06/2017

%Reference: Juusola et al 2017 ELife
%parameters for the HH model for the cell membrane. depending on light
%conditions, the parameters need to change accordingly. 

%param = [-70 4 -53.2 6*0.585e-3 -85 0.5*0.85e-2 0 0 -5 0 5e-3 2e-3 0.3e-3 0]; %Droso BG0
param = [-70 4 -57.1 6*0.585e-3 -85 0.4*0.85e-2 0 0 -5 0 3e-3 0.8e-3 0.11e-3 0]; %Droso BG1
%param = [-70 4 -57.1 1.8*0.585e-3 -85 0.15*0.85e-2 0 0 -5 0 3e-3 0.8e-3 0.11e-3 0]; %Droso BG3

Current = Macro_C;

I =  Current*0.001; % in unit of nA

samprate = 1000;

niter = 100; % iteration of voltage feedback calculation

[Volt Inaca Inak Ileak Inew Ica Icleak Icl Ishaker Ishab Isti]= wt_cc_model_pump(I,param,samprate);
Vol = Volt(1:end,1); % in unit of mV

eval(['VolP_' num2str(tim) '=Vol;']);
save(SaveDataFile, 'VolP_*','-append');

for ii = 1:niter
    ii;
    LICsum = 0;
    for i = 0:99
        j = i+1;
        eval(['load ' SaveDataFile ' BumpS_state' num2str(j) ';']);
        eval(['bss = BumpS_state' num2str(j) ';']);
        Isum = 0;
        clear BumpS_state*
        for k = 1:size(bss,1)
            if size(bss{k,1},2)==2
                BumpS = bss{k,1}(:,1);
                TRPrev = 0.001*20*ones(length(BumpS),1);
                
            else
                BumpS = 0;
                TRPrev = 0.001*20*ones(length(BumpS),1);
            end
            
            TRPc = 8;
            df = TRPrev*1000-Vol;
            drive = max(df,0);
            Iextra = BumpS.* 8.*drive*0.001; %pA
            
            Isum = Isum + Iextra;
        end
        LICsum = LICsum + Isum;
        clear bss
    end
    
    [VoltD InacaD InakD IleakD InewD IcaD IcleakD IclD IshakerD IshabD IstiD]= wt_cc_model_pump(LICsum*0.001,param,samprate);
    Vol = VoltD(:,1);
    if ii == niter-1
        LICsum1 = LICsum;
    elseif ii == niter
        LICsum2 = LICsum;
    end
end
eval(['VolF_' num2str(tim) '=Vol;']);
save(SaveDataFile, 'VolF_*','LICsum*','-append');
