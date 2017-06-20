% This programme calculate the macroscopic LIC 

% The input is SaveDataFile and Nsec. 
% SaveDataFile is a string containing the name of the file that saves the state from 
% the 30,000 microvilli. Nsec is the number of sections of microvilli

% Example: SaveDataFile = 'ST4_flat2s_30k_PR_KCa1_1'
% Macro_C is calculated from Nsec sections of microvilli, with each section
% containing 300 microvilli, Nsec should be less than 100

% Zhuoyi Song, updated in 06/2017

function Macro_C = Sum_MacroSecFun(SaveDataFile,Nsec)
eval(['load ' SaveDataFile ' MacroSecRev' num2str(1) ';']);
Macro_C = 0;bb = size(MacroSecRev1,1);
for h = 1:Nsec
    eval(['load ' SaveDataFile ' MacroSecRev' num2str(h) ';']);
    eval(['mm = MacroSecRev' num2str(h) ';']);
    length(mm);
    if length(mm)>bb
        mm(aa+1:end)=[];
        Macro_C = Macro_C + mm;
    else
         Macro_C = Macro_C + mm;
    end
    clear MacroSecRev*
end
