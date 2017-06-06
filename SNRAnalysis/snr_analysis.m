function [out,resol,ch,tf]=snr_analysis(datafile,chans,samprate,winsize,win,filt,filtV,snr,vc,stim,fwd,fm,gain_cor)
% Calculates the SNR(f) from a suitable datafile; co. Mikko Juusola 1997
%
% This script is a part of a large software package, BIOSYST; co. Mikko Juusola 1997-2017 
% It has been used in many publications, including: Juusola & Hardie, JGP 2001a,b; 
% Juusola & de Polavieja, JGP 2003; Niven et al., Nature 2003; Zheng et al., JGP 2006; 
% Song et al., CB 2012; Wardill et al., Science 2012; Song & Juusola, JNS 2014;
% Juusola et al., eLife 2017.
% 
% INPUT PARAMETERS:
% "datafile" is a matrix of data (individual responses in columns) to repeated stimulation
% "chans" is a number of channels containing the recorded data (eg. light stimulus, voltage response, temp etc.) 
%       When analysing just one block of data (say voltage responses), "chans" should be 1
% "samprate" is the sampling rate; i.e. "samprate" is 1000 for 1 ms sampling (= 1000 samples/s)
% "winsize" is the number of data points used for each FFT calculation
% "win" is the selected windowing function (see lines 61-70 below); e.g. use 1 for BlackmanHarris window
% "filt" is optional filtering (see lines 130-141 below); e.g. 1 = no filtering (default), 2 = median filtering
% "filtV" is an optional filter vector (see lines 130-141 below); the default setting is [3,3]
% "snr" is a BIOSYST tag for tracking the system's recording state; here "snr" should be 1
% "vc" is a BIOSYST tag for tracking the recording mode: 1 is c-clamp, 2 v-clamp; here "snr" should be 1
% "stim" is the BIOSYST default for the stimulus being on channel 2 ("stim=2") 
% "fwd" defaut here is for forward correlation ("fwd=1") 
% "fm" is your selection for the output type ("out"); the default is SNR "fm=1" (see lines 154-161 below)
% "gain_cor" here should be "gain_cor=0" 
%
% OUTPUT PARAMETERS:
% "out" is the requested analytical result; e.g. for "fm" 1, "out" is SNR(f) 
% "resol" is the frequency resolution of the requested analytical result 
% "ch" is a BIOSYST tag for the number of channels (subfigures) used for plotting the result
% "tf" is a BIOSYST tag for indicating whether to use time or frequency axis in the plots

% HOW TO TEST THIS FUNCTION ON YOUR DATA?
% Say, you have a matrix (data) of ten responses in columns, each 2000 points long.
% Your data was collected with 1 ms sampling (samprate=1000).
% You want to use 500 point BlackmanHarris windowing (winsize=500; win=1).
% You do not want spectral filtering (filt=1; filtV=0).
%           But for medial filtering, you could try for example (filt=2; filtV=[3,3]).
% Your data is in the correct (columnal) format (snr=1).
% You have data from just one input channel (vc=1).
% The default is responses from channel 1 and stimulus from channel 2 (stim=2).
% The default is SNR (fm=1) from forward correlation (fwd=1) over frequency without gain correction (gain-cor=0).
% These choices give you these parameter settings:
% [out,resol,ch,tf]=snr_analysis(data,1,1000,500,1,1,0,1,1,2,1,1,0)
% 
% When plotting the result in Matlab, use: plot(resol,out)
%
% ############################################## initialize ##########################################
if snr ==0, disp('This data is not suitable for SNR-analysis!!!'); [out]=[]; [resol]=[]; [ch]=chans; tf=1;
else 

   % ********************************** Defining the number of overlaps *********************************
    dur=size(datafile,1);                                               % dur = file length in points
    no_runs = size(datafile,2)/chans;                                   % the number of runs
    if winsize*2 < dur, slides=round(dur/(winsize/2))-1;                % in the normal case, finds
        while slides*winsize/2+winsize > dur, slides = slides-1; end;   % the number of 50% overlaps  
    else winsize=dur/2; slides = 2;                                     % if not, halves the window
    end, fprintf('%d overlapping %d point windows for a FFT of a %d point datafile!!',...
       slides+1, winsize, dur);

    %********************************** Preparing the WINDOW functions **********************************   
    if win == 1,window=blackmanharris(winsize); end % BlackmanHarris windowing of the data
    if win == 2,window=blackman(winsize); end 		% Blackman windowing of the data
    if win == 3,window=hanning(winsize); end 		% Hanning windowing of the data
    if win == 4,window=hamming(winsize); end 		% Hamming windowing of the data
    if win == 5,window=bohmanwin(winsize); end 		% Bohman windowing of the data
    if win == 6,window=nuttallwin(winsize); end 	% Nuttall windowing of the data
    if win == 7,window=gausswin(winsize); end 		% Gaussian windowing of the data
    if win == 8,window=bartlett(winsize); end 		% Bartlett windowing of the data
    if win == 9,window=barthannwin(winsize); end 	% Barthann windowing of the data
    if win == 10,window=boxcar(winsize); end 		% No windowing of the data
   
    %********************************** Defining the key parameters for the calculations *****************  
    nyquist=1/2;    n_size=winsize;	freq=(1:n_size/2)/(n_size/2)*nyquist;	
    SIGNAL_PWR = []; 	XY_SIGNAL_PWR =[];  	NOISE_PWR=[];	SNR=[];

    %********************************** Calculating the average responses ********************************
    for n=1:chans, noise_var{1,n}=[]; signal_var{1,n}=[]; signal_mean{1,n}=[]; end; 
    OPUT=zeros(dur,chans); pos=1;
    for i = 1:no_runs, 
        for n=1:chans, oput=datafile(:,pos); OPUT(:,n)=OPUT(:,n)+oput; pos=pos+1; end;
    end; SIGNAL=OPUT/no_runs;
        
    %*************************** Calculating the average noise power spectra *****************************        
    sumXX=zeros(winsize,chans); XX=zeros(winsize,chans); NOISE_PWR=zeros(winsize,chans);  
    %------------------------------- The loop for segmenting the data and taking fft ---------------------
    for b = 1:no_runs
        for n=1:chans, 
        noise_pwr(:,n)=zeros(winsize,1); sumXX(:,n)=zeros(winsize,1);     
        noise = datafile(:,chans*b-(chans-n))-SIGNAL(:,n); % selects the noise trace for FFT   
            for i=0:slides                         % splits the noise trace to window-sized samples
                pos=i*(winsize/2)+1;
           	    data=noise(pos:pos+(winsize-1)); 
                data=data-mean(data);	   	
                noise_var{1,n}=[noise_var{1,n}, var(data)];
                data=data.*window;		      % windows each sample
                x=fft(data);
                XX=conj(x).*x;			  % calculates the power spect of the segment
                sumXX(:,n)=sumXX(:,n)+XX;   
            end
        end; 
        windowcorrect=mean(window.^2);                              % power spectum scaling
        noise_pwr=sumXX/slides/winsize/samprate*2/windowcorrect;    % power spectum scaling
        NOISE_PWR=NOISE_PWR+noise_pwr; 
    end; NOISE_PWR=NOISE_PWR/no_runs;                               % the mean noise power spectrum
    
    % ****************************** Calculating the signal power spectra ********************************    
    sumXX=zeros(winsize,chans); sumXY= zeros(winsize,1); SIGNAL_PWR=zeros(winsize,chans);
    XX=zeros(winsize,chans); XY=zeros(winsize,chans);
    %------------------------------- The loop for segmenting the data and taking fft ---------------------
    for i=0:slides, pos=i*(winsize/2)+1; spect=[];
        for n=1:chans, signal=SIGNAL(:,n);
            data=signal(pos:pos+(winsize-1)); 		% Takes i overlaping samples of winsize		
            signal_mean{1,n}=[signal_mean{1,n}, mean(data)];
            data=data-mean(data);	   						% data=channel data, mean removed	    
            signal_var{1,n}=[signal_var{1,n}, var(data)];
            data=data.*window;		      					% windows each sample					
            spect(:,n)=fft(data); 						    % takes the FFT out of it				
            XX(:,n)=conj(spect(:,n)).*spect(:,n);
            sumXX(:,n)=sumXX(:,n)+XX(:,n);                              % in fwd correlation        
        end;                                                            % XY=conj(stim).*spect(resp)
        if chans ==1, XY=XX; sumXY=sumXX;                               % when data from one channel 
        else
            if fwd==1, XY=conj(spect(:,stim)).*spect(:,vc); sumXY=sumXY+XY; end; % in rev corr.  
            if fwd==0, XY=conj(spect(:,vc)).*spect(:,stim); sumXY=sumXY+XY; end; % vice versa    
        end
    end; 
    windowcorrect=mean(window.^2);                                  % power spectum scaling
    SIGNAL_PWR=sumXX/slides/winsize/samprate*2/windowcorrect;       % the mean signal power spectrum
    XY_SIGNAL_PWR=sumXY/slides/winsize/samprate*2/windowcorrect;    % the mean cross spectrum
    
    %****************************** OPTIONAL SPECTRAL FILTERING ************************************                                          
    if filt(1) == 2,                                    % 1 =no filtering, 2 = median filtering
        for n=1:size(SIGNAL_PWR,2), SIGNAL_PWR(:,n)=medfilt1(SIGNAL_PWR(:,n),filtV(1)); end
        for n=1:size(XY_SIGNAL_PWR,2), XY_SIGNAL_PWR(:,n)=medfilt1(XY_SIGNAL_PWR(:,n),filtV(1)); end
    elseif filt(1) == 3, b=ones(1,filtV(1))/filtV(1);	% 2 =FFTFILT filtering
        for n=1:size(SIGNAL_PWR,2),	SIGNAL_PWR(:,n)=FFTFILT(b,SIGNAL_PWR(:,n)); end
        for n=1:size(XY_SIGNAL_PWR,2), XY_SIGNAL_PWR(:,n)=FFTFILT(b,XY_SIGNAL_PWR(:,n)); end
    elseif filt(1) == 4, b=ones(1,filtV(1))/filtV(1); 	% 3 =FILTFILT filtering
        for n=1:size(SIGNAL_PWR,2), SIGNAL_PWR(:,n)=FILTFILT(b,1,SIGNAL_PWR(:,n)); end
        for n=1:size(XY_SIGNAL_PWR,2), XY_SIGNAL_PWR(:,n)=FILTFILT(b,1,XY_SIGNAL_PWR(:,n)); end
    end

    % ****************************** Calculates SNR(f) ******************************************
    if gain_cor==1,  % the response pwr corrected by the stimulus pwr (not for single dataset use) 
        gain_norm=(SIGNAL_PWR(:,stim)/mean(SIGNAL_PWR(1:100, stim))); 
        SIGNAL_PWR(:,vc)=SIGNAL_PWR(:,vc)./gain_norm;
    end
    for n=1:chans, SNR(:,n) = SIGNAL_PWR(:,n)./NOISE_PWR(:,n); end % calculates the SNR(f) 
        x_parameters=[mean(signal_mean{1,1}), mean(signal_var{1,1}),...
        length(signal_mean{1,1}), mean(noise_var{1,1}),...
        length(noise_var{1,1}), mean(signal_var{1,1})/mean(noise_var{1,1})];
        disp('signal mean, variance & n; noise variance & n; snr(t)'), num2str(x_parameters)
    
    switch fm,
        case 1, [out]=SNR;            % SNR
        case 2, [out]=log2(SNR+1);    % Information
        case 3, [out]=SIGNAL_PWR;     % Signal Power Spectrum
        case 4, [out]=NOISE_PWR;      % Noise power Spectrum
        case 5, [out]=SNR./(SNR+1);     % Coherence from SNR
        case 6, [out]=log2(1./(1-(SNR./(SNR+1))));  % Info from Coherence    
    end
    [resol]=(freq*samprate)'; out=out(1:winsize/2,:);    % set the plot frequency resolution !!
    [ch]=chans; [tf]=2;    
end % ends the analysis if snr == 1 

