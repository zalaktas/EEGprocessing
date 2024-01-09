%% This code is for analyzing EEG data obtained from HUMANE AI project

% Contents :
% Epoching data using event markers
% Time-frequency analysis
% Erp analysis and topographical mapping

% Written by Zal Aktaş , June 2023 
% during internship at TUBITAK BILGEM B3LAB

% Acknowladgements :
% Some parts of the code is directly taken from :
% Michael X Cohen 
% From UDEMY Course : 
% Complete neural signal processing and analysis: Zero to hero

%% Load the data and spesify parameters
% this is for getting topoplot function working
eeglab 
% if you get the error "'eeglab' is not found in the current folder"
% simply add its folder to the MATLAB path.

close all
clear all
clc
load ICAprocessedEEG.mat


srate = EEG.srate;
time = EEG.times/1000; %in sec

%reduce data to 32 relevant channels

eegdata = EEG.data;
npnts = length(eegdata);

%load the channel location data
arr={EEG.chanlocs.labels};
chanlocs=[];
for x=arr
    chanlocs=[chanlocs; string(x{1})];
end


%% Plot the temporal and spectral characteristics of the data

chan2plot = 1;

% plot the time-domain signal
figure(1)
subplot(211)
plot(time,eegdata(chan2plot,:))
xlabel('Time (s)'), ylabel('Voltage (\muV)') , 
title(['Filtered EEG data from channel ' , num2str(chan2plot) , ' in time domain'])
zoom on


% static spectral analysis
hz = linspace(0,srate/2,floor(npnts/2)+1);
ampl = 2*abs(fft(eegdata(chan2plot,:))/npnts);

subplot(212)
plot(hz,ampl(1:length(hz)),'r','linew',2)

xlabel('Frequency (Hz)') ,ylabel('Amplitude')
title(['Filtered EEG data from channel ' , num2str(chan2plot) , ' in frequency domain'])

%% Epoch the data based on event markers
epstart = -0.5;
epfinish = 2.5;

mdata = readmatrix("Psychopy_1_2022-11-28T134251.192511_EPOCFLEX_6649_2022.11.28T13.42.51+03.00_intervalMarker.csv");

eptidx(:,1) = dsearchn(time',mdata(:,1)+epstart);
eptidx(:,2) = dsearchn(time',mdata(:,1)+epfinish);

epochtimeintervals(:,1) = time(eptidx(:,1));
epochtimeintervals(:,2) = time(eptidx(:,2));

epochdata = zeros(length(mdata),32,385);
for epi = 1:length(mdata)
    for chi = 1:32
        epochdata(epi,chi,:) = eegdata(chi,eptidx(epi,1):eptidx(epi,2));
    end
end
%% Plot indivual epochs
load eventtitles.mat
eventtitles = strrep(eventtitles, '_' , ' ');

ep2plot = 135;
chan2plot = 7;

figure(3)
subplot(211)
plot(time(eptidx(ep2plot,1):eptidx(ep2plot,2)) - mdata(ep2plot,1) ,squeeze(epochdata(ep2plot,chan2plot,:)));
title("Epoch " + num2str(ep2plot) + ": " + eventtitles(ep2plot) + ...
    "\newline Filtered EEG data in time domain from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
xlabel('Time (s)'); ylabel('Voltage (a.u)');
xline(0,'--r');set(gca,'xlim',[-0.5 2.5]);

% static spectral analysis
npnts = eptidx(ep2plot,2) - eptidx(ep2plot,1);
hz = linspace(0,srate/2,floor(npnts/2)+1);
ampl = 2*abs(fft(squeeze(epochdata(ep2plot,chan2plot,:)))/npnts);

subplot(212)
plot(hz,ampl(1:length(hz)));
title("Epoch " + num2str(ep2plot) + ": " + eventtitles(ep2plot) + ...
    "\newline Filtered EEG data in frequency domain from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
xlabel('Frequency (Hz)'); ylabel('Amplitude (a.u)');



%% Time-frequency analysis

% soft-coded parameters
freqrange  = [5 32]; % extract only these frequencies (in Hz)
numfrex    = 50;       % number of frequencies between lowest and highest
whichEpoch = 333;
baseline_window = [ -0.4 -0.2]; 
baseidx = dsearchn(time',baseline_window' + 1 + time(eptidx(whichEpoch,1))) - eptidx(whichEpoch,1) +1;

% set up convolution parameters
wavtime = (0:2*srate)/srate;
wavtime = wavtime - mean(wavtime);
frex    = linspace(freqrange(1),freqrange(2),numfrex);
nData   = length(eegdata(eptidx(whichEpoch,1):eptidx(whichEpoch,2)));
nKern   = length(wavtime);
nConv   = nData + nKern - 1;
halfwav = (length(wavtime)-1)/2;

% number of cycles
numcyc = linspace(3,15,numfrex);


% create wavelets
cmwX = zeros(numfrex,nConv);
for fi=1:numfrex
    
    % create time-domain wavelet
    twoSsquared = 2 * (numcyc(fi)/(2*pi*frex(fi))) ^ 2;
    cmw = exp(2*1i*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) / twoSsquared );
    
    % compute fourier coefficients of wavelet and normalize
    cmwX(fi,:) = fft(cmw,nConv);
    cmwX(fi,:) = cmwX(fi,:) ./ max(cmwX(fi,:));
end

% initialize time-frequency output matrix
tf = zeros(32,numfrex,nData);

% loop over channels
for chani=1:32
    
    % compute Fourier coefficients of EEG data (doesn't change over frequency!)
    eegX = fft( eegdata(chani,eptidx(whichEpoch,1):eptidx(whichEpoch,2)) ,nConv,2);
    
    % loop over frequencies
    for fi=1:numfrex
        
        % second and third steps of convolution
        as = ifft( cmwX(fi,:).*eegX ,nConv );
        
        % cut wavelet back to size of data
        as = as(halfwav+1:end-halfwav);
        
        % extract power and phase
        tf(chani,fi,:) = abs(as).^2;        
    end % end frequency loop
end % end channel loop

%baseline normalize
tfDB = 10*log10(bsxfun(@rdivide, tf, mean(tf(:,:,baseidx(1):baseidx(2)),3)));

%% plotting time-frequency results

chan2plot = 28;

figure(50)
subplot(211)
plot(time(eptidx(whichEpoch,1):eptidx(whichEpoch,2)) - mdata(whichEpoch,1) ,squeeze(epochdata(whichEpoch,chan2plot,:)));
title("Epoch " + num2str(whichEpoch) + ": " + eventtitles(whichEpoch) + ...
    "\newline Filtered EEG data in time domain from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
xlabel('Time (s)'); ylabel('Voltage (a.u)');
xline(0,'--r');set(gca,'xlim',[-0.5 2.5]);


subplot(212)
contourf(time(eptidx(whichEpoch,1):eptidx(whichEpoch,2)) - mdata(whichEpoch,1) ,frex,squeeze(tfDB(chan2plot,:,:)),40,'linecolor','none')
xlabel('Time (s)'); ylabel('Frequencies (Hz)');
title("Epoch " + num2str(whichEpoch) + ": " + eventtitles(whichEpoch) + ...
    "\newline Baseline Normalized Time-frequency analysis from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
xline(0,'--r','LineWidth',1.5);
set(gca,'clim',[-max(clim) max(clim)] , 'xlim' , [-0.5 2.5]); colorbar
a = colorbar; a.Label.String = 'Power (dB)';
%% make an animation of all epochs

%  for epi = 1:size(epochdata)
%     % soft-coded parameters
%     freqrange  = [5 35]; % extract only these frequencies (in Hz)
%     numfrex    = 50;       % number of frequencies between lowest and highest
%     whichEpoch = epi;
%     baseline_window = [ -0.4 -0.2];
%     baseidx = dsearchn(time',baseline_window' + 1 + time(eptidx(whichEpoch,1))) - eptidx(whichEpoch,1) +1;
% 
%     % set up convolution parameters
%     wavtime = -2:1/srate:2;
%     frex    = linspace(freqrange(1),freqrange(2),numfrex);
%     nData   = length(eegdata(eptidx(whichEpoch,1):eptidx(whichEpoch,2)));
%     nKern   = length(wavtime);
%     nConv   = nData + nKern - 1;
%     halfwav = (length(wavtime)-1)/2;
% 
%     % number of cycles
%     numcyc = linspace(3,15,numfrex);
% 
% 
%     % create wavelets
%     cmwX = zeros(numfrex,nConv);
%     for fi=1:numfrex
% 
%         % create time-domain wavelet
%         twoSsquared = 2 * (numcyc(fi)/(2*pi*frex(fi))) ^ 2;
%         cmw = exp(2*1i*pi*frex(fi).*wavtime) .* exp( (-wavtime.^2) / twoSsquared );
% 
%         % compute fourier coefficients of wavelet and normalize
%         cmwX(fi,:) = fft(cmw,nConv);
%         cmwX(fi,:) = cmwX(fi,:) ./ max(cmwX(fi,:));
%     end
% 
%     % initialize time-frequency output matrix
%     tf = zeros(32,numfrex,nData);
% 
%     % loop over channels
%     for chani=1:32
% 
%         % compute Fourier coefficients of EEG data (doesn't change over frequency!)
%         eegX = fft( eegdata(chani,eptidx(whichEpoch,1):eptidx(whichEpoch,2)) ,nConv);
% 
%         % loop over frequencies
%         for fi=1:numfrex
% 
%             % second and third steps of convolution
%             as = ifft( cmwX(fi,:).*eegX ,nConv );
% 
%             % cut wavelet back to size of data
%             as = as(halfwav+1:end-halfwav);
% 
%             % extract power and phase
%             tf(chani,fi,:) = abs(as).^2;        
%         end % end frequency loop
%     end % end channel loop
% 
%     %baseline normalize
%     tfDB = 10*log10(bsxfun(@rdivide, tf, mean(tf(:,:,baseidx(1):baseidx(2)),3)));
% 
%     chan2plot = 1;
% 
%     figure(4)
%     subplot(211)
%     contourf(time(eptidx(whichEpoch,1):eptidx(whichEpoch,2)) - mdata(whichEpoch,1),frex,squeeze(tf(chan2plot,:,:)),40,'linecolor','none')
%     xlabel('Time (s)'); ylabel('Frequencies (Hz)');
%     title("Epoch " + num2str(whichEpoch) + ": " + eventtitles(whichEpoch) + ...
%         "\newline Time-frequency analysis from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
%     xline(0,'--r','LineWidth',1.5);
%     set(gca,'clim',clim , 'xlim' , [-0.5 2.5]); colorbar
% 
%     subplot(212)
%     contourf(time(eptidx(whichEpoch,1):eptidx(whichEpoch,2)) - mdata(whichEpoch,1) ,frex,squeeze(tfDB(chan2plot,:,:)),40,'linecolor','none')
%     xlabel('Time (s)'); ylabel('Frequencies (Hz)');
%     title("Epoch " + num2str(whichEpoch) + ": " + eventtitles(whichEpoch) + ...
%         "\newline Baseline Normalized Time-frequency analysis from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
%     xline(0,'--r','LineWidth',1.5);
%     set(gca,'clim',[-max(clim) max(clim)] , 'xlim' , [-0.5 2.5]); colorbar
% end

%% erp analysis

%% Average epochs with same event types

[epochNames, epstartidx] = unique(eventtitles,'stable');
epochIntervals = [epstartidx [epstartidx(2:end); epstartidx(end)+1 ]-1 ];

epochErp = zeros(length(epochNames),size(epochdata ,2),size(epochdata ,3));

for edfi = 1:length(epochNames)
    epochErp(edfi,:,:) = squeeze(mean(epochdata(epochIntervals(edfi,1):epochIntervals(edfi,2),:,:),1));
end

%first take the fourier transform of each epoch, then average for stable
%spectral analysis

% calculate power spectrum for each epoch
npnts = size(epochdata,3);
hz = linspace(0,srate/2,floor(npnts/2)+1);
epochPow = zeros(size(epochdata,1),size(epochdata,2),size(epochdata ,3));

for epi = 1:size(epochdata,1)
    for chi = 1:size(epochdata,2)
    epochPow(epi,chi,:) = (2*abs(fft(squeeze(epochdata(epi,chi,:)))/npnts)).^2;
    end
end

% average over same events
epochPowErp = zeros(length(epochNames),size(epochdata ,2),size(epochdata ,3));
for edfi = 1:length(epochNames)
    epochPowErp(edfi,:,:) = squeeze(mean(epochPow(epochIntervals(edfi,1):epochIntervals(edfi,2),:,:),1));
end

%% Plot erp's of epochs
timeInt = linspace(-0.5,2.5,length(epochErp));
epoch2plot = 9;
chan2plot = 1;

figure(6)
subplot(211)
plot(timeInt,squeeze(epochErp(epoch2plot,chan2plot,:)))
title("ERP of Epoch: " + epochNames(epoch2plot) + ...
    "\newline Filtered EEG data in time domain from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
xlabel('Time (s)'); ylabel('Voltage (a.u)');
xline(0,'--r');set(gca,'xlim',[-0.5 2.5]);

subplot(212)
plot(hz,squeeze(epochPowErp(epoch2plot,chan2plot,1:length(hz))))
title("Averaged power spectrum of Epoch: " + epochNames(epoch2plot) + ...
    "\newline Filtered EEG data in frequency domain from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
xlabel('Frequency (Hz)'); ylabel('Power (a.u)');

%% Real lifting vs Imaginary Lifting

imepochs2plot = [12 14 16 ];
reepochs2plot = [11 13 15 ];
chan2plot = 28;

figure(7)
for i = 1:3

    subplot(3,2,2*i-1)
    plot(timeInt,squeeze(epochErp(imepochs2plot(i),chan2plot,:)))
    title("ERP of Epoch: " + epochNames(imepochs2plot(i)) + ...
        "\newline Filtered EEG data in time domain from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
    xlabel('Time (s)'); ylabel('Voltage (a.u)');
    xline(0,'--r');set(gca,'xlim',[-0.5 2.5]);
    
    subplot(3,2,2*i)
    plot(timeInt,squeeze(epochErp(reepochs2plot(i),chan2plot,:)))
    title("ERP of Epoch: " + epochNames(reepochs2plot(i)) + ...
        "\newline Filtered EEG data in time domain from channel " + num2str(chan2plot) + ": " +chanlocs(chan2plot))
    xlabel('Time (s)'); ylabel('Voltage (a.u)');
    xline(0,'--r');set(gca,'xlim',[-0.5 2.5]);
end

%% Topographical analysis of erp

epoch2plot =9;

% time points for topographies
times2plot = -0.5:0.1:2.5; % in sec
tidx = dsearchn(timeInt',times2plot');

% define subplot geometry
subgeomR = ceil(sqrt(length(tidx)));
subgeomC = ceil(length(tidx)/subgeomR);


figure(8)
t = tiledlayout( subgeomR,subgeomC);
for i=1:length(times2plot)
    nexttile;
    topoplot(squeeze(epochErp(epoch2plot,:,tidx(i))),EEG.chanlocs,'electrodes','off' , 'shading','interp');
    set(gca,'clim',[-max(abs(epochErp(epoch2plot,:,:)),[],"all") max(abs(epochErp(epoch2plot,:,:)),[],"all")])
    title([ num2str(times2plot(i)) ' sec' ]);
end

title(t,"Topographical mapping of " + epochNames(epoch2plot) + " around time zero")
