%% This code is for analyzing EEG data obtained from HUMANE AI project

% Contents :
% Epoching data using event markers
% Comparison of short time power spectrum
% ERD/ERS computation with Power method and topographical mapping
% ERD/ERS computation with Intertrial Variance method and topographical mappin

% Written by Zal Aktaş , June 2023 
% during internship at TUBITAK BILGEM B3LAB


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


%% Epoch the data based on trials
epstart = -0.5;
epfinish = 2.5;

mdata = readmatrix("Psychopy_1_2022-11-28T134251.192511_EPOCFLEX_6649_2022.11.28T13.42.51+03.00_intervalMarker.csv");

eptidx(:,1) = dsearchn(time',mdata(:,1)+epstart);
eptidx(:,2) = dsearchn(time',mdata(:,1)+epfinish);

epochdata = zeros(length(mdata),32,385);
for epi = 1:length(mdata)
    for chi = 1:32
        epochdata(epi,chi,:) = eegdata(chi,eptidx(epi,1):eptidx(epi,2));
    end
end
timeInt = linspace(epstart,epfinish,size(epochdata,3));

load eventtitles.mat
eventtitles = strrep(eventtitles, '_' , ' ');

[epochNames, epstartidx] = unique(eventtitles,'stable');
epochIntervals = [epstartidx [epstartidx(2:end); epstartidx(end)+1 ]-1 ];

%% Determine the subject spesific frequency band by comparing the short time power spectra

epoch2test = 333;
chan2test = 28;

%define reference and activity period
refInt = [-0.5 -0.2];
actInt = [0.8 1.1];

refIdx = dsearchn(timeInt', refInt');
actIdx = dsearchn(timeInt',actInt');

refData = epochdata(:,:,refIdx(1):refIdx(2));
actData = epochdata(:,:,actIdx(1):actIdx(2));

%calculate the power spectrum for reference and activity periods

actnpnts = size(actData,3);
acthz = linspace(0,srate/2,floor(actnpnts/2)+1);
actpow = 2*abs(fft(squeeze(actData(epoch2test,chan2test,:)))/actnpnts).^2;

refnpnts = size(refData,3);
refhz = linspace(0,srate/2,floor(refnpnts/2)+1);
refpow = 2*abs(fft(squeeze(refData(epoch2test,chan2test,:)))/refnpnts).^2;

plot(refhz,10*log10(refpow(1:length(refhz))));
hold on;
plot(acthz,10*log10(actpow(1:length(acthz))));
legend(['reference period (' , num2str(refInt(1)) , ' to ' , num2str(refInt(2)) , ' sec)'],... 
    ['activity period (' , num2str(actInt(1)) , ' to ' , num2str(actInt(2)) , ' sec)']);
xlabel('Frequency (Hz)'); ylabel('Band power (dB)'); xlim([0.5 32]); 
title("Analysis to determine the frequency band for ERD power analysis" + ...
"\newline Channel " + num2str(chan2test) + ": " +chanlocs(chan2test) +...
"  Epoch " + num2str(epoch2test) + ": " + eventtitles(epoch2test))

%% ERD / ERS Analysis 

%% Power Method
ch2plot = 7;
epseg2plot = 8;
freqband = [18 22];

[epochNames, epstartidx] = unique(eventtitles,'stable');
epochIntervals = [epstartidx [epstartidx(2:end); epstartidx(end)+1 ]-1 ];
ep2plot = epochIntervals(epseg2plot,1);

subplot(511)
plot(timeInt,squeeze(epochdata(ep2plot,ch2plot,:)))
title(['Raw EEG data (channel ',num2str(ch2plot),', trial ', num2str(ep2plot),')'])
xline(0,'--r'); xlabel('time (sec)'); ylabel('Amplitude (V)')

%1. bandpass filter with 
erdFilt = designfilt('bandpassfir','CutoffFrequency1',freqband(1),'CutoffFrequency2',freqband(2),'SampleRate',srate,'FilterOrder',128);

erdData = zeros(size(epochdata));
for epi = 1:size(epochdata,1)
    for ci = 1:size(epochdata,2)
    erdData(epi,ci,:) = filtfilt(erdFilt,squeeze(epochdata(epi,ci,:)));
    end
end
subplot(512)
plot(timeInt,squeeze(erdData(ep2plot,ch2plot,:)))
title(['Bandpass filtering 8-12 Hz (channel ',num2str(ch2plot),', trial ', num2str(ep2plot),')'])
xline(0,'--r'); xlabel('time (sec)'); ylabel('Amplitude (V)')

%2. squaring of the amplitude samples to obtain power samples
erdData = erdData.^2;
subplot(513)
plot(timeInt,squeeze(erdData(ep2plot,ch2plot,:)))
title(['Squaring (channel ',num2str(ch2plot),', trial ', num2str(ep2plot),')'])
xline(0,'--r'); xlabel('time (sec)'); ylabel('Power (V^2)')

%3. averaging of power samples across all trials with same event

avData = zeros(length(epochIntervals),size(epochdata,2),size(epochdata,3));
for ei = 1: length(epochIntervals)
    avData(ei,:,:) = squeeze(mean(erdData(epochIntervals(ei,1):epochIntervals(ei,2),:,:),1));
end

subplot(514)
plot(timeInt,squeeze(avData(epseg2plot,ch2plot,:)))
title(['Averaging over trials (channel ',num2str(ch2plot),')'])
xline(0,'--r'); xlabel('time (sec)'); ylabel('Power (V^2)')

%4. averaging over time samples to smooth the data and reduce the variability.
timeWin = 0.2;
winIdx = ceil(timeWin*length(timeInt)/3);
avData = movmean(avData,winIdx,3);

subplot(515)
plot(timeInt,squeeze(avData(epseg2plot,ch2plot,:)))
title(['Averaging over time samples to smooth the data (channel ',num2str(ch2plot),')'])
xline(0,'--r'); xlabel('time (sec)'); ylabel('Power (V^2)')

%5. averaging reference interval and ERD/ERS calculation
R = mean(avData(:,:,refIdx(1):refIdx(2)),3);
ERD = 100*(avData-R)./R;

%Plot ERD/ERS
figure()
plot(timeInt,squeeze(ERD(epseg2plot,ch2plot,:)))
title("ERD/ERS of Epoch: " + epochNames(epseg2plot) +...
    " From channel " + num2str(ch2plot) + ': ' +chanlocs(ch2plot) + ...
    "\newline Relative power change to baseline in the frequency band " +num2str(freqband(1)) + "-" +num2str(freqband(2)) + " Hz")
yline(0,'--r');xline(0,'--r');
xlabel('time (sec)'); ylabel('Relative power (%)');

%% topoplot ERD

epseg2plot =10;

% time points for topographies
times2plot = -0.5:0.1:2.5; % in sec
tidx = dsearchn(timeInt',times2plot');

% define subplot geometry
subgeomR = ceil(sqrt(length(tidx)));
subgeomC = ceil(length(tidx)/subgeomR);

figure()
t = tiledlayout( subgeomR,subgeomC);
for i=1:length(times2plot)
    nexttile;
    topoplot(squeeze(ERD(epseg2plot,:,tidx(i))),EEG.chanlocs,'electrodes','off' , 'shading','interp');
    set(gca,'clim',[-max(abs(ERD(epseg2plot,:,:)),[],"all") max(abs(ERD(epseg2plot,:,:)),[],"all")])
    title([ num2str(times2plot(i)) ' sec' ]);
end
cb = colorbar;cb.Layout.Tile = 'east';
cb.Ticks = [-max(abs(ERD(epseg2plot,:,:)),[],"all") max(abs(ERD(epseg2plot,:,:)),[],"all")];
cb.TickLabels = {'ERD','ERS',};cb.FontSize =20;
title(t,"ERD/ERS Analysis for " + epochNames(epseg2plot) + "\newline " + num2str(freqband(1)) + "-" + num2str(freqband(2)) + " Hz")

%% I have not used the results of IV method but code below is written for future use

%% Intertrial variance method

ch2plot = 7;
epseg2plot = 8;
freqband = [8 12];

[epochNames, epstartidx] = unique(eventtitles,'stable');
epochIntervals = [epstartidx [epstartidx(2:end); epstartidx(end)+1 ]-1 ];
ep2plot = epochIntervals(epseg2plot,1);

subplot(311)
plot(timeInt,squeeze(epochdata(ep2plot,ch2plot,:)))
title(['Raw EEG data (channel ',num2str(ch2plot),', trial ', num2str(ep2plot),')'])
xline(0,'--r'); xlabel('time (sec)'); ylabel('Amplitude (V)')

%1. bandpass filter with 8-12 Hz
erdFilt = designfilt('bandpassfir','CutoffFrequency1',freqband(1),'CutoffFrequency2',freqband(2),'SampleRate',srate,'FilterOrder',128);

erdData = zeros(size(epochdata));
for epi = 1:size(epochdata,1)
    for ci = 1:size(epochdata,2)
    erdData(epi,ci,:) = filtfilt(erdFilt,squeeze(epochdata(epi,ci,:)));
    end
end
subplot(312)
plot(timeInt,squeeze(erdData(ep2plot,ch2plot,:)))
title(['Bandpass filtering ', num2str(freqband(1)),'-',num2str(freqband(2)) ,' Hz (channel ',num2str(ch2plot),', trial ', num2str(ep2plot),')'])
xline(0,'--r'); xlabel('time (sec)'); ylabel('Amplitude (V)')

%2. calculation of the point-to-point intertrial variance

%calculate trial averages
avData = zeros(length(epochIntervals),size(epochdata,2),size(epochdata,3));
for ei = 1: length(epochIntervals)
    avData(ei,:,:) = squeeze(mean(erdData(epochIntervals(ei,1):epochIntervals(ei,2),:,:),1));
end

%substract each point from trial average and take the squared mean of the
%difference
IV = zeros(size(avData));
for ai = 1:size(avData,1)
    for ii = epochIntervals(ai,1):epochIntervals(ai,2)
        IV(ai,:,:) = ((epochdata(ii,:,:) - avData(ai,:,:)).^2)/length(epochIntervals(ai,1):epochIntervals(ai,2));
    end
end

subplot(313)
plot(timeInt,squeeze(IV(epseg2plot,ch2plot,:)))
title("Intertrial varience (channel "+num2str(ch2plot)+", "+epochNames(epseg2plot)+")")
xline(0,'--r'); xlabel('time (sec)'); ylabel('IV')

%3. averaging over time
timeWin = 0.15;
winIdx = ceil(timeWin*length(timeInt)/3);
IV = movmean(IV,winIdx,3);

%4. averaging reference interval and ERD/ERS calculation
R = mean(IV(:,:,refIdx(1):refIdx(2)),3);
ERD = 100*(IV-R)./R;

%Plot ERD/ERS
figure()
plot(timeInt,squeeze(ERD(epseg2plot,ch2plot,:)))
title("ERD/ERS of Epoch: " + epochNames(epseg2plot) +...
    " From channel " + num2str(ch2plot) + ': ' +chanlocs(ch2plot) + ...
    "\newline Non-phase locked activity in the frequency band " +num2str(freqband(1)) + "-" +num2str(freqband(2)) + " Hz")
yline(0,'--r');xline(0,'--r');
xlabel('time (sec)'); ylabel('Relative IV (%)');

%% Plot topograph

epseg2plot =16;

% time points for topographies
times2plot = -0.5:0.1:2.5; % in sec
tidx = dsearchn(timeInt',times2plot');

% define subplot geometry
subgeomR = ceil(sqrt(length(tidx)));
subgeomC = ceil(length(tidx)/subgeomR);

figure()
t = tiledlayout( subgeomR,subgeomC);
for i=1:length(times2plot)
    nexttile;
    topoplot(squeeze(-ERD(epseg2plot,:,tidx(i))),EEG.chanlocs,'electrodes','off' , 'shading','interp');
    set(gca,'clim',[-max(abs(ERD(epseg2plot,:,:)),[],"all") max(abs(ERD(epseg2plot,:,:)),[],"all")])
    title([ num2str(times2plot(i)) ' sec' ]);
end
cb = colorbar;cb.Layout.Tile = 'east';
cb.Ticks = [-max(abs(ERD(epseg2plot,:,:)),[],"all") max(abs(ERD(epseg2plot,:,:)),[],"all")];cb.TickLabels = {'ERD','ERS',};cb.FontSize =20;
title(t,"ERD/ERS (IV) Analysis for " + epochNames(epseg2plot))