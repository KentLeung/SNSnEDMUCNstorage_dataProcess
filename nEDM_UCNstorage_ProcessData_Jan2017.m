clc; clear all; set(0,'defaulttextinterpreter','none'); linecolors = lines(20); set(0,'DefaultAxesFontSize', 14); set(0,'DefaultTextFontSize', 14);
close all;
%---Files to load---
dataPath  = '/Users/leung/Documents/docs/2013-2015 NCSU post-doc/2016-02 UCNstorage data/data/';
%[~,~,rawData] = xlsread('./dataProcessing/nEDM_UCNStorage_Feb2016_92K_day1_dataProcess.xls');
[~,~,rawData] = xlsread('./dataProcessing/nEDM_UCNStorage_Feb2016_92K2_day1_dataProcess.xls');
%[~,~,rawData] = xlsread('./dataProcessing/nEDM_UCNStorage_Feb2016_50K_day2_dataProcess.xls');
%[~,~,rawData] = xlsread('./dataProcessing/nEDM_UCNStorage_Feb2016_22K_day3_dataProcess.xls');
%[~,~,rawData] = xlsread('./dataProcessing/nEDM_UCNStorage_Feb2016_15KwExchangeGas_day3_dataProcess');
%[~,~,rawData] = xlsread('./dataProcessing/nEDM_UCNStorage_Feb2016_15KwNoGas_day3_dataProcess');
%[~,~,rawData] = xlsread('./dataProcessing/nEDM_UCNStorage_Feb2016_35K_day4_dataProcess');
%[~,~,rawData] = xlsread('./dataProcessing/nEDM_UCNStorage_Feb2016_58K_day4_dataProcess.xls');

plotDetectorSpectra = 1;
plotDetectorTime = 0;
plotEmptying = 0; % 0 = whole emptying period, 1 = just emptying UCN counts, 2 = just last 100s for BG, 3 = determining drain edge, 4 = start of emptying, 5 = observe drain peak integral region

normDetector = 2; %0=UCNdet, 1=GVmon1, 2=GVmon2, 3=GVmon3, 4=GuideDrainPeakIntegral

%some variables common to all
fillTime = rawData{1,2}; %[s]
GVclosedTime = rawData{2,2}; %[s]
emptyTime = rawData{3,2}; %[s]
UCNdetLLD = rawData{4,2}; UCNdetULD = rawData{5,2}; %channels to make cuts in the UCN detector spectrum
GVmonLLD = rawData{6,2}; GVmonULD = rawData{7,2};
UCNdetPlateauRateNormalization = rawData{8,2}; %[\s] 
GVmon1PlateauRateNormalization = rawData{9,2}; %[\s] 
GVmon2PlateauRateNormalization = rawData{10,2}; %[\s] 
GVmon3PlateauRateNormalization = rawData{11,2}; %[\s]
GuideEmptyPeakIntegralNormalization = rawData{12,2}; %[\s]
lastEmptyTailCorrectionHoldTime = rawData{13,2};
UCNdetAdc = rawData{14,2};
GVmon1Adc = rawData{15,2};
GVmon2Adc = rawData{16,2};
GVmon3Adc = rawData{17,2};
plotYmin = rawData{18,2};
plotYmax = rawData{19,2};
plotXmin = rawData{20,2};
plotXmax = rawData{21,2};

drainEdgeShift = 3; % 3s delay (estimated from looking at data) draining edge and when UCN open
fillingPlateauDelayTime1 = 10+120; %counted backwards
fillingPlateauDelayTime2 = 10;
emptyCountingTime = 400; %[s] time after cell valve open to count emptied UCNs
minHoldTimeForBeforeEmptyBG = 600; %[s] time after cell closed before using data for BG
drainPeakIntegralEnd = 20;

useData = rawData(:,4); % see 'doc cell array' : round brackets to get multiple data from a cell array
useData{22} = NaN;  % set the 'UseData' text string to NaN, considered a number
useData = cell2mat(useData); % convert cell array to normal array
useBG = cell2mat(rawData(useData==1,5));
fileName = rawData(useData==1,3);
fileName = lower(fileName); %make lower case
for i=1:length(fileName)
    runNumber(i) = str2num(strrep(fileName{i}, 'run', ''));
end
holdTime = cell2mat(rawData(useData==1,9));
drainEdgeTime = cell2mat(rawData(useData==1,11));

%initialize some variables
prevHoldTime=0; uniqueHoldTime = 0; filesUsedString=''; numFilesThisHoldTime = 0; j=0; m=1; n=1;
% i is file name index
% j is index for number of files with the same hold time
% n is index if background is used. Need new index because of BG from before cell valve open
% m is when correcting for guide draining tail is needed
% k is for cycling through guide draining tail correction

for i=1:length(fileName)
    currentHoldTime = holdTime(i);
    currentDrainEdgeTime = drainEdgeTime(i);
    cellValveOpenTime = currentDrainEdgeTime+currentHoldTime-drainEdgeShift;
    
    plateauDuration = fillingPlateauDelayTime1 - fillingPlateauDelayTime2;
    plateauStart = currentDrainEdgeTime-fillingPlateauDelayTime1;
    plateauEnd = currentDrainEdgeTime-fillingPlateauDelayTime2;
    
    %Import data
    A = importdata([dataPath fileName{i} '.txt'],' ',1);
    rtc = A.data(:,1);
    UCNdetChan = A.data(:,UCNdetAdc+2);
    GVmon1Chan = A.data(:,GVmon1Adc+2);
    GVmon2Chan = A.data(:,GVmon2Adc+2);
    GVmon3Chan = A.data(:,GVmon3Adc+2);
    GVmon1Chan(GVmon1Chan> 32000) = GVmon1Chan(GVmon1Chan> 32000)-2^15;
    
    [UCNdetCounts,time] = hist(rtc(UCNdetChan>=UCNdetLLD & UCNdetChan<=UCNdetULD),[1:max(rtc(rtc<1e4))]);
    [GVmon1Counts,time] = hist(rtc(GVmon1Chan>GVmonLLD & GVmon1Chan<GVmonULD),[1:max(rtc(rtc<1e4))]);
    [GVmon2Counts,time] = hist(rtc(GVmon2Chan>GVmonLLD & GVmon2Chan<GVmonULD),[1:max(rtc(rtc<1e4))]);
    [GVmon3Counts,time] = hist(rtc(GVmon3Chan>GVmonLLD & GVmon3Chan<GVmonULD),[1:max(rtc(rtc<1e4))]);
    clear A;
    
    [UCNdetSpectrum,channel] = hist(UCNdetChan(UCNdetChan>15),[1:2000]);
    [GVmon1Spectrum,channel] = hist(GVmon1Chan(GVmon1Chan>15),[1:2000]);
    [GVmon2Spectrum,channel] = hist(GVmon2Chan(GVmon2Chan>15),[1:2000]);
    [GVmon3Spectrum,channel] = hist(GVmon3Chan(GVmon3Chan>15),[1:2000]);

    %---Plot detector spectra
    if plotDetectorSpectra==1
        figure(i+1000); clf;
        set(gcf,'units','centimeters');
        pos = get(gcf,'position');
        set(gcf,'position',[5 0 20 30]);
        subplot(2,1,1); hold all; box on; grid on;
        title(fileName{i}, 'Interpreter', 'none')
        xlabel('channel')
        ylabel('counts')
        set(gca,'Xlim',[9 500])
        plot(channel,UCNdetSpectrum,'-');
        plot([UCNdetLLD UCNdetLLD],get(gca,'Ylim'),'k-','linewidth',2)
        plot([UCNdetULD UCNdetULD],get(gca,'Ylim'),'k-','linewidth',2)
        subplot(2,1,2); hold all; box on; grid on;
        xlabel('channel')
        ylabel('counts')
        set(gca,'Xlim',[20 1000])
        set(gca,'Ylim',[1 1E3])
        %set(gca,'Yscale','log')
        h1=plot(channel,GVmon1Spectrum,'-');
        h2=plot(channel,GVmon2Spectrum,'-');
        h3=plot(channel,GVmon3Spectrum,'-');
        legend([h1 h2 h3],'GVmon1', 'GVmon2', 'GVmon3')
        plot([GVmonLLD GVmonLLD],get(gca,'Ylim'),'k-','linewidth',2)
        plot([GVmonULD GVmonULD],get(gca,'Ylim'),'k-','linewidth',2)
    end
    
    %--Plot UCN detector counts vs time raw
    if plotDetectorTime ==1
        figure(i); clf;
        set(gcf,'units','centimeters');
        pos = get(gcf,'position');
        set(gcf,'position',[30 0 20 30]);
        
        subplot(2,1,1); hold all; box on; grid on;
        title(fileName{i}, 'Interpreter', 'none')
        xlabel('time [s]')
        ylabel('count rate [/s]')
        h1=plot(time,UCNdetCounts,'-o');
        h2=plot(time,GVmon1Counts,'-');
        h3=plot(time,GVmon2Counts,'-');
        h4=plot(time,GVmon3Counts,'-');
        legend([h1 h2 h3 h4],'UCN det','GVmon1','GVmon2','GVmon3')
        set(gca,'Xlim',[0 600],'Ylim',[1 1000])
        set(gca,'Yscale','log');
        plot([plateauStart plateauStart],[1E-2 max(UCNdetCounts)],'k--','lineWidth',2,'Color',linecolors(5,:));
        plot([plateauEnd plateauEnd],[1E-2 max(UCNdetCounts)],'k--','lineWidth',2,'Color',linecolors(5,:));
        plot([currentDrainEdgeTime-2 currentDrainEdgeTime-2],[1E-2 max(UCNdetCounts)],'--','lineWidth',2,'Color',linecolors(6,:));
        plot([currentDrainEdgeTime+drainPeakIntegralEnd currentDrainEdgeTime+drainPeakIntegralEnd],[1E-2 max(UCNdetCounts)],'--','lineWidth',2,'Color',linecolors(6,:));
        plot([cellValveOpenTime cellValveOpenTime],[1E-2 max(UCNdetCounts)],'--','lineWidth',2,'Color',linecolors(7,:));
        plot([cellValveOpenTime+emptyCountingTime cellValveOpenTime+emptyCountingTime],[1E-2 max(UCNdetCounts)],'--','lineWidth',2,'Color',linecolors(7,:));
        plot([cellValveOpenTime+emptyTime cellValveOpenTime+emptyTime],[1E-2 max(UCNdetCounts)],'--','lineWidth',2,'Color',linecolors(7,:));
        
        subplot(2,1,2); hold all; box on; grid on;
        title(fileName{i}, 'Interpreter', 'none')
        xlabel('time [s]')
        ylabel('UCN detector counts [/s]')
        plot(time,UCNdetCounts,'-o');
        %plot(time,GVmon1Counts,'-');
        %plot(time,GVmon2Counts,'-');
        %plot(time,GVmon3Counts,'-');
        set(gca,'Ylim',[0.5 5]);
        %set(gca,'Yscale','log');
        plot([plateauStart plateauStart],[1E-2 max(UCNdetCounts)],'k--','lineWidth',2,'Color',linecolors(5,:));
        plot([plateauEnd plateauEnd],[1E-2 max(UCNdetCounts)],'k--','lineWidth',2,'Color',linecolors(5,:));
        plot([currentDrainEdgeTime-2 currentDrainEdgeTime-2],[1E-2 max(UCNdetCounts)],'--','lineWidth',2,'Color',linecolors(6,:));
        plot([currentDrainEdgeTime+drainPeakIntegralEnd currentDrainEdgeTime+drainPeakIntegralEnd],[1E-2 max(UCNdetCounts)],'--','lineWidth',2,'Color',linecolors(6,:));
        plot([cellValveOpenTime cellValveOpenTime],[1E-2 max(UCNdetCounts)],'k--','lineWidth',2,'Color',linecolors(7,:));
        plot([cellValveOpenTime+emptyCountingTime cellValveOpenTime+emptyCountingTime],[1E-2 max(UCNdetCounts)],'k--','lineWidth',2,'Color',linecolors(7,:));
        plot([cellValveOpenTime+emptyTime cellValveOpenTime+emptyTime],[1E-2 max(UCNdetCounts)],'k--','lineWidth',2,'Color',linecolors(7,:));
        
        switch plotEmptying
            case 0
                set(gca,'Xlim',[cellValveOpenTime-100 cellValveOpenTime+emptyTime+100]);
            case 1
                set(gca,'Xlim',[cellValveOpenTime-100 cellValveOpenTime+emptyCountingTime+100]);
            case 2
                set(gca,'Xlim',[cellValveOpenTime+emptyTime-100 cellValveOpenTime+emptyTime+100]);
            case 3
                set(gca,'Xlim',[300 350]);
            case 4
                set(gca,'Xlim',[cellValveOpenTime-20 cellValveOpenTime+20]);
            case 5
                set(gca,'Xlim',[currentDrainEdgeTime-2-20 currentDrainEdgeTime+drainPeakIntegralEnd+20]);
        end
    end
    
    if currentHoldTime ~= prevHoldTime %if next file has different hold time as previous
        if i~=1
            %this is never executed when i = 1 so these values are from below
            %must clear these variables when a new unique hold time is found
            disp(['holding time: ' num2str(prevHoldTime) 's. Files used: ' filesUsedString])
            disp(['emptying peak counts raw:     ' sprintf('%8.0f  ',UCNcounts)])
            disp(['emptying peak counts raw err: ' sprintf('%8.1f  ',sqrt(UCNcounts))])
            disp(['BG rate [mHz]:                ' sprintf('%8.1f  ',UCNBGrate*1000)])
            disp(['BG rate err [mHz]:            ' sprintf('%8.1f  ',UCNBGrate_err*1000)])
            disp(['UCNdetPlateauRate [/s]:       ' sprintf('%8.2f  ',UCNdetPlateauRate)])
            disp(['UCNdetPlateauRate_err [/s]    ' sprintf('%8.2f  ',UCNdetPlateauRate_err)])
            disp(['GVmon1PlateauRate [/s]:       ' sprintf('%8.1f  ',GVmon1PlateauRate)])
            disp(['GVmon1PlateauRate_err [/s]:   ' sprintf('%8.1f  ',GVmon1PlateauRate_err)])
            disp(['GVmon2PlateauRate [/s]:       ' sprintf('%8.1f  ',GVmon2PlateauRate)])
            disp(['GVmon2PlateauRate_err [/s]:   ' sprintf('%8.1f  ',GVmon2PlateauRate_err)])
            disp(['GVmon3PlateauRate [/s]:       ' sprintf('%8.1f  ',GVmon3PlateauRate)])
            disp(['GuideEmptyPeakIntegral:       ' sprintf('%8.1f  ',GuideEmptyPeakIntegral)])
            disp(['GuideEmptyPeakIntegral_err:   ' sprintf('%8.1f  ',GuideEmptyPeakIntegral_err)])
            disp(' ')
            clear UCNcounts UCNBGcounts UCNBGrate UCNBGrate_err UCNdetPlateauCounts UCNdetPlateauRate UCNdetPlateauRate_err GVmon1PlateauCounts GVmon2PlateauCounts GVmon3PlateauCounts GVmon1PlateauRate GVmon2PlateauRate GVmon3PlateauRate GVmon1PlateauRate_err GVmon2PlateauRate_err GVmon3PlateauRate_err GuideEmptyPeakIntegral GuideEmptyPeakIntegral_err
        end
        uniqueHoldTime=uniqueHoldTime+1; %then new hold time
        numFilesThisHoldTime = 1;
        filesUsedString = ['  ' fileName{i}];
        j=1;
    else %same hold time
        numFilesThisHoldTime=numFilesThisHoldTime+1;
        filesUsedString = [filesUsedString '  ' fileName{i}];
        j=j+1;
    end
    prevHoldTime = currentHoldTime;
    
    %--- Calculate count rates and counts
    currentGVmon1PlateauRate = sum(GVmon1Counts(time>=plateauStart & time<=plateauEnd))/plateauDuration;
    currentGVmon2PlateauRate = sum(GVmon2Counts(time>=plateauStart & time<=plateauEnd))/plateauDuration;
    currentGVmon3PlateauRate = sum(GVmon3Counts(time>=plateauStart & time<=plateauEnd))/plateauDuration;
    currentUCNdetPlateauRate = sum(UCNdetCounts(time>=plateauStart & time<=plateauEnd))/plateauDuration;
    currentEmptyPeakIntegral = sum(UCNdetCounts(time>=currentDrainEdgeTime-2 & time<=currentDrainEdgeTime+drainPeakIntegralEnd));
    
    % _all suffix if value stored for all run numbers. For use later.
    UCNdetPlateauRate_all(i) = currentUCNdetPlateauRate;
    GVmon1PlateauRate_all(i) = currentGVmon1PlateauRate;
    GVmon2PlateauRate_all(i) = currentGVmon2PlateauRate;
    GVmon3PlateauRate_all(i) = currentGVmon3PlateauRate;
    EmptyPeakIntegral_all(i) = currentEmptyPeakIntegral;
    
    currentBGcount = sum(UCNdetCounts(time>=cellValveOpenTime+emptyCountingTime & time<=cellValveOpenTime+emptyTime));
    UCNdetBGrateAll(n) = currentBGcount/(emptyTime - emptyCountingTime);
    if currentBGcount == 0
        UCNdetBGrateErrAll(n) = 1/(emptyTime-emptyCountingTime); %if BGcount = 0, set to 1 to get reasonable error (just for plotting)
    else
        UCNdetBGrateErrAll(n) = sqrt(currentBGcount)/(emptyTime-emptyCountingTime);
    end
    runNumberForBG(n) = runNumber(i);
    useBGAll(n) = useBG(i);
    n=n+1;
    
    %Sum UCN counts from start of emptying time for emptyingCountTimes. (cleared above if new holding time)
    UCNcounts(j) = sum(UCNdetCounts(time>=cellValveOpenTime & time<=cellValveOpenTime+emptyCountingTime));
    UCNBGcounts(j) = currentBGcount;
    UCNBGrate(j) = UCNBGcounts(j)/(emptyTime-emptyCountingTime);
    UCNBGrate_err(j) = sqrt(UCNBGcounts(j))/(emptyTime-emptyCountingTime);
    UCNdetPlateauCounts(j)= sum(UCNdetCounts(time>=plateauStart & time<=plateauEnd));
    UCNdetPlateauRate(j) = UCNdetPlateauCounts(j)/plateauDuration;
    UCNdetPlateauRate_err(j) = sqrt(UCNdetPlateauCounts(j))/plateauDuration;
    GVmon1PlateauCounts(j) = sum(GVmon1Counts(time>=plateauStart & time<=plateauEnd));
    GVmon2PlateauCounts(j) = sum(GVmon2Counts(time>=plateauStart & time<=plateauEnd));
    GVmon3PlateauCounts(j) = sum(GVmon3Counts(time>=plateauStart & time<=plateauEnd));
    GVmon1PlateauRate(j) = GVmon1PlateauCounts(j)/plateauDuration;
    GVmon2PlateauRate(j) = GVmon2PlateauCounts(j)/plateauDuration;
    GVmon3PlateauRate(j) = GVmon3PlateauCounts(j)/plateauDuration;
    GVmon1PlateauRate_err(j) = sqrt(GVmon1PlateauCounts(j))/plateauDuration;
    GVmon2PlateauRate_err(j) = sqrt(GVmon2PlateauCounts(j))/plateauDuration;
    GVmon3PlateauRate_err(j) = sqrt(GVmon3PlateauCounts(j))/plateauDuration;
    GuideEmptyPeakIntegral(j) = currentEmptyPeakIntegral;
    GuideEmptyPeakIntegral_err(j) = sqrt(GuideEmptyPeakIntegral(j));

    % Need to normalize after subtracting off BG so can't do it here.
    
    %These are the variables to use later for calculations. 
    %_vec for variables with multiple values per holdTime
    holdTimeUnique(uniqueHoldTime) = currentHoldTime;
    numPasses(uniqueHoldTime) = length(UCNcounts);
    UCNcountsPerEmptyRaw_vec{uniqueHoldTime} = UCNcounts;
    BGrate_vec{uniqueHoldTime} = UCNBGrate;
    BGrate_err_vec{uniqueHoldTime} = UCNBGrate_err;
    UCNdetPlateauRate_vec{uniqueHoldTime} = UCNdetPlateauRate;
    UCNdetPlateauRate_err_vec{uniqueHoldTime}=UCNdetPlateauRate_err;
    GVmon1PlateauRate_vec{uniqueHoldTime} =  GVmon1PlateauRate;
    GVmon1PlateauRate_err_vec{uniqueHoldTime} = GVmon1PlateauRate_err;
    GVmon2PlateauRate_vec{uniqueHoldTime} =  GVmon2PlateauRate;
    GVmon2PlateauRate_err_vec{uniqueHoldTime} = GVmon2PlateauRate_err;
    GVmon3PlateauRate_vec{uniqueHoldTime} =  GVmon3PlateauRate;
    GVmon3PlateauRate_err_vec{uniqueHoldTime} = GVmon3PlateauRate_err;
    GuideEmptyPeakIntegral_vec{uniqueHoldTime} =  GuideEmptyPeakIntegral;
    GuideEmptyPeakIntegral_err_vec{uniqueHoldTime} =  GuideEmptyPeakIntegral_err;
    
    %---calculate the emptying tail correction as a function of holdTime
    if currentHoldTime >= 500  %which data to use
        dummyHoldTime = 10:5:lastEmptyTailCorrectionHoldTime;
        for k = 1:length(dummyHoldTime);
            emptyTailCorrectionUCNdetPlateauRate(m) = UCNdetPlateauRate(j);
            emptyTailCorrectionUCNdetPlateauRate_err(m) = UCNdetPlateauRate_err(j);
            emptyTailCorrectionGVmon1PlateauRate(m) = GVmon1PlateauRate(j);
            emptyTailCorrectionGVmon2PlateauRate(m) = GVmon2PlateauRate(j);
            emptyTailCorrectionGVmon3PlateauRate(m) = GVmon3PlateauRate(j);
            emptyTailCorrectionGVmon1PlateauRate_err(m) = GVmon1PlateauRate_err(j);
            emptyTailCorrectionGVmon2PlateauRate_err(m) = GVmon2PlateauRate_err(j);
            emptyTailCorrectionGVmon3PlateauRate_err(m) = GVmon3PlateauRate_err(j);
            emptyTailCorrectionGuideEmptyPeakIntegral(m) = GuideEmptyPeakIntegral(j);
            emptyTailCorrectionGuideEmptyPeakIntegral_err(m) = GuideEmptyPeakIntegral_err(j);
            
            dummyCellValveOpenTime(k) = currentDrainEdgeTime+dummyHoldTime(k)-5; % 5s delay (estimated from looking at data) draining edge and when UCN open
            emptyTailCorrectionRaw(m,k) = sum(UCNdetCounts(time>=dummyCellValveOpenTime(k) & time<=dummyCellValveOpenTime(k)+emptyCountingTime));
        end
        m=m+1;
    end
    
    %---Use long holding times for BG determination
    if currentHoldTime >= minHoldTimeForBeforeEmptyBG  %which data to use
        currentBGcount = sum(UCNdetCounts(time<=cellValveOpenTime-10 & time>=cellValveOpenTime-10-(currentHoldTime-400)));
        UCNdetBGrateAll(n) = currentBGcount/(currentHoldTime-400-10); %Go back 400s
        if currentBGcount == 0
            UCNdetBGrateErrAll(n) = 1/(currentHoldTime-400-10); %if BGcount = 0, set to 1 to get reasonable error (just for plotting)
        else
            UCNdetBGrateErrAll(n) = sqrt(currentBGcount)/(currentHoldTime-400-10);
        end
        runNumberForBG(n) = runNumber(i)+0.5;
        useBGAll(n) = 1;
        n=n+1;
    end
    
end

%---output to screen the last holding time
disp(['holding time: ' num2str(prevHoldTime) 's. Files used: ' filesUsedString])
disp(['emptying peak counts raw:     ' sprintf('%8.0f  ',UCNcounts)])
disp(['emptying peak counts raw err: ' sprintf('%8.1f  ',sqrt(UCNcounts))])
disp(['BG rate [mHz]:                ' sprintf('%8.1f  ',UCNBGrate*1000)])
disp(['BG rate err [mHz]:            ' sprintf('%8.1f  ',UCNBGrate_err*1000)])
disp(['UCNdetPlateauRate [/s]:       ' sprintf('%8.2f  ',UCNdetPlateauRate)])
disp(['UCNdetPlateauRate_err [/s]    ' sprintf('%8.2f  ',UCNdetPlateauRate_err)])
disp(['GVmon1PlateauRate [/s]:       ' sprintf('%8.1f  ',GVmon1PlateauRate)])
disp(['GVmon1PlateauRate_err [/s]:   ' sprintf('%8.1f  ',GVmon1PlateauRate_err)])
disp(['GVmon2PlateauRate [/s]:       ' sprintf('%8.1f  ',GVmon2PlateauRate)])
disp(['GVmon2PlateauRate_err [/s]:   ' sprintf('%8.1f  ',GVmon2PlateauRate_err)])
disp(['GVmon3PlateauRate [/s]:       ' sprintf('%8.1f  ',GVmon3PlateauRate)])
disp(['GVmon3PlateauRate_err [/s]:   ' sprintf('%8.1f  ',GVmon3PlateauRate_err)])
disp(['GuideEmptyPeakIntegral:       ' sprintf('%8.1f  ',GuideEmptyPeakIntegral)])
disp(['GuideEmptyPeakIntegral_err:   ' sprintf('%8.1f  ',GuideEmptyPeakIntegral_err)])
disp(' ')
clear UCNcounts UCNBGcounts UCNBGrate UCNBGrate_err UCNdetPlateauCounts UCNdetPlateauRate UCNdetPlateauRate_err GVmon1PlateauCounts GVmon2PlateauCounts GVmon3PlateauCounts GVmon1PlateauRate GVmon2PlateauRate GVmon3PlateauRate GVmon1PlateauRate_err GVmon2PlateauRate_err GVmon3PlateauRate_err GuideEmptyPeakIntegral GuideEmptyPeakIntegral_err

%----------------------------------
% Check the monitoring detectors
%---------------------------------
figure(1E5); clf; set(gcf,'units','centimeters'); set(gcf,'position',[15 0 20 20]);
hold all; box on; grid on;
%sort according to runNumber
[~, sortIndex] = sort(runNumber);
disp('----monitoring detectors----')
xlabel('run number')
ylabel('plateau rates normalized')
h1=plot(runNumber(sortIndex),UCNdetPlateauRate_all(sortIndex)/UCNdetPlateauRateNormalization,'-o');
h2=plot(runNumber(sortIndex),GVmon1PlateauRate_all(sortIndex)/GVmon1PlateauRateNormalization,'-o');
h3=plot(runNumber(sortIndex),GVmon2PlateauRate_all(sortIndex)/GVmon2PlateauRateNormalization,'-o');
h4=plot(runNumber(sortIndex),GVmon3PlateauRate_all(sortIndex)/GVmon3PlateauRateNormalization,'-o');
h5=plot(runNumber(sortIndex),EmptyPeakIntegral_all(sortIndex)/GuideEmptyPeakIntegralNormalization,'-o');
legend([h1 h2 h3 h4 h5],'UCNdetPlateauRate','GVmon1PlateauRate','GVmon2PlateauRate','GVmon3PlateauRate','drain peak integral')

figure(1E5+1); clf; set(gcf,'units','centimeters'); set(gcf,'position',[15 0 15 30]);

subplot(3,1,1); hold all; box on; grid on;
xlabel('run number')
ylabel('GVmon2 vs GVmon3')
plot(runNumber(sortIndex),GVmon2PlateauRate_all(sortIndex)./GVmon3PlateauRate_all(sortIndex)/(mean(GVmon2PlateauRate_all)/mean(GVmon3PlateauRate_all)),'-o');

subplot(3,1,2); hold all; box on; grid on;
xlabel('run number')
ylabel('GVmon2 vs drain peak')
plot(runNumber(sortIndex),GVmon2PlateauRate_all(sortIndex)./EmptyPeakIntegral_all(sortIndex)/(mean(GVmon2PlateauRate_all)/mean(EmptyPeakIntegral_all)),'-o');

subplot(3,1,3); hold all; box on; grid on;
xlabel('run number')
ylabel('GVmon3 vs drain peak')
plot(runNumber(sortIndex),GVmon3PlateauRate_all(sortIndex)./EmptyPeakIntegral_all(sortIndex)/(mean(GVmon3PlateauRate_all)/mean(EmptyPeakIntegral_all)),'-o');

disp(['runNumber:            ' sprintf('%8.0f  ',runNumber)])
disp(['UCNdetPlateauRate_all:' sprintf('%8.3f  ',UCNdetPlateauRate_all)])
disp(['GVmon1PlateauRate_all:' sprintf('%8.3f  ',GVmon1PlateauRate_all)])
disp(['GVmon2PlateauRate_all:' sprintf('%8.3f  ',GVmon2PlateauRate_all)])
disp(['GVmon3PlateauRate_all:' sprintf('%8.3f  ',GVmon3PlateauRate_all)])
disp(['EmptyPeakIntegral_all:' sprintf('%8.0f  ',EmptyPeakIntegral_all)])

disp(['UCNdetPlateauRate_mean:' sprintf('%8.3f  ',mean(UCNdetPlateauRate_all))])
disp(['GVmon1PlateauRate_mean:' sprintf('%8.3f  ',mean(GVmon1PlateauRate_all))])
disp(['GVmon2PlateauRate_mean:' sprintf('%8.3f  ',mean(GVmon2PlateauRate_all))])
disp(['GVmon3PlateauRate_mean:' sprintf('%8.3f  ',mean(GVmon3PlateauRate_all))])
disp(['EmptyPeakIntegral_mean:' sprintf('%8.3f  ',mean(EmptyPeakIntegral_all))])

%---- GVmon detector ratios (absolute) to study UCN source spectrum change

GVmon2vsGVmon1 = GVmon2PlateauRate_all(sortIndex)./GVmon1PlateauRate_all(sortIndex);
GVmon3vsGVmon1 = GVmon3PlateauRate_all(sortIndex)./GVmon1PlateauRate_all(sortIndex);
GVmon3vsGVmon2 = GVmon3PlateauRate_all(sortIndex)./GVmon2PlateauRate_all(sortIndex);

figure(1E5+2); clf; set(gcf,'units','centimeters'); set(gcf,'position',[30 0 15 30]);
subplot(3,1,1); hold all; box on; grid on; 
xlabel('run number')
ylabel('GVmon2/GVmon1')
%set(gca,'Ylim',[5.4 7.2])
plot(runNumber(sortIndex),GVmon2vsGVmon1,'-o');
subplot(3,1,2); hold all; box on; grid on;
xlabel('run number')
ylabel('GVmon3/GVmon1')
%set(gca,'Ylim',[2.8 4.0])
plot(runNumber(sortIndex),GVmon3vsGVmon1,'-o');
subplot(3,1,3); hold all; box on; grid on;
%set(gca,'Ylim',[0.51 0.58])
xlabel('run number')
ylabel('GVmon3/GVmon2')
plot(runNumber(sortIndex),GVmon3PlateauRate_all(sortIndex)./GVmon2PlateauRate_all(sortIndex),'-o');

disp(' ')
disp('----monitoring detector plateau rate ratio means----')
disp(['GVmon2/GVmon1:' sprintf('%8.3f %8.3f ',[mean(GVmon2vsGVmon1) std(GVmon2vsGVmon1)/sqrt(length(GVmon2vsGVmon1))])]);
disp(['GVmon3/GVmon1:' sprintf('%8.3f %8.3f ',[mean(GVmon3vsGVmon1) std(GVmon3vsGVmon1)/sqrt(length(GVmon3vsGVmon1))])]);
disp(['GVmon3/GVmon2:' sprintf('%8.3f %8.3f ',[mean(GVmon3vsGVmon2) std(GVmon3vsGVmon2)/sqrt(length(GVmon3vsGVmon2))])]);

%-------------------------------------
% Make a plot of the background rate with time
%-------------------------------------

figure(1E5+3); clf; hold all; box on; grid on;
set(gca,'Ylim',[0 +Inf])
ylabel('Background rate [mHz]')

h1=errorbar(runNumberForBG(useBGAll~=1),UCNdetBGrateAll(useBGAll~=1)*1000, UCNdetBGrateErrAll(useBGAll~=1)*1000, 'o','Color',linecolors(1,:));
h2=errorbar(runNumberForBG(useBGAll==1),UCNdetBGrateAll(useBGAll==1)*1000, UCNdetBGrateErrAll(useBGAll==1)*1000, 'o','Color',linecolors(2,:));
xlabel('run number')
legend([h2 h1],'used in average','not used')

%{
disp('----BG with holding time-----')
%errorbar(holdTime(useBGAll==1),UCNdetBGrateAll(useBGAll==1)*1000, UCNdetBGrateErrAll(useBGAll==1)*1000, 'o','Color',linecolors(2,:));
%xlabel('hold time [s]')
disp(['hold time [s]  BG rate [mHz]   BGrateErr [mHz]' ]);
disp(num2str([holdTime(useBGAll==1) UCNdetBGrateAll(useBGAll==1)'*1000 UCNdetBGrateErrAll(useBGAll==1)'*1000],'%8.0f  %8.2f  %8.2f ; '));
%}

disp(' ')
disp('----averaged background-----')

weightedMean = @(t,p) p(1);
[BGrate, BGrate_err, nchi2] = nonlinft(weightedMean,runNumberForBG(useBGAll==1),UCNdetBGrateAll(useBGAll==1), UCNdetBGrateErrAll(useBGAll==1),[10E-3],[1]);

disp(['BGrate [mHz] =     ' num2str(BGrate*1000)]);
disp(['BGrate_err [mHz] = ' num2str(BGrate_err*1000)]);
disp(['chi2nu =           ' num2str(nchi2)]);

plot(get(gca,'Xlim'),[BGrate-BGrate_err]*1000*ones(1,2),'r')
plot(get(gca,'Xlim'),[BGrate+BGrate_err]*1000*ones(1,2),'r')

%-------------------------------------
% Calculate the effect of the UCN guide emptying tail 
% (uses emptyTailCorrectionRaw from for-loop earlier)
%-------------------------------------

disp(' ')
disp('----Emptying tail counts-----')
disp([ 'dummyHoldTime [s]            : ' sprintf('%7.0f',dummyHoldTime)]);
disp(' ');
emptyTailCorrectionNumPasses = length(emptyTailCorrectionRaw(:,1));
for m=1:emptyTailCorrectionNumPasses
    disp([ 'emptyTailCorrectionRaw       : ' sprintf('%7.0f',emptyTailCorrectionRaw(m,:))]);
    switch normDetector
        case 0
            emptyTailCorrectionNorm(m,:) = emptyTailCorrectionRaw(m,:)*UCNdetPlateauRateNormalization/emptyTailCorrectionUCNdetPlateauRate(m);
        case 1
            emptyTailCorrectionNorm(m,:) = emptyTailCorrectionRaw(m,:)*GVmon1PlateauRateNormalization/emptyTailCorrectionGVmon1PlateauRate(m);
        case 2
            emptyTailCorrectionNorm(m,:) = emptyTailCorrectionRaw(m,:)*GVmon2PlateauRateNormalization/emptyTailCorrectionGVmon2PlateauRate(m);
        case 3
            emptyTailCorrectionNorm(m,:) = emptyTailCorrectionRaw(m,:)*GVmon3PlateauRateNormalization/emptyTailCorrectionGVmon3PlateauRate(m);
        case 4
            emptyTailCorrectionNorm(m,:) = emptyTailCorrectionRaw(m,:)*GuideEmptyPeakIntegralNormalization/emptyTailCorrectionGuideEmptyPeakIntegral(m);
    end
    emptyTailCorrectionNormBGcorr(m,:) = emptyTailCorrectionNorm(m,:) - emptyCountingTime*BGrate;
    disp([ 'emptyTailCorrectionNorm      : ' sprintf('%7.0f',emptyTailCorrectionNorm(m,:))]);
    disp([ 'emptyTailCorrectionNormBGcorr: ' sprintf('%7.0f',emptyTailCorrectionNormBGcorr(m,:))]);
    disp(' ');
end

%--Take the average          
for k = 1:length(dummyHoldTime)
    emptyTailCorrection(k) = sum(emptyTailCorrectionNormBGcorr(:,k))/emptyTailCorrectionNumPasses;
    emptyTailCorrection_err(k) = std(emptyTailCorrectionNormBGcorr(:,k))/sqrt(emptyTailCorrectionNumPasses);
end

disp(' ');
disp('-------------Average Emptying Correction -----------')
disp(['dummyHoldTime [s]            : ' sprintf('%7.0f',dummyHoldTime)]);
disp(['emptyTailCorrection          : ' sprintf('%7.1f',emptyTailCorrection)])
disp(['emptyTailCorrection_err      : ' sprintf('%7.1f',emptyTailCorrection_err)])
disp(' ');

%-------------------------------------
% Correcting the data: BG subtraction, Normalization, and Emptying Tail Correction
%-------------------------------------

for i=1:length(holdTimeUnique)
    %recall: _vec keeps the different passes separate.
    
    %-just for output
    UCNcountsPerEmptyRaw(i) = sum(UCNcountsPerEmptyRaw_vec{i})./numPasses(i);
    UCNcountsPerEmptyRaw_err(i)= sqrt(sum(UCNcountsPerEmptyRaw_vec{i}))./numPasses(i);
    
    %-subtract BG using the weighted mean BG
    countsPerEmptyBGcorrected_vec{i} = UCNcountsPerEmptyRaw_vec{i}-BGrate*emptyCountingTime;
    countsPerEmptyBGcorrected_err_vec{i} = sqrt( UCNcountsPerEmptyRaw_vec{i}+(BGrate_err*emptyCountingTime).^2);

    %weighted mean (for output)
    countsPerEmptyBGcorrected(i) = sum(countsPerEmptyBGcorrected_vec{i}./countsPerEmptyBGcorrected_err_vec{i}.^2)./sum(1./countsPerEmptyBGcorrected_err_vec{i}.^2);
    countsPerEmptyBGcorrected_err(i) = sqrt(1./sum(1./countsPerEmptyBGcorrected_err_vec{i}.^2));
    
    %-normalize
    switch normDetector
        case 0
            countsPerEmptyBGcorrectedNormalized_vec{i} = countsPerEmptyBGcorrected_vec{i}.*UCNdetPlateauRateNormalization./UCNdetPlateauRate_vec{i};
            countsPerEmptyBGcorrectedNormalized_err_vec{i} = sqrt((countsPerEmptyBGcorrected_err_vec{i}./countsPerEmptyBGcorrected_vec{i}).^2 + (UCNdetPlateauRate_err_vec{i}./UCNdetPlateauRate_vec{i}).^2).*countsPerEmptyBGcorrectedNormalized_vec{i};
        case 1
            countsPerEmptyBGcorrectedNormalized_vec{i} = countsPerEmptyBGcorrected_vec{i}.*GVmon1PlateauRateNormalization./GVmon1PlateauRate_vec{i};
            countsPerEmptyBGcorrectedNormalized_err_vec{i} = sqrt((countsPerEmptyBGcorrected_err_vec{i}./countsPerEmptyBGcorrected_vec{i}).^2 + (GVmon1PlateauRate_err_vec{i}./GVmon1PlateauRate_vec{i}).^2).*countsPerEmptyBGcorrectedNormalized_vec{i};
        case 2
            countsPerEmptyBGcorrectedNormalized_vec{i} = countsPerEmptyBGcorrected_vec{i}.*GVmon2PlateauRateNormalization./GVmon2PlateauRate_vec{i};
            countsPerEmptyBGcorrectedNormalized_err_vec{i} = sqrt((countsPerEmptyBGcorrected_err_vec{i}./countsPerEmptyBGcorrected_vec{i}).^2 + (GVmon2PlateauRate_err_vec{i}./GVmon2PlateauRate_vec{i}).^2).*countsPerEmptyBGcorrectedNormalized_vec{i};
        case 3
            countsPerEmptyBGcorrectedNormalized_vec{i} = countsPerEmptyBGcorrected_vec{i}.*GVmon3PlateauRateNormalization./GVmon3PlateauRate_vec{i};
            countsPerEmptyBGcorrectedNormalized_err_vec{i} = sqrt((countsPerEmptyBGcorrected_err_vec{i}./countsPerEmptyBGcorrected_vec{i}).^2 + (GVmon3PlateauRate_err_vec{i}./GVmon3PlateauRate_vec{i}).^2).*countsPerEmptyBGcorrectedNormalized_vec{i};
        case 4
            countsPerEmptyBGcorrectedNormalized_vec{i} = countsPerEmptyBGcorrected_vec{i}.*GuideEmptyPeakIntegralNormalization./GuideEmptyPeakIntegral_vec{i};
            countsPerEmptyBGcorrectedNormalized_err_vec{i} = sqrt((countsPerEmptyBGcorrected_err_vec{i}./countsPerEmptyBGcorrected_vec{i}).^2 + (GuideEmptyPeakIntegral_err_vec{i}./GuideEmptyPeakIntegral_vec{i}).^2).*countsPerEmptyBGcorrectedNormalized_vec{i};
    end
    
    %weighted mean:
    countsPerEmptyBGcorrectedNormalized(i) = sum(countsPerEmptyBGcorrectedNormalized_vec{i}./countsPerEmptyBGcorrectedNormalized_err_vec{i}.^2)./sum(1./countsPerEmptyBGcorrectedNormalized_err_vec{i}.^2);
    countsPerEmptyBGcorrectedNormalized_err(i) = sqrt(1./sum(1./countsPerEmptyBGcorrectedNormalized_err_vec{i}.^2));
    
    %----Correcting for Guide Emptying peak-----
    %if holdTimeUnique(i) < 0 % If not applying emptying peak correction
    if holdTimeUnique(i) <= max(dummyHoldTime)
        if sum(dummyHoldTime==holdTimeUnique(i))==0
            disp('holding time not found in dummyHoldTime!')
            break;
        else
        countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected(i) = countsPerEmptyBGcorrectedNormalized(i)-emptyTailCorrection(dummyHoldTime==holdTimeUnique(i));
        countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected_err(i) = sqrt(countsPerEmptyBGcorrectedNormalized_err(i)^2+emptyTailCorrection_err(dummyHoldTime==holdTimeUnique(i))^2);
        end
    else
        countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected(i) = countsPerEmptyBGcorrectedNormalized(i);
        countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected_err(i) = countsPerEmptyBGcorrectedNormalized_err(i);
    end
end

disp('------Combined data-----')
disp(['UCNPlateauRateNormalization:       ' sprintf('%3.1f ', UCNdetPlateauRateNormalization) ' Hz'] )
disp(['GVmon1PlateauRateNormalization:    ' sprintf('%3.1f ', GVmon1PlateauRateNormalization) ' Hz'] )
disp(['GVmon2PlateauRateNormalization:    ' sprintf('%3.1f ', GVmon2PlateauRateNormalization) ' Hz'] )
disp(['GVmon3PlateauRateNormalization:    ' sprintf('%3.1f ', GVmon3PlateauRateNormalization) ' Hz'] )
disp(['holding time [s]:                  ' sprintf('%6.0f ', holdTimeUnique)])
disp(['number of Passes:                  ' sprintf('%6.0f ', numPasses)])
disp(['UCNcountsPerEmptyRaw:              ' sprintf('%6.1f ', UCNcountsPerEmptyRaw)])
disp(['UCNcountsPerEmptyRaw_err:          ' sprintf('%6.1f ', UCNcountsPerEmptyRaw_err)])
disp(['countsPerEmptyBG:                  ' sprintf('%6.1f ', countsPerEmptyBGcorrected)])
disp(['countsPerEmptyBG_err:              ' sprintf('%6.1f ', countsPerEmptyBGcorrected_err)])
disp(['countsPerEmptyBGNorm:              ' sprintf('%6.1f ', countsPerEmptyBGcorrectedNormalized)])
disp(['countsPerEmptyBGNorm_err:          ' sprintf('%6.1f ', countsPerEmptyBGcorrectedNormalized_err)])
disp(['countsPerEmptyBGNormEmptyTail:     ' sprintf('%6.1f ', countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected)])
disp(['countsPerEmptyBGNormEmptyTail_err: ' sprintf('%6.1f ', countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected_err)])
disp(' ');

figure(1E5+3); clf; hold all; grid on; box on; set(gca,'Yscale','log')
h1=errorbar(holdTimeUnique,countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected,countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected_err,'o');

singleExponential = @(t,p) p(1)*exp(-t/p(2));
[pbest, perror, nchi2] = nonlinft(singleExponential,holdTimeUnique,countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected,countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected_err,[500 300],[1 1]);
h2=plot(linspace(0,1400), singleExponential(linspace(0,1400),pbest),'LineWidth',1.2);
disp('------single exponential--------')
disp(['A =   ' sprintf('%3.0f',pbest(1)) ' +/- ' sprintf('%3.0f',perror(1))])
disp(['tau = ' sprintf('%3.0f',pbest(2)) ' +/- ' sprintf('%3.0f',perror(2)) ' s'])
disp(['norm. chi2= ' sprintf('%3.3f',nchi2)])
disp(' ')

doubleExponential = @(t,p) p(1)*exp(-t/p(2)) + p(3)*exp(-t/p(4));
[pbest, perror, nchi2] = nonlinft(doubleExponential,holdTimeUnique,countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected,countsPerEmptyBGcorrectedNormalizedEmptyTailCorrected_err,[2000 200 200 700],[1 1 1 1]);
h3=plot(linspace(0,1400,5000), doubleExponential(linspace(0,1400,5000),pbest));
xlabel('holding time [s]')
ylabel('UCN per empty')
disp('-------double exponential------------------')
disp(['A_short =   ' sprintf('%3.0f',pbest(1)) ' +/- ' sprintf('%3.0f',perror(1))])
disp(['tau_short = ' sprintf('%3.0f',pbest(2)) ' +/- ' sprintf('%3.0f',perror(2)) ' s'])
disp(['A_long =    ' sprintf('%3.0f',pbest(3)) ' +/- ' sprintf('%3.0f',perror(3))])
disp(['tau_long =  ' sprintf('%3.0f',pbest(4)) ' +/- ' sprintf('%3.0f',perror(4)) ' s'])
disp(['norm. chi2= ' sprintf('%3.3f',nchi2)])

legend([h1 h2 h3],'data','single-exp','two-exp')

%--- Figure stuff----
figure(1E5+3);
set(gca,'Xlim',[plotXmin plotXmax])
set(gca,'Ylim',[plotYmin plotYmax])
