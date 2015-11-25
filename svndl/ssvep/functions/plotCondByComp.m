function [figNums] = plotCondByComp(mainRcaData,mainNoiseData,rcaSettings,plotSettings,comparisonRcaData,comparisonNoiseData)
% [figNums] = plotFreqByComp(mainRcaData,mainNoiseData,rcaSettings,[plotSettings],[comparisonRcaData],[comparisonNoiseData])
%
% Create a numberConditions x numberRCs multipanel figure where the
% amplitude in microVolts is plotted against the bin levels. Optionally
% include an extra column for the comparison channel's data (e.g. for Oz).
% 
% Individual plots are modelled after PowerDiva with unfilled squares to
% indicate noise estimates for each bin.
%


if nargin<2, error('You must provide both signal and noise data.'); end
if nargin<3, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if (nargin<4 || isempty(plotSettings)), useSpecialSettings = false; else useSpecialSettings = true; end
if nargin<6, plotComparison = false; else plotComparison = true; end

poolOverBins = false;
nCond = length(rcaSettings.condsToUse);
separateConds = true;
allCondsInd = 1:nCond;
errorType = plotSettings.errorType;

figNums = [];

nFreqs = length(rcaSettings.freqsToUse);
nBins = length(rcaSettings.binsToUse);
nBinsHalf = floor(nBins/2);
binsToTick = [1 nBinsHalf nBins];

avgRcaData = aggregateData(mainRcaData,rcaSettings,separateConds,errorType);
avgNoise1Data = aggregateData(mainNoiseData.lowerSideBand,rcaSettings,separateConds,errorType);
avgNoise2Data = aggregateData(mainNoiseData.higherSideBand,rcaSettings,separateConds,errorType);

[~,noiseLevsMain] = computeSnr(avgRcaData,avgNoise1Data,avgNoise2Data,poolOverBins);

if plotComparison
    avgCompData = aggregateData(comparisonRcaData,rcaSettings,separateConds,errorType);
    avgCompNoise1Data = aggregateData(comparisonNoiseData.lowerSideBand,rcaSettings,separateConds,errorType);
    avgCompNoise2Data = aggregateData(comparisonNoiseData.higherSideBand,rcaSettings,separateConds,errorType);
    
    [~,noiseLevsCompare] = computeSnr(avgCompData,avgCompNoise1Data,avgCompNoise2Data,poolOverBins);
end

figure
for condNum = allCondsInd  
    for rc=1:rcaSettings.nComp
        
        if plotComparison
            subplot(nCond,rcaSettings.nComp+1,(condNum-1)*(rcaSettings.nComp+1)+rc);
        else
            subplot(nCond,rcaSettings.nComp,(condNum-1)*(rcaSettings.nComp)+rc);
        end
        
        for f=1:nFreqs    
            sc = 1-f/nFreqs;
            if ~mod(f,2)
                plot(avgRcaData.ampBins(:,f,rc,condNum),'k--','LineWidth',1.5,'Color',sc.*plotSettings.conditionColors(condNum,:));
            else
                plot(avgRcaData.ampBins(:,f,rc,condNum),'k-','LineWidth',1.5,'Color',sc.*plotSettings.conditionColors(condNum,:));
            end
            hold on
        end
        if (rc==1), legend(rcaSettings.freqLabels,'Location','NorthWest'); end
        for f=1:nFreqs
            sc = 1-f/nFreqs;
            %plot(noiseLevsMain(:,f,rc,condNum),'ks','Color',sc.*plotSettings.conditionColors(condNum,:));
            if ~isempty(errorType)                
                lb = avgRcaData.ampBins(:,f,rc,condNum) - avgRcaData.ampErrBins(:,f,rc,condNum,1);
                ub = avgRcaData.ampErrBins(:,f,rc,condNum,2) - avgRcaData.ampBins(:,f,rc,condNum);
                for b = 1:nBins
%                     if ~mod(f,2)
%                         errorbar(b,avgRcaData.ampBins(b,f,rc,condNum),lb(b),ub(b),...
%                             'Color',sc.*plotSettings.conditionColors(condNum,:),'LineWidth',1);
%                     else
%                         errorbar(b,avgRcaData.ampBins(b,f,rc,condNum),lb(b),ub(b),...
%                             'Color',sc.*plotSettings.conditionColors(condNum,:),'LineWidth',1.75);
%                     end                        
                end
            end
        end
        set(gca,'XTick',binsToTick);
        set(gca,'XTickLabel',round(rcaSettings.binLevels{condNum}(binsToTick)*100)./100); % ### assumes all conditions plotted have same binLevels
        xlim([0 nBins+1]);
        
        if rc==1, ylabel(['Condition ',num2str(condNum),' Ampl.']); end
        if condNum==1, title(['RC' num2str(rc)']); end 
        
        if plotComparison
            subplot(nCond,rcaSettings.nComp+1,condNum*(rcaSettings.nComp+1));
            
            for f=1:nFreqs
                sc = 1-f/nFreqs;
                if ~mod(f,2)
                    plot(avgCompData.ampBins(:,f,1,condNum),'r--','LineWidth',1.5,'Color',sc.*plotSettings.conditionColors(condNum,:));
                else
                    plot(avgCompData.ampBins(:,f,1,condNum),'r-','LineWidth',1.5,'Color',sc.*plotSettings.conditionColors(condNum,:));
                end
                hold on
                %plot(noiseLevsCompare(:,f,1,condNum),'rs','Color',sc.*plotSettings.conditionColors(condNum,:));
                if ~isempty(errorType)
                    lb = avgCompData.ampBins(:,f,1,condNum) - avgCompData.ampErrBins(:,f,1,condNum,1);
                    ub = avgCompData.ampErrBins(:,f,1,condNum,2) - avgCompData.ampBins(:,f,1,condNum);
                    for b = 1:nBins
%                         if ~mod(f,2)
%                             errorbar(b,avgCompData.ampBins(b,f,1,condNum),lb(b),ub(b),...
%                                 'Color',sc.*plotSettings.conditionColors(condNum,:),'LineWidth',1);
%                         else
%                             errorbar(b,avgCompData.ampBins(b,f,1,condNum),lb(b),ub(b),...
%                                 'Color',sc.*plotSettings.conditionColors(condNum,:),'LineWidth',1.75);                            
%                         end
                    end
                end
            end
            set(gca,'XTick',binsToTick);
            set(gca,'XTickLabel',round(rcaSettings.binLevels{condNum}(binsToTick)*100)./100); % ### assumes all conditions plotted have same binLevels
            xlim([0 nBins+1]);
            xlim([0 nBins+1]);
            if condNum==1
                if useSpecialSettings
                    title(plotSettings.comparisonName)
                else
                    title('Comparison Chan.')
                end
            end
        end
        
    end
end
figNums = [figNums,gcf];





