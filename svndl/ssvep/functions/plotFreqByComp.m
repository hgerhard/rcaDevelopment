function [figNums] = plotFreqByComp(mainRcaData,mainNoiseData,rcaSettings,plotSettings,comparisonRcaData,comparisonNoiseData)

% ### add functionality for condition separation
if nargin<2, error('You must provide both signal and noise data.'); end
if nargin<3, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if nargin<4, useSpecialSettings = false; end
if (nargin>=4 && ~isempty(plotSettings)), useSpecialSettings = true; else useSpecialSettings = false; end
if nargin<6, plotComparison = false; else plotComparison = true; end

poolOverBins = false;

figNums = [];

nFreqs = length(rcaSettings.freqsToUse);

avgRcaData = aggregateData(mainRcaData,rcaSettings);
avgNoise1Data = aggregateData(mainNoiseData.lowerSideBand,rcaSettings);
avgNoise2Data = aggregateData(mainNoiseData.higherSideBand,rcaSettings);

[~,noiseLevsMain] = computeSnr(avgRcaData,avgNoise1Data,avgNoise2Data,poolOverBins);

if plotComparison
    avgCompData = aggregateData(comparisonRcaData,rcaSettings);
    avgCompNoise1Data = aggregateData(comparisonNoiseData.lowerSideBand,rcaSettings);
    avgCompNoise2Data = aggregateData(comparisonNoiseData.higherSideBand,rcaSettings);
    
    [~,noiseLevsCompare] = computeSnr(avgCompData,avgCompNoise1Data,avgCompNoise2Data,poolOverBins);
end

figure
for f=1:nFreqs
    for rc=1:rcaSettings.nComp
        
        if plotComparison
            subplot(nFreqs,rcaSettings.nComp+1,(f-1)*(rcaSettings.nComp+1)+rc);
        else
            subplot(nFreqs,rcaSettings.nComp,(f-1)*(rcaSettings.nComp)+rc);
        end
        
        plot(avgRcaData.ampBins(:,f,rc),'k-','LineWidth',1.5);
        hold on
        plot(noiseLevsMain(:,f,rc),'ks');
        set(gca,'XTickLabel',round(rcaSettings.binLevels([1 5 end])*100)./100); % assumes Matlab puts 3 ticks on x-axis ###
        
        if f==1, title(['RC' num2str(rc)']); end
        
        if rc==1, ylabel(rcaSettings.freqLabels{f}); end
        
        if plotComparison
            subplot(nFreqs,rcaSettings.nComp+1,f*(rcaSettings.nComp+1));
            
            plot(avgCompData.ampBins(:,f,1),'r-','LineWidth',1.5);
            hold on
            plot(noiseLevsCompare(:,f,1),'rs');
            set(gca,'XTickLabel',round(rcaSettings.binLevels([1 5 end])*100)./100); % assumes Matlab puts 3 ticks on x-axis ###
            if f==1
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





