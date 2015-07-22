function plotSnr(mainRcaData,mainNoiseData,rcaSettings,plotSettings,comparisonRcaData,comparisonNoiseData)

% ### add functionality for condition separation

if nargin<2, error('You must provide both signal and noise data.'); end
if nargin<3, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if nargin<4, useSpecialSettings = false; end
if (nargin>=4 && ~isempty(plotSettings)), useSpecialSettings = true; else useSpecialSettings = false; end
if nargin<6, plotComparison = false; else plotComparison = true; end

poolOverBins = true;

nFreqs = length(rcaSettings.freqsToUse);

avgRcaData = aggregateData(mainRcaData,rcaSettings);
avgNoise1Data = aggregateData(mainNoiseData.lowerSideBand,rcaSettings);
avgNoise2Data = aggregateData(mainNoiseData.higherSideBand,rcaSettings);

snrMain = computeSnr(avgRcaData,avgNoise1Data,avgNoise2Data,poolOverBins);

if plotComparison
    avgCompData = aggregateData(comparisonRcaData,rcaSettings);
    avgCompNoise1Data = aggregateData(comparisonNoiseData.lowerSideBand,rcaSettings);
    avgCompNoise2Data = aggregateData(comparisonNoiseData.higherSideBand,rcaSettings);
    
    snrComparison = computeSnr(avgCompData,avgCompNoise1Data,avgCompNoise2Data,poolOverBins);
end

for rc=1:rcaSettings.nComp
    figure;
    set(gca,'Color','w');
    for f=1:nFreqs
        subplot(nFreqs,1,f); hold on
        plot(snrMain(:,f,rc),'-ok','MarkerFaceColor','k');
        dataLabels = sprintf('RC%d',rc);
        if plotComparison
            plot(snrComparison(:,f,1),'-or','MarkerFaceColor','r'); 
            if useSpecialSettings
                dataLabels = {dataLabels,plotSettings.comparisonName};
            else
                dataLabels = {dataLabels,'Comparison'};
            end
        end
        title(rcaSettings.freqLabels{f});

        ylabel('SNR');
        set(gca,'XTickLabel',round(rcaSettings.binLevels*100)./100);
        if f==1, hlg=legend(dataLabels,'Location','NorthWest'); end
    end
    set(hlg,'box','off');
end