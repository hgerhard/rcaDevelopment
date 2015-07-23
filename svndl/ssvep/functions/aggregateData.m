function [avgData] = aggregateData(rcaData,rcaSettings,keepConditions)
% [avgData] = aggregateBins(rcaData,rcaSettings,[keepConditions])
%
% rcaData: ###
% rcaSettings: ###
%
% avgData is a struct containing subject and trial averaged data that
% contains the following fields:
%
%   realBins: bin-by-harmonic-by-component array of REAL RC coefficients
%   imagBins: bin-by-harmonic-by-component array of IMAGINARY RC coefficients
%   ampBins: bin-by-harmonic-by-component array of RC amplitudes
%   phaseBins: bin-by-harmonic-by-component array of RC phases
%
%
% if keepConditions is true (default = false), condition number adds a 4th
%   dimension to each array.
%%
if nargin<2, error('You must provide the rcaSettings struct created when your rca data were created.'); end
if nargin<4, keepConditions=false; end

nConditions = size(rcaData,1);
nSubjects = size(rcaData,2);
nFreqs = length(rcaSettings.freqsToUse);
nBins = length(rcaSettings.binsToUse);
nCompFromInputData = size(rcaData{1,1},2);

% convert to real/imaginary
[rcaDataReal,rcaDataImag] = getRealImag(rcaData);

if keepConditions        
    realBins = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    imagBins = nan(nBins,nFreqs,nCompFromInputData,nConditions);
    muRcaDataReal = nan(nBins*nFreqs,nCompFromInputData,nConditions);
    muRcaDataImag = nan(nBins*nFreqs,nCompFromInputData,nConditions);
    
    for condNum = 1:nConditions        
        % average over trials and subjects for each condition
        muRcaDataReal(:,:,condNum) = nanmean(cat(3,rcaDataReal{condNum,:}),3);
        muRcaDataImag(:,:,condNum) = nanmean(cat(3,rcaDataImag{condNum,:}),3);
        
        for rc=1:nCompFromInputData
            for f=1:nFreqs
                realBins(:,f,rc,condNum) = muRcaDataReal(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc,condNum);
                imagBins(:,f,rc,condNum) = muRcaDataImag(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc,condNum);
            end
        end        
    end
else    
    % average over trials, subjects AND conditions
    muRcaDataReal = nanmean(cat(3,rcaDataReal{:}),3);
    muRcaDataImag = nanmean(cat(3,rcaDataImag{:}),3);
    
    realBins = zeros(nBins,nFreqs,nCompFromInputData);
    imagBins = zeros(nBins,nFreqs,nCompFromInputData);
    for rc = 1:nCompFromInputData
        for f = 1:nFreqs
            realBins(:,f,rc)=muRcaDataReal(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc);
            imagBins(:,f,rc)=muRcaDataImag(rcaSettings.freqIndices==rcaSettings.freqsToUse(f),rc);
        end
    end    
end

ampBins = sqrt(realBins.^2+imagBins.^2);
phaseBins = atan(imagBins./realBins);

avgData.realBins = realBins;
avgData.imagBins = imagBins;
avgData.ampBins = ampBins;
avgData.phaseBins = phaseBins;
