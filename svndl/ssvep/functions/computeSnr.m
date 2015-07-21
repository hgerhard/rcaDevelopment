function [snr] = computeSnr(mainData,noise1Data,noise2Data)

snr = 2*mainData.ampBins./(noise1Data.ampBins+noise2Data.ampBins);