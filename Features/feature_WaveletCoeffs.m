function [WCData, WCNames, WCPars] = feature_WaveletCoeffs(V, ttChannelValidity, ~)
% MClust
% [WCData, WCNames] = feature_WaveletCoeffs(V, ttChannelValidity)
% Calculate wavelet coefficients for each channel
%
% INPUTS
%    V = TT tsd
%    ttChannelValidity = nCh x 1 of booleans
%
% OUTPUTS
%    Data - nSpikes x (NC*nCh) wavelet coeffs
%    Names - "WavCoeff: ch#_coeff#"
%
% ADR April 1998
% BMH June 2016 - modified to return wavelet features
% See standard disclaimer in Contents.m

function [D] = ksdist(S)
    % Get the kolmogorov-smirnov distance between a distribution and
    % a matching normal distribution
    x = linspace(-5, 5, 100); %10 bins per standard deviation, go out to 5 stds
    hS = hist((S-mean(S))/std(S), x); %histogram of z-scored values in S
    cS = cumsum(hS/sum(hS)); %cumulative histogram of S
    cI = cdf('Normal', x, 0, 1); %idealized cdf of a normal dist
    D = max(abs(cS-cI)); %largest difference between 2 dists
end

% Settings
NC = 4; %number of wavelet coeffs to use per TT (uses the NC least normal)
SS = 1; %skip this many scales (1 = skip the highest-freq wavelet scale but use all others)

% Input
TTData = V.data();  %Nspikes-by-Nchannels-by-NsamplesPerSpike array of waveform data
[nSpikes, ~, nSamp] = size(TTData);
NCA = bitshift(nSamp,-SS); %number of ALL coefficients after skipping
f = find(ttChannelValidity); %which TT channels to use

% Allocate output
WCData = zeros(nSpikes, NC*length(f));  %to store the features to return
WCNames = cell(NC*length(f), 1);  %names of each feature
WCPars = {};

% Find per-channel wavelet coeffs
for iCh = 1:length(f)
    
    % Get the wavelet coeffs of each spike on this channel
    tTTData = squeeze(TTData(:,iCh,:))'; %this channel's data (Nsamples-by-Nspikes)
    WC = daub4(tTTData); %discrete wavelet transform each spike
    WC = WC(1:NCA,:); %but only use scales after the ones to skip
    
    % Find which wavelet coeffs are the least normally-distributed
    d = zeros(NCA, 1); %k-s test statistic on each wavelet coeff
    for iWc = 1:NCA %for all coeffs after skipped scales,
        d(iWc) = ksdist(WC(iWc,:)); %get the k-s stat
    end
    [~, lnd] = sort(d, 'descend'); %which WCs are the least normally distributed
        
    % Save the NC least-normally distributed WCs as features
    for iWc = 1:NC %for each coeff to use from this channel,
        WCData(:,(iCh-1)*NC+iWc) = WC(lnd(iWc),:)';  %the iWc-th least normal Wcoeff for each spike
        WCNames{(iCh-1)*NC+iWc} = ['_WavCoeff: ch' num2str(f(iCh)) '_coeff' num2str(iWc)];
    end
    
end

end
