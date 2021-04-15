%  [dBratio, rawValsOut_cl, rawValsOut_clNot] = calcExpectedCLtoCLnotRatio...
%                               (N, cl, clNot, freqs, tokenLen, fs)
%
%
% This script estimates the ratio in power between the CL and CLnot
% components of the stimulus.  This is done by computing the RMS value of
% the summed sinusoids (all with amplitude = 1) across 'N' number of
% iterations. NOTE: each sinusoid is added with *random starting phase*.
%
%  This function is required for the 'osullivanTokenIntensityDeviant_.m
%  stimuli.
%
%  INPUTS:
%   N           Integer number of iterations to be computed for the
%               generation of the RMS estimate.
%   cl          Integer number of frequencies to be used in the computation
%               of one stimulus iteration for the CL tokens.
%   clNot       Integer number of frequencies to be used in the computation
%               of one stimulus iteration for the CLnot tokens.
%   freqs       Array of possible frequencies to be used in the computation
%               of all stimulus interations.
%   tokenLen    Duration in SECONDS of tonal tokens.
%   fs          Sampling frequency of stimuli to be generated.
%
%
%  OUTPUTS: 
%   dBratio             CL-to-CLnot power ratio where a positive value
%                       indicates that the CL component is *greater* in
%                       power than the CLnot component.
%   rawValsOut_cl       RMS value for all 'N' iterations of the stimulus.
%   rawValsOut_clNot    RMS value for all 'N' iterations of the stimulus.
%
%  Created by Darrin K. Reed (1/12/15) from calcExpectedCLamp.m

function [dBratio, rawValsOut_cl, rawValsOut_clNot] = calcExpectedCLtoCLnotRatio(N, cl, clNot, freqs, tokenLen, fs)


t = 0:1/fs:tokenLen;

%--------------------------------------------------------------------------
% Calculation for cl
%--------------------------------------------------------------------------
rawValsOut_cl = zeros(N,1);  % INITIALIZE: output array of raw RMS values for each stimulus iteration

for n = 1:N;
    
    % Randomly select 'cl' number of different frequencies for the
    % iteration 
    tmp = randperm(length(freqs));
    freqsKeep = freqs(tmp(1:cl));

    tmpSig = zeros(1,length(t));  % Initialize array for the signal iteration
    
    % Generate one iteration of the sum-of-sinusoids stimulus
    for i = 1:cl
        randPhase = rand(1)*2*pi;

        t1 = (1 * sin(2*pi*freqsKeep(i).*t + randPhase));
        tmpSig =  t1 + tmpSig;
    end
    
    rawValsOut_cl(n,1) = sqrt(mean(tmpSig.^2));  % Take RMS of sum-of-sinusoids stimulus iteration
    
    
end

rmsEstOut_cl = mean(rawValsOut_cl);  % take mean as the estimated amplitude 


%--------------------------------------------------------------------------
% Calculation for clNot
%--------------------------------------------------------------------------

rawValsOut_clNot = zeros(N,1);  % INITIALIZE: output array of raw RMS values for each stimulus iteration

for n = 1:N;
    
    % Randomly select 'cl' number of different frequencies for the
    % iteration 
    tmp = randperm(length(freqs));
    freqsKeep = freqs(tmp(1:clNot));

    tmpSig = zeros(1,length(t));  % Initialize array for the signal iteration
    
    % Generate one iteration of the sum-of-sinusoids stimulus
    for i = 1:clNot
        randPhase = rand(1)*2*pi;

        t1 = (1 * sin(2*pi*freqsKeep(i).*t + randPhase));
        tmpSig =  t1 + tmpSig;
    end
    
    rawValsOut_clNot(n,1) = sqrt(mean(tmpSig.^2));  % Take RMS of sum-of-sinusoids stimulus iteration
    
    
end

rmsEstOut_clNot = mean(rawValsOut_clNot);  % take mean as the estimated amplitude 


%--------------------------------------------------------------------------
% Ratio calculation
%--------------------------------------------------------------------------

dBratio = 20*log10(rmsEstOut_cl/rmsEstOut_clNot);


