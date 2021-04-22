%**************************************************************************
% [sigOut,stimTypeArray, frameSamps, deviantIndices, stimParams]...
%               = genSNRdeviantStimulusBlock_training(params, FGtypeFLAG))
%
% This function generates the stimulus used for training in the 
% 'intensityFG_training' part of the expeirment. Depending on the
% "FGtypeFLAG" input arg, either all FOUR stimulus types are used (i.e. (1) figure
% standard, (2) background standard, (3) figure deviant, and (4) background
% deviant) or only TWO (either the two figure or the two background stimuli).
%
% IMPORTANT: In order to preserve the meaning of params.trialNo that other
% parts of the functions might depend on, if only two types of stimuli were
% requested (when "FGtypeFLAG" is 0 or 1), we generate DOUBLE the usual
% number of them, so that params.trialNo remains 4*params.N.
%
%
% INPUTS:
%   params          STRUCT, output of a parameters function, 
%                   e.g. "params_intensityFG". Needs to include fields for
%                   all input args of "osullivanTokenIntensityDeviantFig"
%                   and "osullivanTokenIntensityDeviantback".
%  FGtypeFLAG       Numeric value, one of 0:2. Determines the types of 
%                   stimuli requested for the current training block. 
%                   0 = deviant and non-deviant figure stimuli;
%                   1 = deviant and non-deviant background stimuli;
%                   2 = all types of stimuli;
%
% OUTPUTS:
%   sigOut          3D numeric array holding the generated audio signals. 
%                   Sized [no. of stimuli X samples X channels].  
%   stimTypeArray   Column vector of length 4*N that indicates the order
%                   of the presented test conditions.
%   frameSamps      The duration in SAMPLES of each stimulus. 
%   deviantIndices  2D numeric array, sized [no. of stimuli X no. of
%                   deviant tones]. Holds the indices of deviant notes for
%                   each stimulus.
%   stimParams      Struct, holds the parameters of each stimulus as
%                   as returned by the "osullivanTokenIntensityDeviant..."
%                   functions.
%                   
% USAGE:
%   params = params_intensityFG_training;
%  [sigOut, conditionList, frameSamps, deviantIndices, stimParams]...
%               = genSNRdeviantStimulusBlock_training(params, FGtypeFLAG)
%
%
%  Created: February ??, 2016 by Darrin K. Reed
%  Last Modified:  2021.04 by AB, adapting it to the
%               Intensity-deviant-detection FG task with children at SOTE.
%**************************************************************************

function [sigOut, stimulusTypeArray, frameSamps, deviantIndices, stimParams] = genSNRdeviantStimulusBlock_training(params, FGtypeFLAG)

%% Input checks

if nargin ~= 2
    error('Function "genSNRdeviantStimulusBlock_training" requires input args "params" and "FGtypeFLAG"!');
end


%% Settings, preallocation

plotMe = 0;  % Flag for plotting spectrogram of each individual sound segment... NOT SUGGESTED TO PLOT!

% Generate stimulus type array, depending on FGtypeFLAG, in random order
if FGtypeFLAG == 0  % only figure stimuli
    stimulusTypeArray = Shuffle([params.figStandardLabel*ones(2*params.N, 1); params.figDeviantLabel*ones(2*params.N, 1)]);
elseif FGtypeFLAG == 1  % only background stimuli
    stimulusTypeArray = Shuffle([params.backStandardLabel*ones(2*params.N, 1); params.backDeviantLabel*ones(2*params.N, 1)]);
elseif FGtypeFLAG == 2
    stimulusTypeArray = Shuffle([params.figStandardLabel*ones(params.N, 1); params.backStandardLabel*ones(params.N, 1); params.figDeviantLabel*ones(params.N, 1); params.backDeviantLabel*ones(params.N, 1)]);  
else
    error('Wrong value for "FGtypeFLAG", should be one of 0:2!');
end

% No. of samples per stimulus
frameSamps = round(params.fs*params.lenNote)*params.numNotesPerToken;  % No. of samples per stimulus

% preallocate output arrays
sigOut = nan(length(stimulusTypeArray), frameSamps, 2);  % preallocate variable holding all stimuli, sized [no. of stimuli X samples X channels]
deviantIndices = nan(length(stimulusTypeArray), params.deviantDuration);  % preallocate variable holding the indices of deviant notes


%% Loop over stimuli, generating them one-by-one

for i = 1:length(stimulusTypeArray)

    % Create the (random) indices for the deviant
    tmp = randi([params.deviantIndMin, params.deviantIndMax],1);
    deviantInd = tmp:(tmp+params.deviantDuration-1);

    switch stimulusTypeArray(i)
        case params.figStandardLabel  % Generate NON-DEVIANT FIGURE
            [sigOut(i, :, :), trialParams] = osullivanTokenIntensityDeviantFig(params.fs, params.clFig, params.clNotFig,...
                round(params.lenNote*params.fs), round(params.lenRamp*params.fs), params.numNotesPerToken,...
                params.minSemitoneStep, params.clITD, params.clNotITD, params.clTraj, [params.baseFreq params.upFreq],...
                params.maxSemitoneStep, params.clStart, [1 1], plotMe,...
                params.baseSNR, deviantInd, params.NONdeviantSNR);   

        case params.backStandardLabel  % Generate NON-DEVIANT BACKGROUND
            [sigOut(i, :, :), trialParams] = osullivanTokenIntensityDeviantBack(params.fs, params.clBack, params.clNotBack, ...
                round(params.lenNote*params.fs), round(params.lenRamp*params.fs), params.numNotesPerToken,...
                params.minSemitoneStep, params.clITD, params.clNotITD, params.clTraj, [params.baseFreq params.upFreq],...
                params.maxSemitoneStep, params.clStart, [1 1], plotMe,...
                params.deviantNum, deviantInd, params.NONdeviantSNR);

        case params.figDeviantLabel  % Generate DEVIANT FIGURE
            [sigOut(i, :, :), trialParams] = osullivanTokenIntensityDeviantFig(params.fs, params.clFig, params.clNotFig,...
                round(params.lenNote*params.fs), round(params.lenRamp*params.fs), params.numNotesPerToken,...
                params.minSemitoneStep, params.clITD, params.clNotITD, params.clTraj, [params.baseFreq params.upFreq],...
                params.maxSemitoneStep, params.clStart, [1 1], plotMe,...
                params.baseSNR, deviantInd, params.deviantSNR);   

        case params.backDeviantLabel  % Generate DEVIANT BACKGROUND
            [sigOut(i, :, :), trialParams] = osullivanTokenIntensityDeviantBack(params.fs, params.clBack, params.clNotBack, ...
                round(params.lenNote*params.fs),  round(params.lenRamp*params.fs),  params.numNotesPerToken,...
                params.minSemitoneStep, params.clITD, params.clNotITD, params.clTraj, [params.baseFreq params.upFreq],...
                params.maxSemitoneStep, params.clStart, [1 1], plotMe,...
                params.deviantNum, deviantInd, params.deviantSNR); 

        otherwise
            error('DKR: Invalid stimulus type.')
    end
    
    % collect indices of deviant notes
    deviantIndices(i, :) = deviantInd;
    % collect parameters for each stimulus
    if i == 1
        stimParams = trialParams;
    else
        stimParams(i) = trialParams;
    end

end



return
