function intensityFG_main(subjID, deviantType, varargin)
%% Intensity-deviant-detection Figure-ground experiment, main part
%
% USAGE: intensityFG_main(subjID, deviantType, blockNumber=1, serialTrigger='notrigger') 
% 
% Main experimental function for the Intensity-deviant-detection
% Figure-ground experiment. 
% The task is to detect an intensity deviant either in the "figure" or the
% "background" portion of the SFG stimuli. There are four types of trials:
%   - figure and no deviant; 
%   - figure and deviant figure; 
%   - no figure and no deviant; 
%   - no figure and deviant background.
% These trials are mixed randomly and are organized into blocks.
%
% Main parameters are defined in "params_intensityFG".
%
% Mandatory inputs:
% subjID            - Numeric value, one of 1:99. Subject ID.
% deviantType       - Char array, one of {'figure', 'background'}. Determines
%                   the task: deviant detection in the figure or the
%                   background.
%
% Optional inputs:
% blockNumber       - Numeric value, one of 1:99. Number (ID) of block to
%                   start from. Defaults to 1 (=start from first block).
%                   The overall number of blocks for the study is determined by
%                   determined by "params.blockNo" in the params struct
%                   (see "params_intensityFG"). This input arg is for
%                   restarting the experiment at a certain block number.
% serialTrigger     - Char array, one of {'notrigger', 'trigger'}. Controls
%                   if the function sends triggers to the serial port 
%                   for EEG recordings or not. Serial port usage is 
%                   determined by "params.serial" and "params.trig", 
%                   from the output of "params_intensityFG.m". 
%                   Defaults to 'notrigger'.
%
% Outputs:
% All outputs are saved out into a subject- and run-specific .mat file. See
% "saveFcn_intensityFG" for details of filename generation.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Based mostly on previous code written by:
% - Darrin K. Reed (currently a Senior Research Audiologist at Apple)
% Credit also to:
% - Brigitta Tóth, at RCNS, Budapest
% - Tamas Kurics, Zsuzsanna Kocsis, Botond Hajdu (all ex-RCNS)
%
% Regarding the experiment, please contact the PI:
% Orsolya Szalardy, szalardy.orsolya@med.semmelweis-univ.hu
%
% Written / edited from pieces by Adam Boncz
% 2021.04.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Input checks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check number of inputs
if ~ismember(nargin, 2:4)    
    error(['Wrong number of input args! Input args "subjID" ',...
        'and "deviantType" are mandatory while input args "blockNumber" ',...
        'and "serialTrigger" are optional!']);
end

% check mandatory inputs
if ~isnumeric(subjID) || ~ismember(subjID, 1:99)
    error('Input arg "subjID" should be an integer in range 1:99!');
end
if ~ischar(deviantType) || ~ismember(deviantType, {'figure', 'background'})
    error('Input arg "deviantType" should be either "figure" or "background"!');
end

% check optional inputs
if ~isempty(varargin)
    for v = 1:numel(varargin)
        if isnumeric(varargin{v}) && ismember(varargin{v}, 1:99) && ~exist('blockNumber', 'var')
            blockNumber = varargin{v};
        elseif ischar(varargin{v}) && ismember(varargin{v}, {'trigger', 'notrigger'}) && ~exist('serialTrigger', 'var')
            serialTrigger = varargin{v};
        else
            error('At least one input arg could not be mapped nicely to "blockNumber" or "serialTrigger"!');
        end
    end  % for v
end  % if ~isempty

% assign default values
if ~exist('blockNumber', 'var')
    blockNumber = 1;
end
if ~exist('serialTrigger', 'var')
    serialTrigger = 'notrigger';
end

% get logical flag and text for user message from "serialTrigger"
if strcmp(serialTrigger, 'trigger')
    triggerFlag = true;
    triggerText = 'Yes';
elseif strcmp(serialTrigger, 'notrigger')
    triggerFlag = false;
    triggerText = 'No';
end

% user message 
disp([char(10), 'Started the Intensity-deviant-detection Figure-ground experiment with inputs:',...
    char(10), 'Subject ID: ', num2str(subjID),...
    char(10), 'Deviant part of stimuli: ', deviantType,...
    char(10), 'Starting with block number ', num2str(blockNumber),...
    char(10), 'Triggers on serial port: ', triggerText]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Basic settings, parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for psychtoolbox installation + unify key names
PsychDefaultSetup(1);
Priority(1);

% load parameters
params = params_intensityFG;

% init serial port if triggering was requested
if triggerFlag
    serialObj = serial(params.serial.portName);
end

% Set response key code (response key is set in "params" struct
detectKeyCode = KbName(params.detectKey);
% Same for block / task start key and the abort key
goKeyCode = KbName(params.goKey);
abortKeyCode = KbName(params.abortKey);
% if there are multiple keys with the name (e.g. 'Return'), only use the
% first one
if numel(detectKeyCode) ~= 1; detectKeyCode = detectKeyCode(1); end
if numel(goKeyCode) ~= 1; goKeyCode = goKeyCode(1); end
if numel(abortKeyCode) ~= 1; abortKeyCode = abortKeyCode(1); end

% Define a flag for the auditory deviant type
if strcmpi(deviantType, 'figure')
    auditoryFigDeviantFLAG = 1;
else
    auditoryFigDeviantFLAG = 0;
end

% restrict keys for KbCheck
RestrictKeysForKbCheck([detectKeyCode, goKeyCode, abortKeyCode]);

% Force costly mex functions into memory to avoid latency later on
GetSecs; WaitSecs(0.1); KbCheck();

% Store the name of the current experimental function, used later when
% results are saved out
expName = mfilename(pwd);

% user message
disp([char(10), 'Basic parameters set / loaded', char(10)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Init PsychPortAudio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load the sound driver - hopefully WASAPI if on Win 10
InitializePsychSound;

% % This part is for selecting a very specific sound device if the default
% % does not work out-of-the-box
% % get correct audio device
device = [];  % system default is our default as well
% targetString = 'ESI Juli@: ICE1724 (hw:2,0)';  % target device name
% % we only change audio device in the lab, when we see the correct audio
% % card
% tmpDevices = PsychPortAudio('GetDevices');
% for i = 1:numel(tmpDevices)
%     if strcmp(tmpDevices(i).DeviceName, targetString)
%         device = tmpDevices(i).DeviceIndex;
%     end
% end

% mode is simple playback
mode = 1;
% reqlatencyclass is set to low-latency
reqLatencyClass = 2;
% 2 channels output
nrChannels = 2;

% open PsychPortAudio device for playback
pahandle = PsychPortAudio('Open', device, mode, reqLatencyClass, params.fs, nrChannels);

% get and display device status
pahandleStatus = PsychPortAudio('GetStatus', pahandle);
disp([char(10), 'PsychPortAudio device status: ']);
disp(pahandleStatus);

% initial start & stop of audio device to avoid potential initial latencies
tmpSound = zeros(2, params.fs/10);  % silence
tmpBuffer = PsychPortAudio('CreateBuffer', [], tmpSound);  % create buffer
PsychPortAudio('FillBuffer', pahandle, tmpBuffer);  % fill the buffer of audio device with silence
PsychPortAudio('Start', pahandle, 1);  % start immediately
PsychPortAudio('Stop', pahandle, 1);  % stop when playback is over

% user message
disp([char(10), 'PsychPortAudio initialized', char(10)]);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Block loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate the stimulus, prepare for trial start

% user message
disp([char(10), 'Done with general settings, preparing for first block...']);

% set flag for aborting experiment
abortFlag = 0;

for blockIdx = blockNumber:params.blockNo
  
    % user message
    disp([char(10), 'Preparing block no. ', num2str(blockIdx), ',',...
        char(10), 'generating stimuli and preallocating / appending response variables...']);
    
    % Generate the stimulus
    [sigOut, trialTypes, ~, deviantIndices, stimParams] = genSNRdeviantStimulusBlock_perTrial(params);
    
    % Scale the output signal.  First normalize to RMS == 1 and then scale
    % by 'params.dbscl'.
    sigOut = sigOut ./ sqrt(mean(sigOut.^2, 2));  % broadcasting
    sigOut = sigOut .* (10^(params.dbscl/20)); 
    
    % Store the deviant note indices and the general stimulus params
    if blockIdx == 1
        allDeviantIndices = nan(params.trialNo, params.deviantDuration, params.blockNo);
        allStimParams = repmat(stimParams, [params.blockNo, 1]);
    else
        allStimParams(blockIdx, :) = stimParams;
    end
    allDeviantIndices(:, :, blockIdx) = deviantIndices;
    
    % Preallocate response variables
    if blockIdx == 1
        outData.rt = nan(params.trialNo, params.blockNo);  % Numeric, response time relative to stimulus onset
        outData.rtDeviant = nan(params.trialNo, params.blockNo);  % Numeric, response time relative to deviant onset
        outData.deviantOnset = nan(params.trialNo, params.blockNo);  % Numeric, deviant onset time relative to stimulus onset
        outData.accuracy = nan(params.trialNo, params.blockNo);  % Binary [1 = correct, 0 = wrong]
        outData.targetPresented = nan(params.trialNo, params.blockNo);  % Binary [1 = target, 0 = no target] 
        outData.trialTypes = nan(params.trialNo, params.blockNo);  % stimulus / trial type, numeric code
    end
    
    % Collect trial / stimulus type information for block
    outData.trialTypes(:, blockIdx) = trialTypes;
    
    % Collect deviant onsets (in time, relative to audio onset)
    outData.deviantOnset(:, blockIdx) = (deviantIndices(:, 1)-1) * params.lenNote + 0.001;
    
    % Send AUDITORY targetMask (i.e. indicies of deviant/target trials)
    % for the i-th block to the output arrays 
    if auditoryFigDeviantFLAG
        deviantMask = trialTypes == 30;     % ------- NOTE  ----------->              This statement is very specific to the conditionID values specified in 'getnSNRdeviantStimulusBlock_perTrial.m'!!!!!!!!!!!!!!!!!
        targetDeviantID = 30;  % use for quick calculation of correctness
        nontargetDeviantID = 40;  % use for quick calculation of correctness
    else
        deviantMask = trialTypes == 40;     % ------- NOTE  ----------->              This statement is very specific to the conditionID values specified in 'getnSNRdeviantStimulusBlock_perTrial.m'!!!!!!!!!!!!!!!!!
        targetDeviantID = 40;  % use for quick calculation of correctness
        nontargetDeviantID = 30;  % use for quick calculation of correctness
    end
    outData.targetPresented(:, blockIdx) = deviantMask;  % Binary [1 = target, 0 = no target]
         
    % user message
    disp('Done, now buffering all audio in advance...');    
    
    % % fill buffers with audio for whole block
    buffer = [];
    for tmpIdx = 1:params.trialNo
        audioData = squeeze(sigOut(tmpIdx, :, :));
        buffer(end+1) = PsychPortAudio('CreateBuffer', [], audioData');
    end
    
    % user message
    disp('Done.');    
    
    % user message
    disp([char(10), 'Ready to start the block!']);
    
    % Define again the task and the response key for user, just to be safe
    disp([char(10), char(10)]);
    disp('-------------------------------------------------------------');
    disp(['THE TASK IS TO DETECT DEVIANTS IN THE   ', upper(deviantType), '']);
    disp('-------------------------------------------------------------');
    disp(['RESPONSE KEY FOR DEVIANT DETECTION:   ', upper(params.detectKey), '']);
    disp('-------------------------------------------------------------');
    disp([char(10), char(10)]);    
    
    
    % start block after user pressed goKey
    if strcmp(params.goKey, 'Return')
        startKey = 'ENTER'; 
    else
        startKey = upper(params.goKey); 
    end
    disp([char(10), char(10)]);
    disp('-------------------------------------------------------------');
    disp(['-----------    PRESS ', startKey, ' TO START THE BLOCK!    -----------']);
    disp('-------------------------------------------------------------');
    disp([char(10), char(10)]);
    
    % wait for key press
    KbReleaseWait;
    while 1
        [keyDown, secs, keyCode] = KbCheck;
        if keyDown && find(keyCode) == goKeyCode
            blockStartTime = secs + 3;
            break;
        elseif keyDown && find(keyCode) == abortKeyCode
            abortFlag = 1;
            break;
        end          
    end
    
    % check for abort
    if abortFlag
        disp([char(10), char(10), 'USER REQUESTED ABORT!', char(10), char(10)]);
        Priority(0);
        RestrictKeysForKbCheck([]);
        PsychPortAudio('Stop', pahandle, 1, 1);
        PsychPortAudio('Close');
        return;
    end   
    
    % block start trigger
    if triggerFlag
        fprintf(serialObj, params.trig.format, params.trig.blockStart);  % block start trigger
        fprintf(serialObj, params.trig.format, params.trig.blockStart + blockIdx);  % block number trigger = block start trigger + block number
    end    
    
    % user message
    disp('Starting...');
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Trial loop
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for trialIdx = 1:params.trialNo
        
        % make sure the subject is not pressing the detection key already
        KbReleaseWait;
        
        % trial start trigger
        if triggerFlag
            fprintf(serialObj, params.trig.format, params.trig.trialStart);  % trial start trigger
            fprintf(serialObj, params.trig.format, params.trig.trialStart + trialIdx);  % trial number trigger = trial start trigger + trial number
            fprintf(serialObj, params.trig.format, trialTypes(trialIdx));  % trial type trigger, from var "trialTypes" 
        end
        
        % trial start time is determined by block start for first trial,
        % otherwise it is last trial + audio length + iti
        if trialIdx == 1
            trialStartTime = blockStartTime;
        else
            trialStartTime = oldTrialStartTime + params.stimLength + params.iti;
        end
        
        % fill audio buffer with next stimuli
        PsychPortAudio('FillBuffer', pahandle, buffer(trialIdx));
        
        % user message
        disp([char(10), 'Trial ', num2str(trialIdx), ', block ', num2str(blockIdx)]);
        
        % wait till we are 50 ms from the start of the playback
        while trialStartTime - GetSecs >= 0.05
            WaitSecs(0.01);
        end
        
        % blocking playback start for precision
        audioOnset = PsychPortAudio('Start', pahandle, 1, trialStartTime, 1);        
        
        % stimulus onset trigger
        if triggerFlag
            fprintf(serialObj, params.trig.format, params.trig.soundOnset);
        end        
        
        % user message
        disp(['Audio started at ', num2str(audioOnset-trialStartTime), ' secs relative to requested start time']);
        if trialIdx ~= 1
            disp(['Audio start relative to last audio onset was ', num2str(audioOnset-oldTrialStartTime), ' secs ',...
                '(should be ', num2str(params.stimLength + params.iti), ' secs)']);
        end
        if deviantMask(trialIdx)
            disp('There is a TARGET in the trial');
            disp(['Deviant onset is at ', num2str(outData.deviantOnset(trialIdx, blockIdx)), ' secs from audio onset']);
        end
    
        % Listen to response until we are only 200 ms before next stimulus
        respFlag = 0;
        while (audioOnset+params.stimLength + params.iti) - GetSecs >= 0.2
            [keyDown, secs, keyCode] = KbCheck;
            % if subject detected deviant
            if keyDown && find(keyCode) == detectKeyCode
                % response trigger
                if triggerFlag
                    fprintf(serialObj, params.trig.format, params.trig.response);
                end      
                % set response flag, store RT
                respFlag = 1;
                outData.rt(trialIdx, blockIdx) = secs-audioOnset;
                outData.rtDeviant(trialIdx, blockIdx) = secs-(audioOnset + outData.deviantOnset(trialIdx, blockIdx));
                break;
            % if experiment is aborted
            elseif keyDown && find(keyCode) == abortKeyCode
                abortFlag = 1;
                break;
            end
        end
        
        % if abort was requested, quit
        if abortFlag
            disp([char(10), char(10), 'USER REQUESTED ABORT!', char(10), char(10)]);
            Priority(0);
            RestrictKeysForKbCheck([]);
            PsychPortAudio('Stop', pahandle, 1, 1);
            PsychPortAudio('Close');
            return;
        end        
        
        % check accuracy
        accFlag = 0;
        if respFlag
            if deviantMask(trialIdx)
                outData.accuracy(trialIdx, blockIdx) = 1;
                accFlag = 1;
            else
                outData.accuracy(trialIdx, blockIdx) = 0;
            end
        elseif ~respFlag
            if ~deviantMask(trialIdx)
                outData.accuracy(trialIdx, blockIdx) = 1;
                accFlag = 1;
            else
                outData.accuracy(trialIdx, blockIdx) = 0;
            end
        end
        
        % Hit / False alarm trigger
        if triggerFlag
            % if response was a HIT
            if respFlag && accFlag
                fprintf(serialObj, params.trig.format, params.trig.hit);
            % if response was a FALSE ALARM
            elseif respFlag && ~accFlag
                fprintf(serialObj, params.trig.format, params.trig.falseAlarm);
            end 
        end  % if triggerFlag
        
        % audio onset becomes the last trial start
        oldTrialStartTime = audioOnset;
        
        % user messages
        if respFlag && accFlag
            disp([char(10), 'Subject pressed detection key. HIT!'])
            disp(['Subject pressed detection key at ', num2str(outData.rt(trialIdx, blockIdx)*1000), ' msec relative to stimulus onset, '...
                'and ', num2str(outData.rtDeviant(trialIdx, blockIdx)*1000), ' msec relative to deviant onset']);
        elseif respFlag && ~accFlag
            disp([char(10), 'Subject pressed detection key. FALSE ALARM!'])
            disp(['Subject pressed detection key at ', num2str(outData.rt(trialIdx, blockIdx)*1000), ' msec relative to stimulus onset']);
        end
        
        
    end  % for trialIdx
        
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Block-level feedback
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % user message
    disp([char(10), char(10), 'Done with trials, providing block-level feedback:', char(10)]);
    
    % Auditory task       
    % Determine TARGET-DEVIANT detection rates
    tmpIDmask = outData.trialTypes(:, blockIdx) == targetDeviantID;
    tmpCorr = sum(outData.accuracy(tmpIDmask, blockIdx));
    tmpTotal = sum(tmpIDmask);
    percentCorrectDetect = round(100*tmpCorr/tmpTotal);
    fprintf(['For **target deviant** stimuli... ' , num2str(percentCorrectDetect), '%% accuracy. \n'])

    % Determine NON-TARGET-DEVIANT false alarm rates
    tmpIDmask = outData.trialTypes(:, blockIdx) == nontargetDeviantID;
    tmpCorr = sum(~outData.accuracy(tmpIDmask, blockIdx));
    tmpTotal = sum(tmpIDmask);
    falseAlarmNonTargetDeviant = round(100*tmpCorr/tmpTotal);
    fprintf(['False alarm for non-target deviant stimuli... ' , num2str(falseAlarmNonTargetDeviant), '%%. \n'])

    % Determine TARGET-STANDARD false alarm rates
    tmpIDmask = outData.trialTypes(:, blockIdx) == (targetDeviantID-20);        % ------- NOTE  ----------->              This statement is very specific to the conditionID values specified in 'getnSNRdeviantStimulusBlock.m'!!!!!!!!!!!!!!!!!
    tmpCorr = sum(~outData.accuracy(tmpIDmask, blockIdx));
    tmpTotal = sum(tmpIDmask);
    falseAlarmTargetStandard = round(100*tmpCorr/tmpTotal);
    fprintf(['False alarm for target standard conditions... ' , num2str(falseAlarmTargetStandard), '%%. \n'])

    % Determine NON-TARGET-STANDARD false alarm rates
    tmpIDmask = outData.trialTypes(:, blockIdx) == (nontargetDeviantID-20);     % ------- NOTE  ----------->              This statement is very specific to the conditionID values specified in 'getnSNRdeviantStimulusBlock.m'!!!!!!!!!!!!!!!!!
    tmpCorr = sum(~outData.accuracy(tmpIDmask, blockIdx));
    tmpTotal = sum(tmpIDmask);
    falseAlarmNonTargetStandard = round(100*tmpCorr/tmpTotal);
    fprintf(['False alarm for non-target standard conditions... ' , num2str(falseAlarmNonTargetStandard), '%%. \n'])

    % Determine OVERALL ACCURACY rates
    tmpCorr = sum(outData.accuracy(:, blockIdx));
    tmpTotal = params.trialNo;
    percentALLaccuracy = round(100*tmpCorr/tmpTotal);
    fprintf(['Across ALL conditions... ' , num2str(percentALLaccuracy), '%% accuracy. \n'])

    % Determine OVERALL FALSE ALARM rates
    tmpIDmask = outData.trialTypes(:, blockIdx) == (nontargetDeviantID-20) | outData.trialTypes(:, blockIdx) == (targetDeviantID-20) | outData.trialTypes(:, blockIdx) == nontargetDeviantID;     % ------- NOTE  ----------->              This statement is very specific to the conditionID values specified in 'getnSNRdeviantStimulusBlock.m'!!!!!!!!!!!!!!!!!
    tmpCorr = sum(~outData.accuracy(tmpIDmask, blockIdx));
    tmpTotal = sum(tmpIDmask);
    falseAlarmALL = round(100*tmpCorr/tmpTotal);
    fprintf(['Overall false alarm rate for non-target-deviant conditions... ' , num2str(falseAlarmALL), '%%. \n'])

    % Determine D-PRIME... this applies the d' correction mentioned in 
    % Doelling and Poeppel (2015) PNAS E6233-E6242 described in their
    % Methods and Materials section under 'Behavioral Analysis'.
    % Ultimately, it is based off of Macmillan and Kaplan (1985)
    % Psychol. Bulletin, vol 98, pp. 185-199 
    if percentCorrectDetect == 100
        pHit = 1-1/(2*sum(deviantMask)); 
    elseif percentCorrectDetect == 0 
        pHit = 1/(2*sum(deviantMask));
    else
        pHit = percentCorrectDetect/100;
    end

    if falseAlarmTargetStandard == 0 
        pFA = 1/(2*sum(deviantMask)); 
    elseif falseAlarmTargetStandard == 100
        pFA = 1-1/(2*sum(deviantMask));
    else
        pFA = falseAlarmTargetStandard/100; 
    end

    dprimeResult = dprime(pHit, pFA);
    fprintf(['The d'' value for relevant "target" conditions... ' , num2str(dprimeResult), '. \n'])

    
    %----------------------------------------------------
    % Save the data between each block... BUT this is a 
    % temporary save without some final output information.
    % This file is overwritten on each save, including the
    % final save of the output data.
    %----------------------------------------------------
    disp([char(10), 'Temporary save...', char(10)]);
    % call script that does the actual saving that is specific to this particular experiment
    saveFcn_intensityFG(subjID, deviantType, expName, 'temporary',... 
        params, outData, allDeviantIndices, allStimParams, blockIdx);
      

end  % for blockIdx



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Save & cleanup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

saveFcn_intensityFG(subjID, deviantType, expName, 'final',... 
    params, outData, allDeviantIndices, allStimParams);

Priority(0);
RestrictKeysForKbCheck([]);
PsychPortAudio('Stop', pahandle, 1, 1);
PsychPortAudio('Close');                
                
fprintf('-----------------------------------------------\n')
fprintf('Experiment completed.  Thanks for your time!\n')
fprintf('-----------------------------------------------\n')


return
