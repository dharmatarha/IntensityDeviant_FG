% [sigOut] = osullivanTokenIntensityDeviantBack(fs, cl, clNot, lenNote, lenNoteRamp, numNotesPerToken, minSemitoneStep, clITD, clNotITD, clTraj, [baseFreq, upFreq], maxSemitoneStep, clStart, [clPhaseFLAG, clNotPhaseFLAG], plotMe, deviantNum, deviantInd, deviantSNR)
%
% [sigOut, params] = osullivanTokenIntensityDeviantBack(...
% [sigOut, params, clOut] = osullivanTokenIntensityDeviantBack(...
% [sigOut, params, clOut, clNotOut] = osullivanTokenIntensityDeviantBack(...
%
%  ***NOTE*** This function can also generate the same stimuli as
%  'osullivanToken.m' by specifying the *exact same* input parameters.  It
%  is, therefore, necessary to use the "optional" parameters to generate
%  the deviant stimulus.
%
%--------------------------------------------------------------------------
% INPUTS:
%   fs              sampling frequency (in SAMPLES/SEC) 
%   cl              coherence level, i.e. # of note frequencies with common
%                   trajectory over time (specified as Integer). This
%                   can also be specified as an array of integers so that
%                   multiple 'cl' frequencies can be positioned at multiple
%                   ITD values... per the appropriate usage of 'clITD'.
%   clNot           number of frequencies (specified as Integer) not
%                   associated with the chord. This can also be specified
%                   as an array of integers so that multiple 'clNot'
%                   frequencies can be positioned at multiple ITD values...
%                   per the appropriate usage of 'clNotITD'.
%   lenNote         length of each individual note frequency (in SAMPLES)
%   lenNoteRamp     length of ramp applied to each note frequency (in
%                   SAMPLES) 
%   numNotesPerToken                        ...        
%                   number of *note-bursts* generated per token(an
%                   INTEGER). This naming convention could be slightly
%                   confusing since the previous input parameters use
%                   "note" as a reference to a given frequency, but "notes"
%                   for this input parameter refers to the number of
%                   temporally successive note-bursts used, i.e. this input
%                   parameter indirectly specifies the duration of the
%                   output signal.  
%   minSemitoneStep                         ...
%                   smallest semitone-step between possible notes.
%   clITD           ITD values (in MICROSECONDS) for each note of the
%                   *chord*. It is important that the dimensions of this
%                   array are correct!! Should be [1 x numNotesPerToken].
%                   If a value of NaN is used, the phase difference between
%                   the left and right channels is random. If 'cl' is
%                   specified as an array of integers, 'clITD' must be
%                   specified as a CELL ARRAY such that in each cell entry
%                   has an array of dimension [1 x numNotesPerToken].
%                   Therefore, different ITDs can be specified for
%                   different 'cl' frequencies. NOTE: each clTone and its
%                   respective frequency trajectory will maintain the same
%                   specified ITD throughout the duration of the token.
%   clNotITD        ITD values (in MICROSECONDS) for each note that is not
%                   associated with the chord. It is important that the
%                   dimensions of this array are correct!! Should be [1 x
%                   numNotesPerToken]. If a value of NaN is used, the phase
%                   difference between the left and right channels is
%                   random. If 'clNot' is specified as an array of
%                   integers, 'clNotITD' must be specified as a CELL ARRAY
%                   such that in each cell entry has an array of dimension
%                   [1 x numNotesPerToken].  Therefore, different ITDs can
%                   be specified for different 'clNot' frequencies.
%   clTraj          Specifies the frequency trajectory of the notes in the
%                   chord. Values are denoted by (INTEGER) steps. A step of
%                   1 indicates all notes in the chord increase by one
%                   'minSemitoneStep' in the next tone burst. Negative
%                   values can be used to decrease the trajectory. This
%                   input parameter can also be specified as an empty
%                   matrix, i.e. [].  This will cause the trajectory to be
%                   random within the bounds of [minSemitoneStep,
%                   maxSemitoneStep].  If the array is not empty, the
%                   dimensions of this array is imporant!! Should be [1 x
%                   numNotesPerToken-1].  Note the (-1) because only N-1
%                   changes in frequency trajectory are required. It is OK
%                   if 'clTraj' is specified as [1 x numNotesPerToken], but
%                   the last value of the array will not be used.
%   [baseFreq, upFreq]                      ...
%                   (OPTIONAL) define the lower and upper frequencies (Hz)
%                   that are available for generating the token. Default
%                   values are [100, 2000].
%   maxSemitoneStep                         ...
%                   (OPTIONAL) define the largest semitone-step between
%                   temporally adjacent notes of a chord. This value is
%                   only important when isempty(clTraj). The default value
%                   is 2.
%   clStart         (OPTIONAL) define the initial frequency indices for the
%                   chord. In this way, the chord can be specified to
%                   always begin at certain frequencies and at certain
%                   frequency separation between components of the chord.
%                   The default value is [] which indicates a randomized
%                   selection of the frequencies of the first chord. IFF
%                   multiple ITD values for 'cl' are specified, i.e.
%                   length(cl) > 1, then this input variable should be
%                   specified as a CELL array where the length of each cell
%                   entry has a length equal to the corresponding 'cl'
%                   value, e.g. if cl = [5,3,4] then length(clStart{1}) ==
%                   5, length(clStart{2})==3, and length(clStart{3})==4.
%   [clPhaseFLAG, clNotPhaseFLAG]                   ...
%                   (OPTIONAL) Specify the *starting phase* for the 'cl'
%                   and the 'clNot' tokens. If ~0, the starting phase
%                   across all tones is random.  If 0, the frequencies all
%                   begin at zero (sine-phase). Default values are [1 1].
%   plotMe          (OPTIONAL) setting this to a non-zero value will result
%                   in the spectrogram of the left channel of the primary
%                   stimulus being plotted, i.e. sigOut(:,1)
%----  Inputs related to the generation of the "deviant"
%   deviantNum      (OPTIONAL) Integer number of tones that are associated
%                   with the deviant. This value must be less than or equal
%                   to sum(clNot).  Note that this input parameter can be
%                   specified as an array where the size of the array
%                   corresponds to the number of ITD values for the 'clNot'
%                   tokens.  
%   deviantInd      (OPTIONAL) This is an array of integers where the
%                   deviant should occur.  This value is restricted to be
%                   less than or equal to 'numNotesPerToken'.
%   deviantSNR      (OPTIONAL) change in deviant - to - non-deviant power
%                   ratio in dB. Positive value indicates the deviant tones
%                   are greater in power than the remaining background
%                   tones. NOTE: although this input parameter maintains
%                   the same name as in 'IntensityDeviantFig.m', the value
%                   represents a different property of the stimulus.
%
% OUTPUTS:
%   sigOut          output token with an RMS of 1.
%   params          (OPTIONAL) structure with the relevant parameters to
%                   the stimulus. 
%   clOut           (OPTIONAL) stereo signal for the notes only part of the
%                   chord. 
%   clNotOut        (OPTIONAL) stereo signal for the notes not associated
%                   with the chord.
%
% USAGE:
%   ----->  Refer to 'osullivanToken.m' for some examples for using the
%           obligatory input parameters.
%
% ******EXAMPLE 1: single ITD used for 'clNot' (with no values for 'cl') 
% fs = 48000;
% numNotesPerToken = 20;
% lenNote = 0.120; %1/13;
% cl = 0; 
% clITD = 0*ones(1, numNotesPerToken);     
% clNot = 10; 
% clNotITD = 0*ones(1, numNotesPerToken);    
% clTraj = 0*ones(1, numNotesPerToken-1);  % this can also be empty [] to generate a random trajectory
% clStart = []; % use [] to specify random starting notes
% baseFreq = 100;
% upFreq = 2000;
% minSemitoneStep = 0.5;
% plotMe = 1;
% 
% deviantNum = 3;  
% deviantInd = [10:10];
% deviantSNR = -25;
% 
% [sigOut] = osullivanTokenIntensityDeviantBack(fs, cl, clNot,  round(lenNote*fs),  round(0.010*fs),  numNotesPerToken, minSemitoneStep, clITD, clNotITD, clTraj, [baseFreq upFreq], 2, clStart, [1 1], plotMe, deviantNum, deviantInd, deviantSNR);   
% 
%
% ******EXAMPLE 2: multiple ITD values used for 'clNot' (with no values for 'cl') 
% fs = 48000;
% numNotesPerToken = 20;
% lenNote = 0.120; %1/13;
% cl = [0 0 0];
% clITD = {nan(1, numNotesPerToken), 400*ones(1, numNotesPerToken), zeros(1,numNotesPerToken)};    %-400*ones(1, numNotesPerToken);    
% clNot = [0 4 6];
% clNotITD = {nan(1, numNotesPerToken), 400*ones(1, numNotesPerToken), zeros(1,numNotesPerToken)};    
% clTraj = 0*ones(1, numNotesPerToken-1);  % this can also be empty [] to generate a random trajectory
% clStart = []; %{[15; 21; 35; 43], [53;59;57], [25; 31]};   
% baseFreq = 100;
% upFreq = 2000;
% minSemitoneStep = 0.5;
% plotMe = 1;
% 
% deviantNum = [0 4 0];  
% deviantInd = [10:14];
% deviantSNR = -25;
% 
% [sigOut] = osullivanTokenIntensityDeviantBack(fs, cl, clNot,  round(lenNote*fs),  round(0.010*fs),  numNotesPerToken, minSemitoneStep, clITD, clNotITD, clTraj, [baseFreq upFreq], 2, clStart, [1 1], plotMe, deviantNum, deviantInd, deviantSNR);   
%
% 
%-----------------------------------------------
% Created by Darrin K. Reed on January 22, 2016.
%       Stimuli is a modification of my stimulus 'osullivanToken.m' where
%       the last modifcation ocurred on October 4, 2015.
% Last modified by Darrin K. Reed on January 22, 2016   
%
%**************************************************************************

function [sigOut, varargout] = osullivanTokenIntensityDeviantBack(fs, cl, clNot, lenNote, lenNoteRamp, numNotesPerToken, minSemitoneStep, clITD, clNotITD, clTraj, varargin)

%--------------------------------------------------------------------------
% Handle variable number of input arguments and ensure appropriate values
%--------------------------------------------------------------------------
switch nargin
    case 10
        baseFreq = 100;  % lower frequency bound (in Hz)
        upFreq = 2000;  % upper frequency bound (in Hz)
        maxSemitoneStep =  2;  % largest semitone step for successive notes in chord (in semitones)
        clStart = [];
        clPhaseFLAG = 1;
        clNotPhaseFLAG = 1;
        plotMe = 0;
        deviantNum = zeros(1, length(clNot));
        deviantInd = [];
        deviantSNR = 0;
        
    case 11
        if size(varargin{1},2) ~= 2
            error(['DKR: argument #', num2str(nargin), ' must be specified as [''baseFreq'', ''upFreq''].'])
        end

        baseFreq = varargin{1}(1);  % lower frequency bound (in Hz)
        upFreq = varargin{1}(2);  % upper frequency bound (in Hz)
        maxSemitoneStep =  2;  % largest semitone step for successive notes in chord (in semitones)
        clStart = [];
        clPhaseFLAG = 1;
        clNotPhaseFLAG = 1;
        plotMe = 0;
        deviantNum = zeros(1, length(clNot));
        deviantInd = [];
        deviantSNR = 0;
    
    case 12
        if size(varargin{1},2) ~= 2
            error(['DKR: argument #', num2str(nargin-1), ' must be specified as [''baseFreq'', ''upFreq''].'])
        end
        baseFreq = varargin{1}(1);  % lower frequency bound (in Hz)
        upFreq = varargin{1}(2);  % upper frequency bound (in Hz)
        maxSemitoneStep =  varargin{2};  % largest semitone step for successive notes in chord (in semitones)
        clStart = [];
        clPhaseFLAG = 1;
        clNotPhaseFLAG = 1;
        plotMe = 0;
        deviantNum = zeros(1, length(clNot));
        deviantInd = [];
        deviantSNR = 0;
        
    case 13
        if size(varargin{1},2) ~= 2
            error(['DKR: argument #', num2str(nargin-1), ' must be specified as [''baseFreq'', ''upFreq''].'])
        end
        baseFreq = varargin{1}(1);  % lower frequency bound (in Hz)
        upFreq = varargin{1}(2);  % upper frequency bound (in Hz)
        maxSemitoneStep =  varargin{2};  % largest semitone step for successive notes in chord (in semitones)
        clStart = varargin{3};
        clPhaseFLAG = 1;
        clNotPhaseFLAG = 1;
        plotMe = 0;
        deviantNum = zeros(1, length(clNot));
        deviantInd = [];
        deviantSNR = 0;
    
    case 14
        if size(varargin{1},2) ~= 2
            error(['DKR: argument #', num2str(nargin-1), ' must be specified as [''baseFreq'', ''upFreq''].'])
        end
        baseFreq = varargin{1}(1);  % lower frequency bound (in Hz)
        upFreq = varargin{1}(2);  % upper frequency bound (in Hz)
        maxSemitoneStep =  varargin{2};  % largest semitone step for successive notes in chord (in semitones)
        clStart = varargin{3};
        clPhaseFLAG = varargin{4}(1);
        clNotPhaseFLAG = varargin{4}(2);
        plotMe = 0;
        deviantNum = zeros(1, length(clNot));
        deviantInd = [];
        deviantSNR = 0;
        
    case 15
        if size(varargin{1},2) ~= 2
            error(['DKR: argument #', num2str(nargin-1), ' must be specified as [''baseFreq'', ''upFreq''].'])
        end
        baseFreq = varargin{1}(1);  % lower frequency bound (in Hz)
        upFreq = varargin{1}(2);  % upper frequency bound (in Hz)
        maxSemitoneStep =  varargin{2};  % largest semitone step for successive notes in chord (in semitones)
        clStart = varargin{3};
        clPhaseFLAG = varargin{4}(1);
        clNotPhaseFLAG = varargin{4}(2);
        plotMe = varargin{5};
        deviantNum = zeros(1, length(clNot));
        deviantInd = [];
        deviantSNR = 0;
        
    case 18
        if size(varargin{1},2) ~= 2
            error(['DKR: argument #', num2str(nargin-1), ' must be specified as [''baseFreq'', ''upFreq''].'])
        end
        baseFreq = varargin{1}(1);  % lower frequency bound (in Hz)
        upFreq = varargin{1}(2);  % upper frequency bound (in Hz)
        maxSemitoneStep =  varargin{2};  % largest semitone step for successive notes in chord (in semitones)
        clStart = varargin{3};
        clPhaseFLAG = varargin{4}(1);
        clNotPhaseFLAG = varargin{4}(2);
        plotMe = varargin{5};
        deviantNum = varargin{6};
        deviantInd = varargin{7};
        deviantSNR = varargin{8};
        
    otherwise
        error('DKR: invalid number of input arguments.')
end


%--------------------------------------------------------------------------
% INITIAL PARAMS CHECK... ensure appropriate starting parameters
%--------------------------------------------------------------------------
if baseFreq >= upFreq
    error('DKR:  ''baseFreq'' must be less than ''upFreq''')
end

if isempty(clTraj) && (round(maxSemitoneStep/minSemitoneStep) ~= maxSemitoneStep/minSemitoneStep)
    error('DKR: need integer step size for note trajectory')
end


if ~iscell(clITD) && length(cl) > 1
    error('DKR: if multiple values are used for ''cl'', a corresponding CELL array must be used for ''clITD''.')
end
if iscell(clITD) && (length(clITD) ~= length(cl))
    error('DKR: if multiple values are used for ''cl'', there must be a corresponding array of ITD values in ''clITD'' for each ''cl'' entry.... and vice versa.')
end
if iscell(clITD) && (size(clITD{1},1) ~= 1 || size(clITD{1},2) ~= numNotesPerToken) 
    error('DKR: each cell in ''clITD'' must be of dimensions [1 x numNotesPerToken]')
end
if ~iscell(clITD) && (size(clITD,1) ~= 1 || size(clITD,2) ~= numNotesPerToken)
    error('DKR: dimensions of ''clITD'' must be of dimensions [1 x numNotesPerToken]')
end



if ~iscell(clNotITD) && length(clNot) > 1
    error('DKR: if multiple values are used for ''clNot'', a corresponding CELL array must be used for ''clNotITD''.')
end
if iscell(clNotITD) && (length(clNotITD) ~= length(clNot))
    error('DKR: if multiple values are used for ''clNot'', there must be a corresponding array of ITD values in ''clNotITD'' for each ''clNot'' entry... and vice versa.')
end
if iscell(clNotITD) && (size(clNotITD{1},1) ~= 1 || size(clNotITD{1},2) ~= numNotesPerToken) 
    error('DKR: each cell in ''clNotITD'' must be of dimensions [1 x numNotesPerToken]')
end
if ~iscell(clNotITD) && (size(clNotITD,1) ~= 1 || size(clNotITD,2) ~= numNotesPerToken)
    error('DKR: dimensions of ''clNotITD'' must be of dimensions [1 x numNotesPerToken]')
end

if ~isempty(clTraj)
    if (size(clTraj,1) ~= 1 || ~(size(clTraj,2) == numNotesPerToken || size(clTraj,2) == numNotesPerToken-1))
        error('DKR: dimensions of ''clTraj'' must be [1 x numNotesPerToken-1] or [1 x numNotesPerToken]')
    end
end

%---------------- Checks and handling for 'clStart'
% CHECK if 'clStart' is *effectively* empty, i.e. are all entries of
% clStart{ind} empty?
if iscell(clStart)
    
    % Loop through to see if all entries of the CELL array are empty matrices
    tmp = [];
    for i = 1:length(clStart)  
        try
            tmp = [tmp, clStart{i}];
        catch
            tmp = [tmp; clStart{i}];
        end
    end
    
    if isempty(tmp)  % if ALL cell arrays are EMPTY
        clStartEMPTY = 1;
        clStart = [];   % Force clStart to a single empty matrix
    else
        clStartEMPTY = 0;
    end  
else  % only a double array
    if isempty(clStart)
        clStartEMPTY = 1;
    else
        clStartEMPTY = 0;
    end
end
   

if clStartEMPTY == 0 % Run checks if 'clStart' is being specified
    if length(cl) == 1
        if ~iscell(clStart) % Force 'clStart' to a cell array if only one ITD is specified for the cl values
            tmp = clStart; clear clStart;
            clStart{1} = tmp;
        end
        
        if length(clStart{1}) ~= cl
            error('DKR: ''clStart'' must have a length equal to ''cl'' OR must be specified as an empty array.')
        end
    else
        if ~iscell(clStart) 
            error('DKR: If multiple ITD values are specified for ''cl'' (i.e. length(cl)~=1), ''clStart'' must be specified as a CELL array or as an EMPTY matrix.')
        end
        
        if length(clStart) ~= length(cl)
            error('DKR: An array of starting values must be specified for *every* ITD used for the cl values.')
        end
        
        for i =1:length(clStart)
            if length(clStart{i}) ~= cl(i)
                 error('DKR: ''clStart'' must have a length equal to ''cl'' for *every* ITD used for the cl values.')
            end            
        end
    end
end

if ~(nargout <= 4)
    error('DKR: the maximum number of output arguments is 4.')
end

% 'maxSemitoneStep' is only relevant when the trajectory is random. Set to
% NaN for less confusion in 'params' output structure 
if ~isempty(clTraj)  
   maxSemitoneStep = nan;  
end


% Ensure proper dimensions for the 'deviantNum' array
if length(deviantNum) ~= length(clNot)
    error('DKR: Dimensions of ''clNot'' and ''deviantNum'' must agree.')
end

% Ensure that there are enought clNot tones for the deviant specification
if sum(clNot<deviantNum)
    error('DKR: insufficient ''clNot'' tones for the specified ''deviantNum''')
end

% Warn users if they specify unusual parameters
if sum(clNot - deviantNum) == 0
    warning('DKR: Because ''clNot'' and ''deviantNum'' are equal, there is effectively NO DEVIANT with the current settings.')
end


%--------------------------------------------------------------------------
% Calculate useful values in advance
%--------------------------------------------------------------------------
winNote = windowMe(lenNoteRamp, lenNote-2*lenNoteRamp, lenNoteRamp);  % generate the window to be applied to each note

tmp = 0:minSemitoneStep:125; % number of semitones used for 'notesPossible' (integer)
notesPossible = baseFreq*2.^(tmp/12);  % possible notes to be used (in Hz)
notesPossible = notesPossible(notesPossible<=upFreq);  % limit the possible notes below an upper frequency

numNotesPossible = length(notesPossible);  % number of possible frequencies that can be selected for a note
t = 0:1/fs:(lenNote-1)/fs;  % time indices for a single note (i.e. not the entire token)
numSampsPerNote = length(t);  


%--------------------------------------------------------------------------
% IFF the user wants to also specify cl tones, need to handle calculate an
% appropriate ratio to sum the cl and clNot signals together
%--------------------------------------------------------------------------
if ~(sum(cl) == 0 || sum(clNot) == 0)
    N = 10;  % Number of iterations to be used for the estimate
    [clRatio, ~, ~] = calcExpectedCLtoCLnotRatio(N, sum(cl), sum(clNot), notesPossible, lenNote/fs, fs);
end



%--------------------------------------------------------------------------
% Handle the condition when user specifies the indicies for the first note
%--------------------------------------------------------------------------
if clStartEMPTY == 0
    % Determine max/min values for 'clStart', noting that this depends
    % on the size of 'cl'
    if length(cl) == 1
        clStartMAX = max(clStart{1});  % recall that clStart is forced to a CELL array in the parameter checks up above
        clStartMIN = min(clStart{1});    
        
    else  % cycle through each CELL of clStart and determine the max/min freq-index values specified
        clStartMAX = max(clStart{1});
        clStartMIN = min(clStart{1}); 
        for i = 2:length(cl)  
            tmpMAX = max(clStart{i});
            tmpMIN = min(clStart{i}); 
            
            if isempty(clStartMAX)  % Need this awkward handling in case multiple empty values are used for 'clStart'
                clStartMAX = tmpMAX;
            elseif ~isempty(tmpMAX)
                if tmpMAX > clStartMAX
                    clStartMAX = tmpMAX;
                end
            end
            
            if isempty(clStartMIN)  % Need this awkward handling in case multiple empty values are used for 'clStart'
                clStartMIN = tmpMIN; 
            elseif ~isempty(tmpMIN)
                if tmpMIN < clStartMIN
                    clStartMIN = tmpMIN; 
                end
            end
            
        end
    end

    
    if clStartMAX > numNotesPossible || clStartMIN <= 0
        error('DKR: The specified starting indices for the first note exceed the frequency limits.')
    end
end



%--------------------------------------------------------------------------
% Generate *random* note trajectory if trajectory is not specified.
% Although the trajectory is random, the largest step size is limited by
% 'maxSemitoneStep'. Note that having 'cjTraj' with dimension
% 'numNotesPerToken' is fine.  
%--------------------------------------------------------------------------
if isempty(clTraj)
    clTraj = randi([-maxSemitoneStep/minSemitoneStep,+maxSemitoneStep/minSemitoneStep], 1, numNotesPerToken-1); % select the trajectory for all cl notes after the starting chord
    
    %----------------------------------------------------------------------
    % If the indices of the first note are specified, then the random
    % trajectory must be adequate, i.e. not exceed the maximum/minimum
    % frequencies available. Check to ensure values are sufficient.
    %----------------------------------------------------------------------
    if clStartEMPTY == 0    
        % Calculate the minimum and maximum deviation as a result of 'clTraj'.
        % This will help determine if the randomly chosen 'clTraj' exceeds
        % the frequency bounds throughout the duration of the token as a
        % result of the specified starting indicies.
        clTraj_maxPos = max(cumsum(clTraj));
        clTraj_minPos = min(cumsum(clTraj));

        if clTraj_maxPos < 0  % force negative max values to zero since values of 'clTraj_maxPos' are interpreted as positive values in subsequent calculations
            clTraj_maxPos = 0;
        end

        if clTraj_minPos > 0  % force positive min values to zero since values of 'clTraj_minPos' are interpreted as negative values in subsequent calculations
            clTraj_minPos = 0;
        end
        
        % Re-calculate the randomized 'clTraj' array until it is sufficient
        % for the specified starting indicies of 'clStart'
        safeCheck = 1;
        while (clStartMAX + clTraj_maxPos > numNotesPossible) || (clStartMIN  + clTraj_minPos <= 0)
            clTraj = randi([-maxSemitoneStep/minSemitoneStep,+maxSemitoneStep/minSemitoneStep], 1, numNotesPerToken-1); % select the trajectory for all cl notes after the starting chord
            
            % Calculate the minimum and maximum deviation as a result of 'clTraj'.
            clTraj_maxPos = max(cumsum(clTraj));
            clTraj_minPos = min(cumsum(clTraj));

            if clTraj_maxPos < 0  % force negative max values to zero since values of 'clTraj_maxPos' are interpreted as positive values in subsequent calculations
                clTraj_maxPos = 0;
            end

            if clTraj_minPos > 0  % force positive min values to zero since values of 'clTraj_minPos' are interpreted as negative values in subsequent calculations
                clTraj_minPos = 0;
            end
           
            if safeCheck > 200000
                fprintf('After 200000 iterations, an adequate randomized ''clTraj'' cannot be established as a result of the values specified in ''clStart''\n.')
                fprintf('As a result, the function was given keyboard command.  If you wish to continue the attempt at finding\n')
                fprintf('a possible randomized ''clTraj'', type ''safeCheck = 0'', press <enter>, and then type ''return''.\n')
                keyboard
            end
            safeCheck = safeCheck+1;
            
        end
    end
end

% Calculate the minimum and maximum deviation as a result of 'clTraj'.
% This will help "revise" the starting tokens used so that the notes will
% not be out of bounds throughout the duration of the token.
clTraj_maxPos = max(cumsum(clTraj));
clTraj_minPos = min(cumsum(clTraj));

if clTraj_maxPos < 0  % force negative max values to zero since values of 'clTraj_maxPos' are interpreted as positive values in subsequent calculations
    clTraj_maxPos = 0;
end

if clTraj_minPos > 0  % force positive min values to zero since values of 'clTraj_minPos' are interpreted as negative values in subsequent calculations
    clTraj_minPos = 0;
end


%--------------------------------------------------------------------------
% Generate and select randomized index values for the starting note.
%--------------------------------------------------------------------------
if clStartEMPTY == 1
    % generate random order for all possible notes
    tmp = randperm(numNotesPossible);

    % lower/upper index bounds (based on 'clTraj' are used to prevent
    % out of bound CL index values.
    tmp = tmp(tmp > abs(clTraj_minPos) & tmp < (numNotesPossible-clTraj_maxPos+1));
    
    
    % select the note *indices* used for the first chord
    if length(cl)==1
        clStartSpacing{1} = tmp(1:cl); 
    else
        % Calculate indicies if multiple ITDs for 'cl' are specified
        for i = 1:length(cl)
            hi = sum(cl(1:i));
            lw = hi - cl(i) + 1;
            clStartSpacing{i} = tmp(lw:hi);      
        end
    end
        
else
    % Use values calculated above to check if maximum/minimum values
    % specified by 'clStart' are in violation of the specified 'clTraj'    
    if (clStartMIN + clTraj_minPos) <= 0
        error('DKR: Values defined by ''clTraj'' and ''clStart'' are invalid... exceed the LOWER limit of possible frequencies.')
    end
    
    if (clStartMAX + clTraj_maxPos) > numNotesPossible
        error('DKR: Values defined by ''clTraj'' and ''clStart'' are invalid... exceed the UPPER limit of possible frequencies.')
    end
    
    
    clStartSpacing = clStart;   % 'clStart' was specified by the user
    % RECALL: clStart was forced to being a CELL array in the checks up
    % above
end


%--------------------------------------------------------------------------
% Generate indices of notes/frequencies used in the stimulus. 
%--------------------------------------------------------------------------

%%% Generate index values for the CL tokens.
if sum(cl) ~= 0  
    clIndexVals = [];  % Initialize an empty array to add the index values to
    for k = 1:length(cl)
        
        if cl(k) > 0  % Only need to do this if actual values are requested for CL, i.e. CL ~= 0
            % initialize temporary array
            tmp = zeros(cl(k),numNotesPerToken); 

            % specify the indices for the first chord of the token for the
            % k-th entry of CL values
            tmp(:,1) = (clStartSpacing{k})';

            % Loop through to get all indices subsequent to starting indices
            for i = 2:numNotesPerToken  
                tmp(:,i) = tmp(:,i-1) + clTraj(1,i-1); 
            end

            % Store the index values generated for the k-th entry of CL values
            clIndexVals = [clIndexVals; tmp]; 
        end
    end
end

%%% Generate index values for the CLnot tokens.
if sum(clNot) ~= 0 
    clNotIndexVals = zeros(sum(clNot), numNotesPerToken);  % Initialize array for non-CL index values
    
    for i = 1:numNotesPerToken  
        tmpPossibleNonCL = randperm(numNotesPossible);
        clNotIndexVals(:,i) = tmpPossibleNonCL(1,1:sum(clNot));   % select a sufficient number of non-CL value to account for ALL non-CL index values
    end
    
    if sum(cl) == 0
        clStartSpacing = nan;  % just for the 'params' output array
    end
end


%--------------------------------------------------------------------------
% Generate a TABLE of interaural phase values corresponding to each of the
% specified ITD values for each possible note frequency
%--------------------------------------------------------------------------
tmp = repmat(notesPossible', 1, numNotesPerToken);  % replicate the 'notesPossible' array for easier computation of the TABLE

% Calculate the TABLE for the CL values
if ~iscell(clITD)  % if only one value for 'cl' is specified
    clITDarray = ones(numNotesPossible, 1) * clITD;
    ipdCL_TABLE = 2*pi*(clITDarray*1e-6) .* tmp;
else  % multiple values for 'clNot' were specified
    for i = 1:length(cl)  % for each 'clNot' value specified, send corresponding ITD values to a CELL ARRAY.
        clITDarray = ones(numNotesPossible, 1) * clITD{i};
        ipdCL_TABLE{i} = 2*pi*(clITDarray*1e-6) .* tmp;
    end
end

% Calculate the TABLE for the CL-NOT values
if ~iscell(clNotITD)  % if only one value for 'clNot' is specified
    clNotITDarray = ones(numNotesPossible, 1) * clNotITD;
    ipdCLnot_TABLE = 2*pi*(clNotITDarray*1e-6) .* tmp;
else  % multiple values for 'clNot' were specified
    for i = 1:length(clNot)  % for each 'clNot' value specified, send corresponding ITD values to a CELL ARRAY.
        clNotITDarray = ones(numNotesPossible, 1) * clNotITD{i};
        ipdCLnot_TABLE{i} = 2*pi*(clNotITDarray*1e-6) .* tmp;
    end
end



%--------------------------------------------------------------------------
% Generate notes for one token
%--------------------------------------------------------------------------
sigOutL = zeros(numSampsPerNote*numNotesPerToken,1);  
sigOutR = zeros(numSampsPerNote*numNotesPerToken,1);  
clOutL = zeros(numSampsPerNote*numNotesPerToken,1);  
clOutR = zeros(numSampsPerNote*numNotesPerToken,1); 
clNotOutL = zeros(numSampsPerNote*numNotesPerToken,1);  
clNotOutR = zeros(numSampsPerNote*numNotesPerToken,1); 

for j = 1:numNotesPerToken    
    
    % Initialize the note-burst for the CL and CLnot bursts to zeros
    tmpNoteL_cl = zeros(1, numSampsPerNote);
    tmpNoteR_cl = zeros(1, numSampsPerNote);
    tmpNoteL_clNot = zeros(1, numSampsPerNote);
    tmpNoteR_clNot = zeros(1, numSampsPerNote);
    tmpNoteL_clNotDEVIANT = zeros(1, numSampsPerNote);
    tmpNoteR_clNotDEVIANT = zeros(1, numSampsPerNote);
    
    %%% Generate notes associated with the chord for the j-th noise burst
    if sum(cl) ~= 0 
        
        if ~iscell(clITD)  % if only one value for 'cl' is specified
            for i = 1:cl  
                noteF = notesPossible(clIndexVals(i,j));
                if clPhaseFLAG ~= 0  %clNotPhaseFLAG = 1;
                    notePhase = rand(1)*2*pi; 
                else
                    notePhase = 0*2*pi; 
                end
                notePhaseDiff = ipdCL_TABLE(clIndexVals(i,j), j);
                if isnan(notePhaseDiff)
                   notePhaseDiff = rand(1)*2*pi; 
                end
                tmpNoteL_cl = tmpNoteL_cl + sin(2*pi*noteF.*t + notePhase);
                tmpNoteR_cl = tmpNoteR_cl + sin(2*pi*noteF.*t + (notePhase + notePhaseDiff));
            end
            
            clEndSpacing = clIndexVals(:,end);

            
        else  % multiple values for 'cl' were specified
            ctr =  1; % need a counter to cycle through each 'cl' value used, i.e. cycle through sum(cl)number of times
            
            % cycle through each 'cl' entry and generate the j-th CL tone burst
            for i = 1:length(cl)
                tmpIPD = ipdCL_TABLE{i};  % the dimensions of ipdCL_TABLE{i} are [numNotesPossible, numNotesPerToken]
                
                for tmp = 1:cl(i)                  
                    noteF = notesPossible(clIndexVals(ctr,j));
                    if clPhaseFLAG ~= 0
                        notePhase = rand(1)*2*pi; 
                    else
                        notePhase = 0*2*pi; 
                    end   
                    notePhaseDiff = tmpIPD(clIndexVals(ctr,j), j); 
                    if isnan(notePhaseDiff)
                       notePhaseDiff = rand(1)*2*pi; 
                    end
                    tmpNoteL_cl = tmpNoteL_cl + sin(2*pi*noteF.*t + notePhase);
                    tmpNoteR_cl = tmpNoteR_cl + sin(2*pi*noteF.*t + (notePhase + notePhaseDiff));
                    
                    ctr = ctr+1;
                end
            end
            
            clEndSpacing = clIndexVals(:,end);
        end
    else
        clEndSpacing = nan;  % just for the 'params' output array
    end
    
    
    %%% Generate notes not associated with the chord for the j-th noise burst
    if sum(clNot) ~= 0  
        
        if ~iscell(clNotITD)  % if only one value for 'clNot' is specified
            
            % Need to generate an array of indicies for the tones that should be
            % designated as deviants
            if deviantNum ~= 0
                deviantNumArray = 1:deviantNum;
            else
                deviantNumArray = 0;
            end
            
            
            for i = 1:clNot
                noteF = notesPossible(clNotIndexVals(i,j));
                if clNotPhaseFLAG ~= 0
                    notePhase = rand(1)*2*pi; 
                else
                    notePhase = 0*2*pi; 
                end  
                notePhaseDiff = ipdCLnot_TABLE(clNotIndexVals(i,j), j); 
                if isnan(notePhaseDiff)
                   notePhaseDiff = rand(1)*2*pi; 
                end
                
                %**********************************************************
                % Handle the segregation of clNot tokens into deviant and
                % non-deviant categories
                %**********************************************************
                if ~sum(j == deviantInd)
                    tmpNoteL_clNot = tmpNoteL_clNot + sin(2*pi*noteF.*t + notePhase);
                    tmpNoteR_clNot = tmpNoteR_clNot + sin(2*pi*noteF.*t + (notePhase + notePhaseDiff));
                else
                    if ~sum(i == deviantNumArray)
                        tmpNoteL_clNot = tmpNoteL_clNot + sin(2*pi*noteF.*t + notePhase);
                        tmpNoteR_clNot = tmpNoteR_clNot + sin(2*pi*noteF.*t + (notePhase + notePhaseDiff));
                    else 
                        tmpNoteL_clNotDEVIANT = tmpNoteL_clNotDEVIANT + sin(2*pi*noteF.*t + notePhase);
                        tmpNoteR_clNotDEVIANT = tmpNoteR_clNotDEVIANT + sin(2*pi*noteF.*t + (notePhase + notePhaseDiff));
                    end
                end
            end
            
        else  % multiple values for 'clNot' were specified
            ctr =  1; % need a counter to cycle through each 'clNot' value used, i.e. cycle through sum(clNot)number of times
            for i = 1:length(clNot)
                tmpIPD = ipdCLnot_TABLE{i};  % cycle through each 'clNot' entry  and generate the j-th CLnot tone burst
                
                % Need to generate an array of indicies for the tones that should be
                % designated as deviants
                if deviantNum(i) ~= 0
                    deviantNumArray = 1:deviantNum(i);
                else
                    deviantNumArray = 0;
                end
                
                for tmp = 1:clNot(i)
                    noteF = notesPossible(clNotIndexVals(ctr,j));  
                    if clNotPhaseFLAG ~= 0
                        notePhase = rand(1)*2*pi; 
                    else
                        notePhase = 0*2*pi; 
                    end  

                    notePhaseDiff = tmpIPD(clNotIndexVals(ctr,j), j); 
                    if isnan(notePhaseDiff)
                       notePhaseDiff = rand(1)*2*pi; 
                    end
                    
                    
                    %**********************************************************
                    % Handle the segregation of clNot tokens into deviant and
                    % non-deviant categories
                    %**********************************************************
                    if ~sum(j == deviantInd)
                        tmpNoteL_clNot = tmpNoteL_clNot + sin(2*pi*noteF.*t + notePhase);
                        tmpNoteR_clNot = tmpNoteR_clNot + sin(2*pi*noteF.*t + (notePhase + notePhaseDiff));
                    else
                        if ~sum(tmp == deviantNumArray)
                            tmpNoteL_clNot = tmpNoteL_clNot + sin(2*pi*noteF.*t + notePhase);
                            tmpNoteR_clNot = tmpNoteR_clNot + sin(2*pi*noteF.*t + (notePhase + notePhaseDiff));
                        else 
                            tmpNoteL_clNotDEVIANT = tmpNoteL_clNotDEVIANT + sin(2*pi*noteF.*t + notePhase);
                            tmpNoteR_clNotDEVIANT = tmpNoteR_clNotDEVIANT + sin(2*pi*noteF.*t + (notePhase + notePhaseDiff));
                        end
                    end
                    
                    ctr = ctr+1;
                    
                end
            end
        end
    end
    
    
    % Normalize tokens so that the RMS values are set to one... assuming
    % that the tokens are not all zeros.
    if sum(cl) ~= 0
        tmpNoteL_cl = tmpNoteL_cl/sqrt(mean(tmpNoteL_cl.^2));
        tmpNoteR_cl = tmpNoteR_cl/sqrt(mean(tmpNoteR_cl.^2));
    end
    
    if sum(clNot) ~= 0 && sum(clNot - deviantNum) > 0
        tmpNoteL_clNot = tmpNoteL_clNot/sqrt(mean(tmpNoteL_clNot.^2));
        tmpNoteR_clNot = tmpNoteR_clNot/sqrt(mean(tmpNoteR_clNot.^2));
    end
      
    % Handle the DEVIANT scaling and summation.
    if sum(j == deviantInd) && (mean(tmpNoteL_clNotDEVIANT.^2) ~= 0)
        % SCALE the deviant tones to an RMS of 1
        tmpNoteL_clNotDEVIANT = tmpNoteL_clNotDEVIANT/sqrt(mean(tmpNoteL_clNotDEVIANT.^2));
        tmpNoteR_clNotDEVIANT = tmpNoteR_clNotDEVIANT/sqrt(mean(tmpNoteR_clNotDEVIANT.^2));
        
        % SCALE the deviant tones by the deviant SNR
        tmpNoteL_clNotDEVIANT = tmpNoteL_clNotDEVIANT * 10^(deviantSNR/20); 
        tmpNoteR_clNotDEVIANT = tmpNoteR_clNotDEVIANT * 10^(deviantSNR/20);
        
        % COMBINE deviant and non-deviant tones
        tmpNoteL_clNot = tmpNoteL_clNotDEVIANT + tmpNoteL_clNot;
        tmpNoteR_clNot = tmpNoteR_clNotDEVIANT + tmpNoteR_clNot;
        
        % RESCALE the summed stimulus so that any subsequent summation with
        % the cl signal can be added at an appropriate ratio
        tmpNoteL_clNot = tmpNoteL_clNot/sqrt(mean(tmpNoteL_clNot.^2));
        tmpNoteR_clNot = tmpNoteR_clNot/sqrt(mean(tmpNoteR_clNot.^2)); 
    end
    
    
    % Scale the CL tokens to achieve desired CL-to-CLnot ratio...
    % assuming both the CL and CLnot stimuli are non-zero.
    if ~(sum(cl) == 0 || sum(clNot) == 0)
        tmpNoteL_cl = tmpNoteL_cl * 10^(clRatio/20);
        tmpNoteR_cl = tmpNoteR_cl * 10^(clRatio/20);
    end
    
    
    % sum the CL and CLnot note-bursts into a single output note-burst
    tmpNoteL = tmpNoteL_cl + tmpNoteL_clNot;
    tmpNoteR = tmpNoteR_cl + tmpNoteR_clNot;
    
    % Ensure the summed signal is always set to an RMS of one
    tmpNoteL = tmpNoteL/sqrt(mean(tmpNoteL.^2));
    tmpNoteR = tmpNoteR/sqrt(mean(tmpNoteR.^2));
    
    
    % Apply the window to the note-bursts
    tmpNoteL = tmpNoteL' .* winNote;  % window the single note
    tmpNoteR = tmpNoteR' .* winNote;  % window the single note
    tmpNoteL_cl = tmpNoteL_cl' .* winNote;  % window the single note
    tmpNoteR_cl = tmpNoteR_cl' .* winNote;  % window the single note
    tmpNoteL_clNot = tmpNoteL_clNot' .* winNote;  % window the single note
    tmpNoteR_clNot = tmpNoteR_clNot' .* winNote;  % window the single note
    
    % Insert the note-burst into the output array(s)
    indLow = numSampsPerNote*j - numSampsPerNote + 1;
    indHigh = numSampsPerNote*j;
    sigOutL(indLow:indHigh,1) = tmpNoteL;
    sigOutR(indLow:indHigh,1) = tmpNoteR;
    
    if nargout == 3 || nargout == 4
        clOutL(indLow:indHigh,1) = tmpNoteL_cl;
        clOutR(indLow:indHigh,1) = tmpNoteR_cl;
        clNotOutL(indLow:indHigh,1) = tmpNoteL_clNot;
        clNotOutR(indLow:indHigh,1) = tmpNoteR_clNot;
    end

        
end

%--------------------------------------------------------------------------
% Final things to send the signal(s) out
%--------------------------------------------------------------------------
% Force RMS of the output token(s) to 1
sigOutL = sigOutL/sqrt(mean(sigOutL.^2));
sigOutR = sigOutR/sqrt(mean(sigOutR.^2));
sigOutL = sigOutL/sqrt(mean(sigOutL.^2));
sigOutR = sigOutR/sqrt(mean(sigOutR.^2));
sigOut = [sigOutL, sigOutR];

if nargout == 3 || nargout == 4
    clOutL = clOutL/sqrt(mean(clOutL.^2));
    clOutR = clOutR/sqrt(mean(clOutR.^2));
    clNotOutL = clNotOutL/sqrt(mean(clNotOutL.^2));
    clNotOutR = clNotOutR/sqrt(mean(clNotOutR.^2));
    
    clOut = [clOutL, clOutR];
    clNotOut = [clNotOutL, clNotOutR];
end

%--------------------------------------------------------------------------
% Plot the left channel of the primary output signal
%--------------------------------------------------------------------------
if plotMe ~= 0
    winLen = pow2(nextpow2(lenNote)-1);
    [~,f,t,p] = spectrogram(sigOut(:,1),winLen,winLen/2,winLen,fs); % Display the spectrogram
    surf(t,f,10*log10(abs(p)),'EdgeColor','none'); 
    axis xy; axis tight; colormap(jet); view(0,90);
    ylim([baseFreq upFreq])
    xlabel('Time (s)');
    ylabel('Frequency (Hz)');
end


%--------------------------------------------------------------------------
%  Put important signal params into a structure for requested output
%--------------------------------------------------------------------------
params.fs = fs;
params.cl = cl;
params.clNot = clNot;
params.lenNote = lenNote;
params.lenNoteRamp = lenNoteRamp;
params.numNotesPerToken = numNotesPerToken;
params.clITD = clITD;
params.clNotITD = clNotITD;
params.clTraj = clTraj;
params.baseFreq = baseFreq;
params.upFreq = upFreq;
params.minSemitoneStep = minSemitoneStep;
params.maxSemitoneStep = maxSemitoneStep;  
params.clPhaseFLAG = clPhaseFLAG;
params.clNotPhaseFLAG = clNotPhaseFLAG;
%-----------
params.notesPossible = notesPossible;
params.clStartSpacing = clStartSpacing;
params.clEndSpacing = clEndSpacing;  % added by DKR 10/4/2015
params.sigOut = sigOut;


%--------------------------------------------------------------------------
%  Handle variable number of output arguments 
%--------------------------------------------------------------------------
switch nargout
    case 2
        varargout{1} = params;
    case 3
        params.clOut = clOut;
        varargout{1} = params;
        varargout{2} = clOut;
    case 4
        params.clOut = clOut;
        params.clNotOut = clNotOut;
        varargout{1} = params;
        varargout{2} = clOut;
        varargout{3} = clNotOut;
end



%**************************************************************************
% [myWin] = windowMe(lenAttack, lenHold, lenDecay)
%      OR
% [myWin] = windowMe(lenAttack, lenHold, lenDecay, winType)
%
% USAGE:
%   [myWin] = windowMe(0.3*fs, 1*fs, 0.1*fs);
%      OR
%   [myWin] = windowMe(3000, 200, 1000, @gausswin);
%
%
% INPUTS:
%   lenAttack   -   Length (in SAMPLES) of attack portion.  
%   lenHold     -   Length (in SAMPLES) of hold portion.  
%   lenDecay    -   Length (in SAMPLES) of decay portion.  
%   winType     -   (optional)  Choose the type of window shape.  Use the
%                   same shapes as the 'window.m' function.  See available
%                   options below.  (Default is @hann)
%
%  Possible winType choices:
%     @bartlett       - Bartlett window.
%     @barthannwin    - Modified Bartlett-Hanning window. 
%     @blackman       - Blackman window.
%     @blackmanharris - Minimum 4-term Blackman-Harris window.
%     @bohmanwin      - Bohman window.
%     @chebwin        - Chebyshev window.
%     @flattopwin     - Flat Top window.
%     @gausswin       - Gaussian window.
%     @hamming        - Hamming window.
%     @hann           - Hann window.
%     @kaiser         - Kaiser window.
%     @nuttallwin     - Nuttall defined minimum 4-term Blackman-Harris window.
%     @parzenwin      - Parzen (de la Valle-Poussin) window.
%     @rectwin        - Rectangular window.
%     @taylorwin      - Taylor window.
%     @tukeywin       - Tukey window.
%     @triang         - Triangular window.

function [myWin] = windowMe(lenAttack, lenHold, lenDecay, varargin)

if nargin == 3
    winType = @hann;
elseif nargin == 4
    winType = varargin{1};
else
    error('DRK error: invalid number of input arguments')
end

% create the full windows for the attack, hold and decay
atk = window(winType, 2*lenAttack);
dec = window(winType, 2*lenDecay);
hld = ones(lenHold,1); 

% select the desired portion of the attack and decay
atk = atk(1:lenAttack,1);
dec = dec(lenDecay+1:end,:);


myWin = [atk; hld; dec];  % combine to make the full window

%eof