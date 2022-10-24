function [ ] = exp_one_condition(NI_session, TPs, trial, pulse_dur, temp_base, ...
    temp_I1, temp_I2, V0_temp, V10_temp, fn_session, subj_name, foreperiod,...
    postperiod, varargin)

% Generate one test condition for the probabilistic learning experiment. 
%
% Inputs:
%  TPs    (1,2) vector containing [p(1|2), p(2|1)]. 
%  fn_session   To save the behavioral data AND laser temp feedback
%
lab_screen = 1 ; 
beep on ; 
opt_beep = 2 ;                  % 1 = matlab beep ; 2 = custom beep
formatted_text = 0 ;            % depends on computer, to avoid cropped text
test_no_LSD = 0 ;               % test the code without triggering stimulation

% TCS2_BNC:
use_foreground = 0 ; add_beep = 1 ; 

ListenChar(0);                  % reset buffer, to avoid warning 
                                % "KbQueue event buffer is full! Maximum capacity..."
two_screens = 1 ;               % participants have a monitor in front of them
VAS_min_max = 1 ;               % ask to subject to encode min and max of VAS
% Parameters %%%%%%%%%%%%%%%%%%%%
isi         = [2.9:0.05:3.1] ;  % 3 pm 0.1
                                % SHOULD BE SUFFICIENT AS:
                                % LSD: need to have time to move the laser
                                % TCS: should be in contact with the skin
                                % since around 1 sec before pulse arrive to
                                % avoid A beta ERP at the time of stim
resptime    = 8 ;               % per question (*2 if 2 questions)
ask_conf    = 1 ;               % also ask confidence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFAULT 
IRI         = 12:18 ;           % inter-response-interval (15 pm 3 trials)
stim_type   = 1 ;               % 1: LSD, 2: TCS2
TCS_USB     = 0 ;               % USB or NO
active_zones= 1:5 ;             % TCS2 (stim_type = 2)
temp_slope  = 200 ;             % °C/s
trigger     = 255 ;             % pins
temp_feedback = 0 ;             % when using TCS_USB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


yPositionIsBaseline = 0 ;       % XPS: only ok if 0 (otherwise cropped text)
% 1: the “y” pen start location defines the base line of drawn text, 
% otherwise it defines the top of the drawn text

nargs_left = length(varargin) ;
if nargs_left > 0
    if ~(round(nargs_left/2) == nargs_left/2)
        error('-- There should be an even number of varargin arguments')
    end
    for i = 1:2:nargs_left
        Name_arg = varargin{i};
        Val_arg = varargin{i+1};
        if ~ischar(Name_arg)
            error('Flag arguments must be strings')
        end
        Name_arg = lower(Name_arg);
        switch Name_arg
            case 'iri'
                IRI = Val_arg ;
            case 'stim_type'
                stim_type = Val_arg ; 
            case 'tcs_usb'
                TCS_USB = Val_arg ; 
            case 'active_zones'
                active_zones = Val_arg ; 
            case 'temp_slope'
                temp_slope = Val_arg ; 
            case 'trigger'
                trigger = Val_arg ; 
            case 'temp_feedback'
                temp_feedback = Val_arg ; 
            case 'vas_min_max'
                VAS_min_max = Val_arg ; 
        end
    end
end

if TCS_USB
    SR = NaN ; 
elseif ~test_no_LSD
    SR = NI_session.Rate ; 
end

if ask_conf
    qu_time = 2*resptime ; 
else
    qu_time = resptime ; 
end

% ======================================================================= %
% == (1) Generate the sequence of stimuli and the timing of the questions. 
% ======================================================================= %
% ==*== Stimuli intensities
s = GenRandSeq(trial, TPs) ;
% s:    sequence with entries in {1,2} of length = # trials = trial.

% ==*== Timing of the questions. 
[trials_resp, tot_trial_time, timings, tot_run_time] = sample_lists(...
    trial, isi, qu_time, IRI) ;

save(fn_session, '-v7.3', 's', 'trials_resp', 'tot_trial_time', ...
    'timings', 'tot_run_time','subj_name','foreperiod','postperiod',...
    'temp_I1', 'temp_I2') ; 
% '-v7.3': to be able to save .mat file > 2 GB

% ======================================================================= %
% == (2) Get output waveforms
% ======================================================================= %
if ~test_no_LSD
    switch stim_type
        case 1
            % ======== LSD ======== %
            wave_I1 = get_waveform_LSD(temp_I1, temp_base, pulse_dur, SR,...
                foreperiod, postperiod, V0_temp, V10_temp) ; 
            wave_I2 = get_waveform_LSD(temp_I2, temp_base, pulse_dur, SR,...
                foreperiod, postperiod, V0_temp, V10_temp) ;
        case 2
            % ======== TCS ======== %
            if ~TCS_USB
                wave_I1 = get_waveform_TCS2(NI_session, temp_I1, temp_base, pulse_dur,...
                    foreperiod, postperiod, active_zones) ; 
                wave_I2 = get_waveform_TCS2(NI_session, temp_I2, temp_base, pulse_dur,...
                    foreperiod, postperiod, active_zones) ; 
            end
    end
end

% ======================================================================= %
% == (3) Setup Psychtoolbox
% ======================================================================= %
% Setup visuals and keyboard
Screen('Preference', 'SkipSyncTests', 1) ; % no checks, for debugging
%Screen('Preference', 'SkipSyncTests', 0) ;
% default settings for setting up Psychtoolbox
PsychDefaultSetup(2) ;
AssertOpenGL ;

% Choosing the display with the highest dislay number
screens = Screen('Screens') ;
screenNumber = max(screens) ;
% rect:
%   XPS:    3840 x 2160
%   Stim:   1366 x 768
screen_dims = [1200,1000] ; % XPS
screen_dims = [600,400] ; % STIM

type = 0 ; % 0 only available for windows 
Screen('Preference','TextRenderer', type) ; 
fontSz = 20 ; % TO ADJUST 
fontSz_cross = 100 ; % fixation cross
Screen('Preference', 'DefaultFontSize',fontSz) ;

if two_screens
    rect = Screen('Rect', screenNumber) ; 
    [w, rect] = Screen('OpenWindow',screenNumber,[0,0,0],rect) ; 
    % OK on MAIN DISPLAY (always) --> need to set S2 as main display
else
    [w, rect] = Screen('OpenWindow', screenNumber, 0,[100,100,screen_dims]) ;
end 

screen_width = rect(3) ; 
screen_height = rect(4) ; 

% TO CLOSE THE psych windows: sca ; 

% Enable alpha blending with proper blend-function. We need it
% for drawing of smoothed points:
Screen('BlendFunction', w, GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA) ;
fps = Screen('FrameRate',w) ;      % frames per second
ifi = Screen('GetFlipInterval', w) ;
if fps==0
    fps=1/ifi ;
end

% Get the size of the on screen window
[screen_width, screen_height] = Screen('WindowSize', w) ;

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(rect) ;

% Define black and white
white = WhiteIndex(w) ;
black = BlackIndex(w) ;
red = [255 0 0] ;
blue = [0 0 255] ;
gray = [125 125 125] ; 

%HideCursor; % Hide the mouse cursor (even on second display!!! --> never)

Priority(MaxPriority(w)) ;

% set keyboard
KbName('UnifyKeyNames') ; % http://psychtoolbox.org/docs/KbName
escapeKey = KbName('q') ;
oneKey = KbName('1!') ;
twoKey = KbName('2@') ;
threeKey = KbName('3#') ;
fourKey = KbName('4$') ;
sevenKey = KbName('7&') ;
%leftKey = KbName('LeftArrow');
%rightKey = KbName('RightArrow');
upKey = KbName('UpArrow');
downKey = KbName('DownArrow');
enterKey = KbName('Return') ;
% To get the name of each keycode: KbName('KeyNames') ; 
% returns a cell array

delKey = KbName('delete') ; % to pause sequence % 'clear'

fiveKey = KbName('5') ;% KbName('5%') ;
twoKey = KbName('2') ;% KbName('2@') ;
zeroKey = KbName('0') ; 

if test_no_LSD
    key_vas_down = downKey ; % twoKey
    key_vas_up = upKey ; % threeKey
    key_vas_enter = enterKey ; %oneKey
else
    key_vas_down = twoKey ; %downKey ; % twoKey
    key_vas_up = fiveKey ; %upKey ; % threeKey
    key_vas_enter = zeroKey ; %enterKey ; %oneKey
end

ListenChar(1);
%ListenChar(2);

% Tell the Psychtoolbox function “GetChar” to start or stop listening for keyboard input
% 0: turn off character listening and reset the buffer which holds the captured characters.
% 1: enable listening
% 2: enable listening, AND any output of keypresses to Matlab windows is suppressed
%  --> Matlab will be left with a dead keyboard until we press CTRL+C to reenable keyboard input.
%      suppress echo to the command line for keypresses
%      --> CANNOT DO ANYTHING IN MATLAB PROMPT!!!

%% MAKE VAS
% ----- Position of the title -----
xpos_text = xCenter-screen_width*0.25 ; %xCenter-230 %'center' % 1370
if lab_screen
    ypos_text = screen_height * 0.03 ; %0.05 ; %512
else
    ypos_text = screen_height * 0.10 ; %512
end

% ----- Rectangle for the VAS ----- (x by y pixels)
vasL = ceil(0.15*screen_height/100)*100 % 300
% to reach 0 and 100% in integer nb of press (1press=1%)
baseRect = [0 0 15 vasL] ;

if ask_conf
    if lab_screen
        % 8*fontSz
        VAS_center = ypos_text + ceil(0.5*vasL) + 5*fontSz  ; % yCenter % y axis: towards bottom
        VAS_c_conf = VAS_center + vasL + 10*fontSz ; % center of conf VAS
        ypos_text_conf = VAS_c_conf - 5*fontSz - ceil(0.5*vasL) ;
        % BEFORE BUG 2019-08-08
        %VAS_center = ypos_text + ceil(0.5*vasL) + 8*fontSz  ; % yCenter % y axis: towards bottom
        %VAS_c_conf = VAS_center + vasL + 16*fontSz ; % center of conf VAS
        %ypos_text_conf = VAS_c_conf - 8*fontSz - ceil(0.5*vasL) ;
    else
        VAS_center = ypos_text + ceil(0.5*vasL) + 10*fontSz  ; % yCenter % y axis: towards bottom
        VAS_c_conf = VAS_center + vasL + 20*fontSz ;
        ypos_text_conf = VAS_c_conf - 10*fontSz - ceil(0.5*vasL) ; 
    end
else
    VAS_center = yCenter ; 
end
vas.center = VAS_center ; 
vas.L = vasL ; 

% Make our vas rectangle coordinates
vasRect = CenterRectOnPointd(baseRect, xCenter, VAS_center);

% make slider
sbaseRect = [0 0 50 6] ; % slider: has a thickness of 6 pixels!!
% Set the intial position of the square to be in the centre of the screen
sliderX = xCenter ;
sliderY = VAS_center ;
sliderRect = CenterRectOnPointd(sbaseRect, sliderX, sliderY) ;
% Set the amount we want our square to move on each button press
pixelsPerPress = vasL/100 ; %3 ;

% tick
tick = [0 0 30 6] ;
tickRect_m = CenterRectOnPointd(tick, xCenter, VAS_center) ; % middle tick
tickRect_b = CenterRectOnPointd(tick, xCenter, VAS_center-(vasL/2)) ; % top tick
tickRect_t = CenterRectOnPointd(tick, xCenter, VAS_center+(vasL/2)) ; % bottom tick
if ask_conf
    tickRect_c_m = CenterRectOnPointd(tick, xCenter, VAS_c_conf) ; % middle tick
    tickRect_c_b = CenterRectOnPointd(tick, xCenter, VAS_c_conf-(vasL/2)) ; % top tick
    tickRect_c_t = CenterRectOnPointd(tick, xCenter, VAS_c_conf+(vasL/2)) ; % bottom tick
    
    vasRect_c = CenterRectOnPointd(baseRect, xCenter, VAS_c_conf);
    sliderRect_c = CenterRectOnPointd(sbaseRect, sliderX, VAS_c_conf) ;
end

%% GET VAS MIN MAX
if VAS_min_max
    % Do initial flip...
    Screen('Flip', w) ; % flip: to display on screen

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND VAS min

    str_title = 'Please select the LOWEST point in the scale and enter' ; 
    add_general_text(w, str_title, xpos_text, ypos_text, formatted_text,...
        1, white, black, yPositionIsBaseline) ; 

    Screen('FillRect', w, gray, tickRect_m);
    Screen('FillRect', w, gray, tickRect_b);
    Screen('FillRect', w, gray, tickRect_t);
    Screen('FillRect', w, gray, vasRect);
    Screen('FillRect', w, red, sliderRect);
    Screen('Flip', w);
    % This is the cue which determines whether we exit the vas
    exitvas = false;

    % Loop the animation until the escape key is pressed
    while exitvas == false

        % Check the keyboard to see if a button has been pressed
        [keyIsDown,secs, keyCode] = KbCheck ;
        %keyIsDown 1 if any key, including modifiers such as <shift>,<control> or <caps lock> is down. 0 otherwise.

        % Depending on the button press, move slider up and down
        % ATTENTION: y axis goes downwards from upper left corner
        if keyCode(key_vas_enter)
            exitvas = true;
        elseif keyCode(key_vas_down)
            sliderY = sliderY + pixelsPerPress;
        elseif keyCode(key_vas_up)
            sliderY = sliderY - pixelsPerPress;
        end

        % We set bounds to make sure our slider doesn't go completely off of
        % the VAS boundaries
        if sliderY < VAS_center-(vasL/2)
            sliderY = VAS_center-(vasL/2);
        elseif sliderY > VAS_center+(vasL/2)
            sliderY = VAS_center+(vasL/2);
        end
        sliderRect = CenterRectOnPointd(sbaseRect, sliderX, sliderY);
        % Draw to the screen
        add_general_text(w, str_title, xpos_text, ypos_text, formatted_text,...
            1, white, black, yPositionIsBaseline) ; 

        Screen('FillRect', w, gray, tickRect_m);
        Screen('FillRect', w, gray, tickRect_b);
        Screen('FillRect', w, gray, tickRect_t);
        Screen('FillRect', w, gray, vasRect);
        Screen('FillRect', w, red, sliderRect); % cursor
        % Indicate percentage
        val_percent = 50 + 50*abs(sliderY-VAS_center)/(vasL/2) ; 
        str_percent = num2str(val_percent) ; 
        add_general_text(w, str_percent, sliderX+0.8*sbaseRect(3), sliderY-sbaseRect(4), formatted_text,...
            0, white, black, yPositionIsBaseline) ; 

        % Flip to the screen
        Screen('Flip', w);
    end

    vas.min = sliderY ; % VAS_center + vasL/2

    Screen('Flip', w) ;
    WaitSecs(0.5) ; % to avoid escaping next VAS due to last key press

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND VAS max
    sliderY = VAS_center ;
    sliderRect = CenterRectOnPointd(sbaseRect, sliderX, sliderY);

    str_title = 'Now select the HIGHEST point in the scale and enter' ; 
    add_general_text(w, str_title, xpos_text, ypos_text, formatted_text,...
        1, white, black, yPositionIsBaseline) ; 

    Screen('FillRect', w, gray, tickRect_m);
    Screen('FillRect', w, gray, tickRect_b);
    Screen('FillRect', w, gray, tickRect_t);
    Screen('FillRect', w, gray, vasRect);
    Screen('FillRect', w, red, sliderRect);
    Screen('Flip', w);

    % This is the cue which determines whether we exit the vas
    exitvas = false;

    % Loop the animation until the escape key is pressed
    while exitvas == false

        % Check the keyboard to see if a button has been pressed
        [keyIsDown,secs, keyCode] = KbCheck;

        % Depending on the button press, move slider up and down
        if keyCode(key_vas_enter)
            exitvas = true;
        elseif keyCode(key_vas_down)
            sliderY = sliderY + pixelsPerPress;
        elseif keyCode(key_vas_up)
            sliderY = sliderY - pixelsPerPress;
        end

        % We set bounds to make sure our slider doesn't go completely off of
        % the VAS boundaries
        if sliderY < VAS_center-(vasL/2)
            sliderY = VAS_center-(vasL/2);
        elseif sliderY > VAS_center+(vasL/2)
            sliderY = VAS_center+(vasL/2);
        end
        sliderRect = CenterRectOnPointd(sbaseRect, sliderX, sliderY);

        % Draw to the screen
        add_general_text(w, str_title, xpos_text, ypos_text, formatted_text,...
            1, white, black, yPositionIsBaseline) ; 

        Screen('FillRect', w, gray, tickRect_m);
        Screen('FillRect', w, gray, tickRect_b);
        Screen('FillRect', w, gray, tickRect_t);
        Screen('FillRect', w, gray, vasRect);
        Screen('FillRect', w, red, sliderRect);

        val_percent = 50 + 50*abs(sliderY-VAS_center)/(vasL/2) ; 
        str_percent = num2str(val_percent) ; 
        add_general_text(w, str_percent, sliderX+0.8*sbaseRect(3), sliderY-sbaseRect(4), formatted_text,...
            0, white, black, yPositionIsBaseline) ; 

        % Flip to the screen
        Screen('Flip', w);
    end

    vas.max = sliderY ;
end
ypos_cross = 0.8*yCenter ; 
Screen(w,'TextSize',fontSz_cross) ;
add_general_text(w, '+', xCenter, ypos_cross, formatted_text,1, white, black) ; 

Screen('Flip', w);


% ======================================================================= %
% == (4) Start the sequence
% ======================================================================= %
if lab_screen
    yshift_txt = 16 ;
else
    yshift_txt = 40 ; % to adapt for the tick labels of the VAS
end

laser_data = repmat(struct(),trial,1) ; % store data and time
all_start_times = NaN(trial,1) ;
n_shots = zeros(trial,1) ; % nb of times each stim had to be tried to be effectively delivered

% Probability estimates
proba_same = NaN(trial,1) ;
all_RT = NaN(trial,1) ;
missed_resp = NaN(trial,1)  ;
prob_type = NaN(trial,1) ;

% confidence ratings
if ask_conf
    conf_proba = NaN(trial,1) ;
    all_RT_conf = NaN(trial,1) ;
    missed_conf = NaN(trial,1)  ;
end

disp(['-- Press enter to start the sequence'])
pause ;

Screen(w,'TextSize',fontSz_cross) ;
add_general_text(w, '+', xCenter, ypos_cross, formatted_text,1, white, black) ; 

tstart = GetSecs;

for idx_stim=1:trial
    
    Screen(w,'TextSize',fontSz_cross) ;
    add_general_text(w, '+', xCenter, ypos_cross, formatted_text,1, white, black) ; 
    Screen('Flip', w);
    
    % ============= PAUSE the squence???
    [keyIsDown,secs, keyCode] = KbCheck ;
    %keyIsDown 1 if any key is pressed, 0 otherwise.
    if keyCode(delKey)
        disp('Sequence paused, press enter to continue...')
        pause ; 
    end
    
    curr_I = s(idx_stim) ; 
    if test_no_LSD
        t = GetSecs;
        WaitSecs(pulse_dur+foreperiod+postperiod) ; % ====== instead of real stim
        play_custom_beep(opt_beep) ; 
    else
        switch stim_type
            case 1
                % =================================================================== %
                % Laser stimulation
                % =================================================================== %
                if curr_I==1
                    outputData = wave_I1 ;                    
                else %curr_I==2
                    outputData = wave_I2 ;
                end
                laser_shooted = 0 ; 
                while ~laser_shooted
                    % ---- Pulse
                    % we do not use [data,time] = stimulate_LSD(NI_session, outputData, 1) ; 
                    % for the timing...
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    queueOutputData(NI_session,outputData) ;
                    prepare(NI_session) ;
                    t = GetSecs;
                    [data, time] = startForeground(NI_session) ;

                    laser_data(idx_stim).data=data ;
                    laser_data(idx_stim).time=time ;

                    % out_slope = 10/(V10_temp-V0_temp) ;
                    % V_one_deg = out_slope ; % tension for one degree
                    cond_laser_ok = (max(data(:,1)) - min(data(:,1)))>1 ; % TO TEST?
                    %%%%%%%%%%% TO DO: check if laser shooted or not!!
                    if cond_laser_ok
                        laser_shooted = 1 ; 
                    else
                        disp('/!\ LASER DID NOT SHOOT')
                        disp('    ---> (1) Set laser to correct height')
                        disp('    ---> (2) Press doctor button')
                        disp('    ---> (3) Press enter to try again...')
                        pause ; 
                    end   
                    n_shots(idx_stim) = n_shots(idx_stim) + 1 ; 
                end
                % ---- Indicate that stimulus has been delivered
                %%%%%%%%%%%%%% NEED HEADPHONES
                play_custom_beep(opt_beep) ; 
            case 2
                % =================================================================== %
                % TCS2 stimulation
                % =================================================================== %
                if TCS_USB
                    if curr_I==1
                        curr_temp = temp_I1 ; 
                    else %curr_I==2
                        curr_temp = temp_I2 ;
                    end
                    flush_waveform_TCS2_USB(NI_session, curr_temp, pulse_dur, ...
                        active_zones, temp_slope, trigger, temp_feedback) ;
                    t = GetSecs;
                    data = stimulate_TCS2_USB(NI_session, temp_feedback) ;                     
                    WaitSecs('UntilTime', t + pulse_dur) ;  
                    play_custom_beep(opt_beep) ;
                    
                    laser_data(idx_stim).data=data ;
                else
                    if curr_I==1
                        outputData = wave_I1 ;
                    else %curr_I==2
                        outputData = wave_I2 ;
                    end
                    % ---- Pulse
                    % we do not use [data,time] = stimulate_TCS2(NI_session, outputData, 1) ; 
                    % for the timing...
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %queueOutputData(NI_session.NI,outputData) ;
                    %prepare(NI_session.NI) ;
                    %t = GetSecs ;
                    %[data,time] = NI_session.NI.startForeground() ;% startForeground(NI_session.NI)
                    %play_custom_beep(opt_beep) ;
                    
                    
                    %%%%%%%%%%%%% using background session
                    [data,time,t] = stimulate_TCS2(NI_session, outputData, ...
                        use_foreground, add_beep) ; 
                    laser_data(idx_stim).data=data ;
                    laser_data(idx_stim).time=time ;
                    
                end                
        end
    end
    
    disp(['=== Stim ',num2str(idx_stim),'/',num2str(trial),' delivered (',num2str(curr_I),') ==='])

    all_start_times(idx_stim) = t ; % only keep the real start time of each trial
    
    % ---- ISI    
    WaitSecs('UntilTime', t + timings(idx_stim,2) ) ; 
   
    % =================================================================== %
    % Behavioral questions
    % =================================================================== %
    if ~isempty(find(idx_stim==trials_resp,1))
        disp('Behavioral question') ; 

        t0 = GetSecs ;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PROBABILITY ESTIMATE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sliderY = VAS_center;
        Screen(w,'TextSize',fontSz) ;
        str_title = 'ESTIMATE PROBABILITY' ;         
        if curr_I == 2 % i.e. H
            prob_type(idx_stim) = 0 ;
            %top_str = 'H > H' ; bottom_str = 'H > L' ; 
            top_str = 'warm -> warm' ; bottom_str = 'warm -> cool' ;
        else
            prob_type(idx_stim) = 1 ;
            %top_str = 'L > L' ; bottom_str = 'L > H' ; 
            top_str = 'cool -> cool' ; bottom_str = 'cool -> warm' ;
        end
        
        % This is the cue which determines whether we exit the vas
        exitvas = false;
        missed_resp(idx_stim,1) = 1;
        
        % Loop the animation until the one key is pressed or time has
        % elapsed
        while exitvas == false
            
            % quit if max time elapses
            clock = GetSecs;
            t_spent = clock - t0 ; 
            if (t_spent > resptime), exitvas = true; end
            
            % Check the keyboard to see if a button has been pressed
            [keyIsDown,keyTime, keyCode] = KbCheck;
            
            if ((keyTime - t0) > resptime), exitvas = true; end
            
            % Depending on the button press, move slider up and down
            if keyCode(key_vas_enter)
                exitvas = true;
                all_RT(idx_stim,1) = keyTime - t0;
                missed_resp(idx_stim,1) = 0;
            elseif keyCode(key_vas_down)
                sliderY = sliderY + pixelsPerPress;
            elseif keyCode(key_vas_up)
                sliderY = sliderY - pixelsPerPress;
            end
            
            
            % We set bounds to make sure our slider doesn't go completely off of
            % the VAS boundaries
            if sliderY < VAS_center-(vasL/2)
                sliderY = VAS_center-(vasL/2);
            elseif sliderY > VAS_center+(vasL/2)
                sliderY = VAS_center+(vasL/2);
            end
            sliderRect = CenterRectOnPointd(sbaseRect, sliderX, sliderY);
            
            % Draw to the screen
            add_general_text(w, str_title, xCenter-screen_width*0.1, ypos_text, formatted_text,1, white, black) ; % was 0.2
            
            add_general_text(w, top_str, xCenter-30, VAS_center - vasL/2 - 2.8*yshift_txt, ...
                formatted_text,1, white, black) ; 
            add_general_text(w, bottom_str, xCenter-30, VAS_center + vasL/2 + yshift_txt, ...
                formatted_text,1, white, black) ; 

            add_general_text(w, '100%', xCenter+30, VAS_center - vasL/2-yshift_txt, ...
                formatted_text,1, gray, black) ; 
            add_general_text(w, '100%', xCenter+30, VAS_center + vasL/2-yshift_txt, ...
                formatted_text,1, gray, black) ; 
            add_general_text(w, '50-50%', xCenter+30, VAS_center - yshift_txt, ...
                formatted_text,1, gray, black) ; 

            Screen('FillRect', w, gray, tickRect_m);
            Screen('FillRect', w, gray, tickRect_b);
            Screen('FillRect', w, gray, tickRect_t);
            Screen('FillRect', w, gray, vasRect);
            Screen('FillRect', w, red, sliderRect);

            str_sec_remain = [num2str(max(0,floor(resptime-t_spent))),' sec'] ; 
            add_general_text(w, str_sec_remain, xCenter+0.2*screen_width, VAS_center, formatted_text,1, red, black) ; 
            
            % Flip to the screen
            Screen('Flip', w);
        end
        proba_same(idx_stim) = 100*((VAS_center+vasL/2)-sliderY)/(vasL) ; % sliderY ;
        %sliderY = sliderY
        %proba_same_tmp = proba_same(idx_stim) 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ASK CONFIDENCE IN THE RATINGS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ask_conf
            % increase time before confidence 
            WaitSecs(2) ; % to avoid escaping the question because of previous press  
            t0 = GetSecs ;  
            
            slider_conf = VAS_c_conf ; % REPLACE sliderY = VAS_center;
            str_tit_conf = 'ARE YOU SURE?' ; 
            top_conf = 'COMPLETELY' ; %'YES' ; %cfr Meyniel2015 fig1 p3.
            bottom_conf = 'NOT AT ALL' ; %'NO' ;

            
            exitvas = false;
            missed_conf(idx_stim,1) = 1 ;
            
            while exitvas == false
                % quit if max time elapses
                clock = GetSecs ;
                t_spent = clock - t0 ; 
                if (t_spent > resptime), exitvas = true; end

                % Check the keyboard to see if a button has been pressed
                [keyIsDown,keyTime, keyCode] = KbCheck;

                if ((keyTime - t0) > resptime), exitvas = true; end

                % Depending on the button press, move slider up and down
                if keyCode(key_vas_enter)
                    exitvas = true;
                    all_RT_conf(idx_stim,1) = keyTime - t0 ;
                    missed_conf(idx_stim,1) = 0;
                elseif keyCode(key_vas_down)
                    slider_conf = slider_conf + pixelsPerPress;
                elseif keyCode(key_vas_up)
                    slider_conf = slider_conf - pixelsPerPress;
                end
                
                if slider_conf < VAS_c_conf-(vasL/2)
                    slider_conf = VAS_c_conf-(vasL/2);
                elseif slider_conf > VAS_c_conf+(vasL/2)
                    slider_conf = VAS_c_conf+(vasL/2);
                end
                sliderRect_c = CenterRectOnPointd(sbaseRect, sliderX, slider_conf);
                
                % DRAW TO THE SCREEN
                
                % === (1) === Still show the proba estimate
                str_title = 'ESTIMATE PROBABILITY' ; 
                add_general_text(w, str_title, xCenter-screen_width*0.1, ypos_text, formatted_text,1, white, black) ; 
                add_general_text(w, top_str, xCenter-30, VAS_center - vasL/2 - 2.8*yshift_txt, ...
                    formatted_text,1, white, black) ; 
                add_general_text(w, bottom_str, xCenter-30, VAS_center + vasL/2 + yshift_txt, ...
                    formatted_text,1, white, black) ;
                add_general_text(w, '100%', xCenter+30, VAS_center - vasL/2-yshift_txt, ...
                    formatted_text,1, gray, black) ; 
                add_general_text(w, '100%', xCenter+30, VAS_center + vasL/2-yshift_txt, ...
                    formatted_text,1, gray, black) ; 
                add_general_text(w, '50-50%', xCenter+30, VAS_center - yshift_txt, ...
                    formatted_text,1, gray, black) ;
                Screen('FillRect', w, gray, tickRect_m);
                Screen('FillRect', w, gray, tickRect_b);
                Screen('FillRect', w, gray, tickRect_t);
                Screen('FillRect', w, gray, vasRect);
                Screen('FillRect', w, red, sliderRect);
                
                % === (2) === VAS for the confidence
                add_general_text(w, str_tit_conf, xCenter-screen_width*0.1, ypos_text_conf, formatted_text,1, white, black) ; 
                add_general_text(w, top_conf, xCenter-30, VAS_c_conf - vasL/2 - 2.8*yshift_txt, ...
                    formatted_text,1, white, black) ; 
                add_general_text(w, bottom_conf, xCenter-30, VAS_c_conf + vasL/2 + yshift_txt, ...
                    formatted_text,1, white, black) ; 
                add_general_text(w, '100%', xCenter+30, VAS_c_conf - vasL/2-yshift_txt, ...
                    formatted_text,1, gray, black) ; 
                add_general_text(w, '0%', xCenter+30, VAS_c_conf + vasL/2-yshift_txt, ...
                    formatted_text,1, gray, black) ; 
                add_general_text(w, '50%', xCenter+30, VAS_c_conf - yshift_txt, ...
                    formatted_text,1, gray, black) ; 

                Screen('FillRect', w, gray, tickRect_c_m);
                Screen('FillRect', w, gray, tickRect_c_b);
                Screen('FillRect', w, gray, tickRect_c_t);
                Screen('FillRect', w, gray, vasRect_c);
                Screen('FillRect', w, red, sliderRect_c);
                
                str_sec_remain = [num2str(max(0,floor(resptime-t_spent))),' sec'] ; 
                add_general_text(w, str_sec_remain, xCenter+0.2*screen_width, VAS_c_conf, formatted_text,1, red, black) ; 
                
                % Flip to the screen
                Screen('Flip', w);                
            end
            
            conf_proba(idx_stim) = 100*((VAS_c_conf+vasL/2)-slider_conf)/(vasL) ; 
            % slider_conf
        end
        Screen(w,'TextSize',fontSz_cross) ;
        add_general_text(w, '+', xCenter, ypos_cross, formatted_text,1, white, black) ; 
        Screen('Flip', w);
        %WaitSecs('UntilTime', (trial_est_startime(idx_stim) + tot_trial_time(idx_stim)) );
        WaitSecs(0.5) ; % before sending next pulse 
    end
    
    % =================================================================== %    
    
    save(fn_session, '-v7.3', 's', 'trials_resp', 'tot_trial_time', ...
        'timings', 'tot_run_time','subj_name','foreperiod','postperiod',...
        'laser_data','proba_same','all_RT','missed_resp','prob_type',...
        'vas','all_start_times','yCenter','xCenter','rect', 'tstart',...
        'conf_proba','missed_conf', 'all_RT_conf', 'temp_I1', 'temp_I2') ; 
    
end

end



