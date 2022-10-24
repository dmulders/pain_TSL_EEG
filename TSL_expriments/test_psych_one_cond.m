function [ ] = test_psych_one_cond()
% Test the interface generated with the Psychtoolbox for the learning of TP
% during sequences of painful stimuli. 
%
% Reproduce the interface used in exp_one_condition (without commanding the
% stimulation device). 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
run_debug = 1 ; 
two_screens = 1 ; 
formatted_text = 0 ; % ok on stim, ko on XPS (text is cropped)
%ListenChar(0); % reset buffer

% Parameters %%%%%%%%%%%%%%%%%%%%
isi         = [1.9:0.05:2.1] ;  % [1.4:0.05:1.6];
IRI         = 4:8 ;             % inter-response-interval
resptime    = 10 ;              %
ask_conf    = 1 ;               % also ask confidence

cond_idx    = 1 ;               % condition (session)
subj_name   = 'test' ;          %
subj_idx    = 0 ;               %

TPs         = [0.5,0.3] ; 
trial       = 10 ;             %
pulse_dur   = 0.1 ;             % duration of laser pulses [s]

SR          = 1000 ;            % samplingrate
V0_temp     = 20 ;              % temp corresp. to 0V, as set in the laser
V10_temp    = 70 ;              % temp corresp. to 10V, as set in the laser
foreperiod  = 0.2 ;             % to see the laser temperature feedback 
                                % before stimulation. 
postperiod  = 0.2 ;             % to be able to record inputs (eg RT) after the stim
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fn_dir = ['./all_data/'] ; 
if ~exist(fn_dir, 'dir') % 0 or 7 if it exists
    mkdir(fn_dir)
end
fn_session = [fn_dir, 'Seq_subj', ...
    num2str(subj_idx),'_cond',num2str(cond_idx),'.mat'] ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ======================================================================= %
% == (1) Generate the sequence of stimuli and the timing of the questions. 
% ======================================================================= %
% ==*== Stimuli intensities
s = GenRandSeq(trial, TPs) ;
% s:    sequence with entries in {1,2} of length = # trials = trial.

% ==*== Timing of the questions.
if ask_conf
    qu_time = 2*resptime ; 
else
    qu_time = resptime ; 
end

[trials_resp, tot_trial_time, timings, tot_run_time] = sample_lists(...
    trial, isi, qu_time, IRI) ;

save(fn_session, '-v7.3', 's', 'trials_resp', 'tot_trial_time', ...
    'timings', 'tot_run_time','subj_name','foreperiod','postperiod') ;

% ======================================================================= %
% == (3) Setup Psychtoolbox
% ======================================================================= %

%% Setup visuals and keyboard
% COMMENT: 
% * positions are given in (xpos,ypos,width,height), where (xpos, ypos) is
% the position of the upper left corner
% * y-axis: goes from top to bottom!

if run_debug
    Screen('Preference', 'SkipSyncTests', 1) ; % no checks, for debugging
else
    Screen('Preference', 'SkipSyncTests', 0) ; % no checks, for debugging
end
% default settings for setting up Psychtoolbox
PsychDefaultSetup(2) ;
AssertOpenGL ;

% Choosing the display with the highest dislay number
% /!\ if 2 screens: should be plugged in BEFORE launching matlab and the
% psychtoolbox
% AND the second screen should be the MAIN DISPLAY!
screens = Screen('Screens') ;
screenNumber = max(screens) ;
% rect:
%   XPS:    3840 x 2160
%   Stim:   1366 x 768
%   Dell screen maxwell: 3200 x 2560
screen_dims = [1200,1000] ; % XPS
screen_dims = [600,400] ; % STIM


type = 0 ; % 0 only available for windows 
Screen('Preference','TextRenderer', type) ; 
fontSz = 20 ; % TO ADJUST 
fontSz_cross = 100 ; % fixation cross
Screen('Preference', 'DefaultFontSize',fontSz) ;

if run_debug
    if two_screens
        rect = Screen('Rect', screenNumber) ; 
        %[w, rect] = Screen('OpenWindow',screenNumber,[0,0,0]) ; % OK on correct screen, but not FULL! 
          
        [w, rect] = Screen('OpenWindow',screenNumber,[0,0,0],rect) ; 
        % OK on MAIN DISPLAY (always) --> need to set S2 as main display
    else
        [w, rect] = Screen('OpenWindow', screenNumber, 0,[100,100,screen_dims]) ; % debugging
    end 
    
    % check which font is used:
    %[oldFontName,oldFontNumber,oldTextStyle]=Screen('TextFont',w)
    
    Screen('TextStyle', w, 1) ;
    %0=normal,1=bold,2=italic,4=underline
    Screen('TextFont',w,'Times new roman') ;    
    % Screen(w,'TextFont','Arial') ; Screen(w,'TextSize',40);
else
    [w, rect] = Screen('OpenWindow', screenNumber, 0) ; % real TESTS
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

% Get the centre coordinate of the window
[xCenter, yCenter] = RectCenter(rect) ;

% Define black and white
white = WhiteIndex(w) ;
black = BlackIndex(w) ;
red = [255 0 0] ;
blue = [0 0 255] ;
gray = [125 125 125] ; 

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

fiveKey = KbName('5') ;% KbName('5%') ;
twoKey = KbName('2') ;% KbName('2@') ;
zeroKey = KbName('0') ; 


key_vas_down = twoKey ; %downKey ; % twoKey
key_vas_up = fiveKey ; %upKey ; % threeKey
key_vas_enter = zeroKey ; %enterKey ; %oneKey

if run_debug
    ListenChar(1);
    %ListenChar(2);
else
    ListenChar(1);
    
    % suppress echo to the command line for keypresses
    %ListenChar(2);
    % --> CANNOT DO ANYTHING IN MATLAB PROMPT!!!
end
% Tell the Psychtoolbox function “GetChar” to start or stop listening for keyboard input
% 0: turn off character listening and reset the buffer which holds the captured characters.
% 1: enable listening
% 2: enable listening, AND any output of keypresses to Matlab windows is suppressed
%  --> Matlab will be left with a dead keyboard until we press CTRL+C to reenable keyboard input.

%% MAKE VAS
% ----- Position of the title -----
xpos_text = xCenter-screen_width*0.25 ; %xCenter-230 %'center' % 1370
ypos_text = screen_height * 0.10 ; %512

% ----- Rectangle for the VAS ----- (x by y pixels)
% vasL = k*100 to reach the min and max of 0 and 100% in integer nb of
% press
vasL = ceil(0.2*screen_height/100)*100 %300 ; %ceil(0.2*screen_height/100)*100
baseRect = [0 0 15 vasL] ;

if ask_conf
    VAS_center = ypos_text + ceil(0.5*vasL) + 10*fontSz  ; % yCenter % y axis: towards bottom
    VAS_c_conf = VAS_center + vasL + 20*fontSz ;
    ypos_text_conf = VAS_c_conf - 10*fontSz - ceil(0.5*vasL) ; 
else
    VAS_center = yCenter ; 
end
vas.center = VAS_center ; 
vas.L = vasL ; 

% Make our vas rectangle coordinates
vasRect = CenterRectOnPointd(baseRect, xCenter, VAS_center);

% make slider
sbaseRect = [0 0 50 6] ; % slider: has a thickness of 6 pixels
% Set the intial position of the square to be in the centre of the screen
sliderX = xCenter ;
sliderY = VAS_center ; %yCenter ;
sliderRect = CenterRectOnPointd(sbaseRect, sliderX, sliderY) ;

% Set the amount we want our square to move on each button press
pixelsPerPress = vasL/100 ; %3 ; %vasL/100 ; % vasL == 100%

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

% Do initial flip...
Screen('Flip', w) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% FIND VAS min

yPositionIsBaseline = 0 ; % XPS: only ok if 0 (otherwise cropped text)
% 1: the “y” pen start location defines the base line of drawn text, 
% otherwise it defines the top of the drawn text
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
    % ATTENTION: y axis goes downwards from upper left croner
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

vas.min = sliderY ; 

Screen('Flip', w);
WaitSecs(0.5);

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

vas.max = sliderY 
Screen(w,'TextSize',fontSz_cross) ;
add_general_text(w, '+', xCenter, yCenter, formatted_text,1, white, black) ; 

Screen('Flip', w);

%%
% ======================================================================= %
% == (4) Start the sequence
% ======================================================================= % 
disp(['-- Press enter to start the sequence'])
pause ;

yshift_txt = 40 ; % to adapt

% if ExG_state ~=0
%     io64(ioObj,address,marker_scanner);   %mark ExG
%     io64(ioObj,address,0);
% end

Screen(w,'TextSize',fontSz_cross) ;
add_general_text(w, '+', xCenter, yCenter, formatted_text,1, white, black) ; 

tstart = GetSecs;

proba_same = NaN(trial,1) ;
all_RT = NaN(trial,1) ;
missed_resp = NaN(trial,1)  ;
prob_type = NaN(trial,1) ;
all_start_times = NaN(trial,1) ;

% confidence ratings
if ask_conf
    conf_proba = NaN(trial,1) ;
    all_RT_conf = NaN(trial,1) ;
    missed_conf = NaN(trial,1)  ;
end

for idx_stim=1:trial
    %%%%%%%%
    Screen(w,'TextSize',fontSz_cross) ;
    add_general_text(w, '+', xCenter, yCenter, formatted_text,1, white, black) ; 
    Screen('Flip', w);
    
%     if ExG_state ~=0
%         io64(ioObj,address,s(t));   %mark ExG
%         io64(ioObj,address,0);
%     end
    %%%%%%%%%
    
    curr_I = s(idx_stim) ; 
    
    % =================================================================== %
    % Stimuli
    % =================================================================== %
    t = GetSecs;
    all_start_times(idx_stim) = t ; 
    % ---- Pulse
    WaitSecs(pulse_dur+foreperiod+postperiod) ; % ====== instead of real stim
    disp(['Stim ',num2str(idx_stim),'/',num2str(trial), ' (',...
        num2str(curr_I),') delivered']) ; 
    
    % ---- Indicate that stimulus has been delivered --> move the laser
    % head
    beep ; %%%%%%%%%%%%%% NEED HEADPHONES
    
    % ---- ISI    
    WaitSecs('UntilTime', t + timings(idx_stim,2) ) ;
   
    % =================================================================== %
    % Behavioral questions
    % =================================================================== %
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(find(idx_stim==trials_resp,1))
        disp('Behavioral question') ; 

        t0 = GetSecs ;        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % PROBABILITY ESTIMATE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sliderY = VAS_center;
        %disp(curr_I) ;
        Screen(w,'TextSize',fontSz) ;
        str_title = 'ESTIMATE PROBABILITY' ;         
        if curr_I == 2 % i.e. H
            prob_type(idx_stim) = 0 ;
            top_str = 'H > H' ; bottom_str = 'H > L' ; 
        else
            prob_type(idx_stim) = 1 ;
            top_str = 'L > L' ; bottom_str = 'L > H' ; 
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
            add_general_text(w, str_title, xCenter-screen_width*0.2, ypos_text, formatted_text,1, white, black) ; 
            
            add_general_text(w, top_str, xCenter-30, VAS_center - vasL/2 - 2.8*yshift_txt, ...
                formatted_text,1, white, black) ; 
            add_general_text(w, bottom_str, xCenter-30, VAS_center + vasL/2 + yshift_txt, ...
                formatted_text,1, white, black) ; 

            add_general_text(w, '100%', xCenter+30, VAS_center - vasL/2-yshift_txt, ...
                formatted_text,1, gray, black) ; 
            add_general_text(w, '100%', xCenter+30, VAS_center + vasL/2-yshift_txt, ...
                formatted_text,1, gray, black) ; 
            add_general_text(w, '50%', xCenter+30, VAS_center - yshift_txt, ...
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
        sliderY = sliderY
        proba_same_tmp = proba_same(idx_stim)   
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % ASK CONFIDENCE IN THE RATINGS
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if ask_conf
            WaitSecs(0.5) ; % to avoid escaping the question because of previous press  
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
                add_general_text(w, '50%', xCenter+30, VAS_center - yshift_txt, ...
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
                Screen('Flip', w) ;                
            end
            
            conf_proba(idx_stim) = 100*((VAS_c_conf+vasL/2)-slider_conf)/(vasL) ;
            slider_conf = slider_conf
            conf_proba_tmp = conf_proba(idx_stim)
        end
        
        Screen(w,'TextSize',fontSz_cross) ;
        add_general_text(w, '+', xCenter, yCenter, formatted_text,1, white, black) ; 
        Screen('Flip', w);
        %WaitSecs('UntilTime', (trial_est_startime(idx_stim) + tot_trial_time(idx_stim)) );
        
    end
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TO ADD:
    % 'screenXpixels', 'screenYpixels', 'xCenter','yCenter','vas'
    
    save(fn_session, '-v7.3', 's', 'trials_resp', 'tot_trial_time', ...
        'timings', 'tot_run_time','subj_name','foreperiod','postperiod',...
        'proba_same', 'all_RT','missed_resp','prob_type','vas','all_start_times') ; 
end

end



