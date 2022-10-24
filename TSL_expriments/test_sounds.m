% Options for the beep sounds

%%%%% (1) default beep
disp('Sound 1 --- press enter') ; pause ; 
beep on;
beep

%%%%% (2) sounds (Cedric) 

[y1,Fs1] = audioread('Beep.wav');
[y2,Fs2] = audioread('Cue.wav');
[y3,Fs3] = audioread('finish.wav');

% === Long sound
% disp('Sound 2 --- press enter') ; pause ; 
% sound(y1,Fs1);

disp('Sound 3 --- press enter') ; pause ; 
sound(y2, Fs2) ;

% === low-pitched sound
% disp('Sound 4 --- press enter') ; pause ; 
% sound(y3,Fs3);

%%%%% (2) ex
disp('Sound 5 --- press enter') ; pause ; 
%%
WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
WarnWave = [sin(1:.6:400)] ; %, sin(1:.7:400), sin(1:.4:400)];

WarnWave = [sin(1:0.2:400)] ; 
% increase the step: decrease time to reach cycles (2pi) and hence increase 
% frequency (higher-pitched sound)
% small step: low-pitched tone
% freq = n_samples_up_to_2pi
Fs = 22050 ; 
Audio = audioplayer(WarnWave, Fs);
play(Audio);
%pause(2) ; play(Audio) ; pause(2) ; play(Audio);
