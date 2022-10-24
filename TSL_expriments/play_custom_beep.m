function play_custom_beep(opt_beep) 
% Produce beep sound.

% for opt_beep = 1: need to call beep on ; before.

switch opt_beep
    case 1
        beep ; 
    case 2
        % cfr test_sounds
        
        WarnWave = sin(1:.2:400) ;
        Audio = audioplayer(WarnWave, 22050) ;
        %play(Audio);
        playblocking(Audio) ; 
end
% pause(0.5) ; % needed if play

end

