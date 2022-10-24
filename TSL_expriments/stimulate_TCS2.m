function [data,time, t] = stimulate_TCS2(TCS2, volt_mat, ...
    use_foreground, add_beep)
%
% In:
%   TCS2            structure with params and initialized daq using 
%                   initialize_TCS2. 
%   volt_mat        \in (n_time_step, 8)
%   use_foreground  binary variable indicating to acquire data in
%                   foreground or background (both can cause unexpected bug
%                   with old versions of matlab, at the end of the
%                   stimulations and before saving the data). 

if nargin<4
    add_beep = 0 ; 
end
opt_beep = 2 ; % 1 = matlab beep ; 2 = custom beep
if nargin<3
   use_foreground = 0 ;  
end

% fetch samplingrate
SR = TCS2.Rate ;

NI_session = TCS2.NI ; 

if use_foreground
    % === Foreground acquisition and generation
    %queue the data to NI
    queueOutputData(NI_session,volt_mat);
    disp('Data sent to NI');

    %prepare LSD
    prepare(NI_session);
    %disp('NI prepared');

    if nargout>2
        t = GetSecs ;
    end
    %start stimulation
    [data,time]=NI_session.startForeground();
    
    if add_beep
        play_custom_beep(opt_beep) ;
    end 
    disp('Trial finished') ;
else
    duration_bins = size(volt_mat,1) ;
    % === Background acquisition and generation
    % continuous background generation: add a listener event to queue 
    % additional data to be outputted.

    %LSD_params.NI.IsContinuous = true ;  % ======= ? 
    fid1 = fopen('log.bin', 'w') ; 

    % add listener for the input (acquired data)
    lh_in=addlistener(NI_session,'DataAvailable',...
        @(src,event) logData(src,event, fid1)) ;
    % the handle function specified for the listener cannot increment a global
    % variable ==> have to write temporary results in a file (fid1)

    % add listener for the output (data sent to NI)
    % LSD_params.NI.IsContinuous = true; % in principle, should not be added...
    lh_out = addlistener(NI_session,'DataRequired', ...
            @(src,event) src.queueOutputData(volt_mat));

    %queue the data to NI
    queueOutputData(NI_session,volt_mat);
    %disp('Data sent to NI');

    %prepare LSD % optional
    prepare(NI_session);
    %disp('NI prepared');
    % start stimulation
    if nargout>2
        t = GetSecs ;
    end
    NI_session.startBackground();
    pause(0.25) ; 
    %disp('Waiting for NI to finish');

    % as long as the NI is busy, pause in the foreground 
    % (busy waiting, to be able to see when trial is finished)
    %while LSD_params.NI.IsRunning;% isDone % isLogging
    t_pause = 0.3 ;
    pause(t_pause) ; 
    if add_beep
        play_custom_beep(opt_beep) ;
    end    
    pause(duration_bins/SR-t_pause)
%     while ~LSD_params.NI.IsDone
%         pause(0.25);
%     end;

    NI_session.stop ; 
    
    % Cannot set IsDone, IsRunning nor IsLogging
   
    % in principle only required when LSD_params.NI.IsContinuous but depends on the bug
    delete(lh_in);
    delete(lh_out);
    disp('NI finished');

    fclose(fid1) ; 

    % Read the data of the current session to output them
    fid2 = fopen('log.bin', 'r') ;
    [raw_data, count] = fread(fid2, [9, Inf], 'double') ;
    % [9, Inf] = size of data to read
    fclose(fid2) ;

    % same format as when foreground acquisition is used
    time = raw_data(1,:) ; 
    time = time' ; 
    data = raw_data(2:end,:) ; 
    data = data' ; 
end

end

