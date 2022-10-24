function TCS2 = initialize_TCS2_NI(samplingrate)
% This function requires two NI devices (Dev{i} and Dev{j}, i \neq j 
% in {1,2,3}).
% Allows to control the 5 stimulation zones of TCS2 separately. 
%
%Analog OUTPUTS on DEV1 : 
%DEV1/AO1 : neutral temperature
%DEV1/AO2 : active temperature zone 1
%DEV1/AO3 : active temperature zone 2
%DEV1/AO4 : active temperature zone 3
%
%Analog OUTPUTS on DEV2 : 
%DEV2/AO1 : active temperature zone 4
%DEV2/AO2 : active temperature zone 5
%DEV2/AO3 : neutral temperature trigger
%DEV2/AO4 : active temperature trigger
%
%Analog INPUTS on DEV1 : 4
%Analog INPUTS on DEV2 : 4

% Parameters %%%%%%%%%%%%             
NI_numbers = {'1', '3'} ;   % Specify the numbers of the 2 NI devices used 
                        % cfr names of the NI (Dev1, Dev2 or Dev3).
% {'4_1', '6'}
if nargin<1                        
    samplingrate = 1000;    % set the daq sampling rate
end
%%%%%%%%%%%%%%%%%%%%%%%%%

%nbstr_NI1 = num2str(NI_numbers(1)) ; 
%nbstr_NI2 = num2str(NI_numbers(2)) ; 
nbstr_NI1 = NI_numbers{1} ; 
nbstr_NI2 = NI_numbers{2} ; 

%create the session
disp('Creating NI session. This may take a few seconds.');
TCS2.NI=daq.createSession('ni');

%device names
dev1 = ['Dev', nbstr_NI1] ; 
dev2 = ['Dev', nbstr_NI2] ;
disp(['ID of first NI device : ' dev1]);
disp(['ID of second NI device : ' dev2]);

%trigger connections
trigger1 = ['Dev', nbstr_NI1, '/PFI4'] ;
trigger2 = ['Dev', nbstr_NI2, '/PFI0'] ;
disp(['Trigger connection on Device 1 : ' trigger1]);
disp(['Trigger connection on Device 2 : ' trigger2]);
disp('Trigger connections must be physically connected using BNC cable');
TCS2.trigger1=trigger1;
TCS2.trigger2=trigger2;

%clock connections
clock1 = ['Dev', nbstr_NI1, '/PFI5'] ;
clock2 = ['Dev', nbstr_NI2, '/PFI1'] ;
disp(['Clock connection on Device 1 : ' clock1]);
disp(['Clock connection on Device 2 : ' clock2]);
disp('Clock connections must be physically connected using BNC cable');
TCS2.clock1=clock1;
TCS2.clock2=clock2;

%set the DAQ samplingrate to 1000 Hz
disp(['Setting sampling rate to ' num2str(samplingrate)]);
TCS2.NI.Rate=samplingrate;
TCS2.Rate=samplingrate;

%display information
disp(['Analog OUTPUTS on DEV', nbstr_NI1]); 
disp(['DEV',nbstr_NI1,'/AO1 : neutral temperature']);
disp(['DEV',nbstr_NI1,'/AO2 : active temperature zone 1']);
disp(['DEV',nbstr_NI1,'/AO3 : active temperature zone 2']);
disp(['DEV',nbstr_NI1,'/AO4 : active temperature zone 3']);
disp('');
disp(['Analog OUTPUTS on DEV', nbstr_NI2]);
disp(['DEV',nbstr_NI2,'/AO1 : active temperature zone 4']);
disp(['DEV',nbstr_NI2,'/AO2 : active temperature zone 5']);
disp(['DEV',nbstr_NI2,'/AO3 : neutral temperature trigger']);
disp(['DEV',nbstr_NI2,'/AO4 : active temperature trigger']);

%Add four analog output channels on 'Dev1'
addAnalogOutputChannel(TCS2.NI,dev1,0,'Voltage');
addAnalogOutputChannel(TCS2.NI,dev1,1,'Voltage');
addAnalogOutputChannel(TCS2.NI,dev1,2,'Voltage');
addAnalogOutputChannel(TCS2.NI,dev1,3,'Voltage');

%Add four analog output channels on 'Dev2'
addAnalogOutputChannel(TCS2.NI,dev2,0,'Voltage');
addAnalogOutputChannel(TCS2.NI,dev2,1,'Voltage');
addAnalogOutputChannel(TCS2.NI,dev2,2,'Voltage');
addAnalogOutputChannel(TCS2.NI,dev2,3,'Voltage');

%Add four analog input channels on 'Dev1'
disp(['Analog INPUTS on DEV',nbstr_NI1,' : 4']);
addAnalogInputChannel(TCS2.NI,dev1,0,'Voltage');
addAnalogInputChannel(TCS2.NI,dev1,1,'Voltage');
addAnalogInputChannel(TCS2.NI,dev1,2,'Voltage');
addAnalogInputChannel(TCS2.NI,dev1,3,'Voltage');

%Add four analog input channels on 'Dev2'
disp(['Analog INPUTS on DEV',nbstr_NI2,' : 4']);
addAnalogInputChannel(TCS2.NI,dev2,0,'Voltage');
addAnalogInputChannel(TCS2.NI,dev2,1,'Voltage');
addAnalogInputChannel(TCS2.NI,dev2,2,'Voltage');
addAnalogInputChannel(TCS2.NI,dev2,3,'Voltage');

disp('Analog OUTPUTS and INPUTS added.');

%synchronize the two devices

%add trigger connection
disp(['Adding trigger connection from ' trigger1 ' to ' trigger2]);
addTriggerConnection(TCS2.NI,trigger1,trigger2,'StartTrigger');

%add clock connection
disp(['Adding clock connection from ' clock1 ' to ' clock2]);
TCS2.NI.addClockConnection(clock1,clock2,'ScanClock');

%trigger voltage for neutral and active control
TCS2.trigger_voltage=3.3;

%voltage_temperature conversion
% 0deg C = 0V
%20deg C = 2V
%40deg C = 4V
%60deg C = 6V
%voltage=(temperature*TCS2.V_slope)+TCS2.V_intercept
TCS2.V0_temp = 0; % temp for V = 0
TCS2.V1_temp = 10;% temp for V = 1
TCS2.V_slope = 1/(TCS2.V1_temp-TCS2.V0_temp);% delta_volt/delta_temp
TCS2.V_intercept = 1-(TCS2.V_slope*TCS2.V1_temp); % V1-slope*V1_temp
TCS2.neutral_V_min = 2;
TCS2.neutral_V_max = 4;
TCS2.active_V_min = 0;
TCS2.active_V_max = 6;

disp('Initialized!');
end

