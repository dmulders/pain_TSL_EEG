function [rejected_indices] = get_rejected_epochs(subj_idx, cond_idx)
% Gives the epoch indices to reject, for the subject subj_name and the
% condition cond_idx.
% TSL-EEG.

rejected_indices = [] ; 

switch subj_idx
    case 1
        % ===*=== Subject 1
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
        
    case 2
        % ===*=== Subject 2
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end        
        
    case 3
        % ===*=== Subject 3
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                rejected_indices = [38] ; % FC3
            case 7
                
            case 8
                
            case 9
                rejected_indices = [20, 21, 27] ; % 20, 21 (TP8), 27 (FC6)
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
        
    case 4
        % ===*=== Subject 4
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                rejected_indices = [2] ; % P6, all chans
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    
    case 5
        % ===*=== Subject 5
        switch cond_idx
            case 1
                rejected_indices = [18] ; % after 1 second, so could be kept...
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                rejected_indices = [1,2] ; 
                % 2 first stim NOT delivered (forgot to place the thermode, though it was the first "enter")
            case 8
                
            case 9
                
            case 10
                rejected_indices = [97] ; 
                % 97 (ou close to?): artifacts (movements)
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
        
    case 6
        % ===*=== Subject 6
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                rejected_indices = [2] ; 
            case 4
                
            case 5
                
            case 6
                rejected_indices = [97] ; 
            case 7
                
            case 8
                
            case 9
                rejected_indices = [69] ; 
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
        
    case 7
        % ===*=== Subject 7
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
        
    case 8
        % ===*=== Subject 8
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                rejected_indices = [4,5] ; 
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
        
    case 9
        % ===*=== Subject 9
        switch cond_idx
            case 1
                rejected_indices = [76] ; 
            case 2
                rejected_indices = [24, 70, 82] ; 
            case 3
                rejected_indices = [77, 89] ; 
            case 4
                rejected_indices = [69] ; 
            case 5
                
            case 6
                
            case 7
                %rejected_indices = [8] ;  % only PZ
            case 8
                
            case 9
                
            case 10
                %rejected_indices = [98] ; % Fp chan, >1sec
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 10
        % ===*=== Subject 10
        switch cond_idx
            %%%%%%%%%%%%% ATTENTION: LF oscillations... due to ANT
            %%%%%%%%%%%%% cables?!? cap a bit tight... because no L ok, M
            % --> some LF IC rejected...
            % a lot of alpha rhythm along almost all blocs, eg on FCZ...
            case 1
                
            case 2
                rejected_indices = [99, 100] ; % alpha ++
            case 3
                rejected_indices = [100] ; 
            case 4
                
            case 5
                
            case 6
                rejected_indices = [78] ; 
            case 7
                
            case 8
                rejected_indices = [73, 74] ; % 7, 11, 12, 13 on AF8...
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
        
    case 11
        % ===*=== Subject 11
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 12
        % ===*=== Subject 12
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                rejected_indices = [18] ; 
            case 8
                rejected_indices = [13, 14, 36] ; 
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 13
        % ===*=== Subject 13
        switch cond_idx
            case 1
                rejected_indices = [7] ; 
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                % Comm. EOG elec felt at the end of cond6 (2 last epochs)
                % --> new elec from cond7
                
            case 7
                rejected_indices = [36] ; 
            case 8
                
            case 9
                rejected_indices = [90] ; 
            case 10
                rejected_indices = [100] ; 
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 14
        % ===*=== Subject 14
        switch cond_idx
            case 1
                rejected_indices = [57] ; % (<-0.5sec)
            case 2
                
            case 3
                
            case 4
                
            case 5
                rejected_indices = [51] ; % (C6)
            case 6
                
            case 7
                
            case 8
                
            case 9
                rejected_indices = [66, 92] ; % C6
            case 10
                rejected_indices = [11, 76] ; % (C6, <-0.5sec)
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 15
        % ===*=== Subject 15
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 16
        % ===*=== Subject 16
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                rejected_indices = [70] ; % (M1++!!!)
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 17
        % ===*=== Subject 17
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                % Comm. gsm rang around ep ~78-80 of cond7.
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 18
        % ===*=== Subject 18
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                %%%%%%%%%% ATTENTION: cond9: 1 elec crappy --> try to
                %%%%%%%%%% reref?! (all chan are affected)
                % --> SEEMS OK after ICA!! (artifacts from M1 removed)
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 19
        % ===*=== Subject 19
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 20
        % ===*=== Subject 20
        switch cond_idx
            %%%%% ATTENTION: M1 crappy
            %%%%% LF artifacts (train, cond1, ...) --> reref?!?
            % --> seems ok...
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 21
        % ===*=== Subject 21
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 22
        % ===*=== Subject 22
        % A LOT OF alpha rhythm!!! --> new ICA with alpha removed (4 ICs)
        % icfilt-STRONG
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                rejected_indices = [25] ; % (C5-TP7-...)
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 23
        % ===*=== Subject 23
        switch cond_idx
            case 1
                
            case 2
                rejected_indices = [2] ; % (T7)
            case 3
                rejected_indices = [6, 100] ; % (6: Fp1, 100: O2)
            case 4
                
            case 5
                %rejected_indices = [5] ; % (FC5, >1sec)
            case 6
                
            case 7
                rejected_indices = [3] ; % (FT7)
            case 8
                
            case 9
                
            case 10
                rejected_indices = [79] ; % (P6)
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 24
        % ===*=== Subject 24
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                % Artifacts ep ~89, 95, ...
                % Maybe reref M1M2??
                % Attention M1 sometimes crappy (lifted off) -> IC1 removed
                rejected_indices = [7] ; % (all chans)
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 25
        % ===*=== Subject 25
        switch cond_idx
            case 1
                %rejected_indices = [9] ; % (C5-FC5, early)
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                rejected_indices = [100] ; % (FC5)
            case 8
                
            case 9
                
            case 10
                rejected_indices = [11] ; % (all chans)
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 26
        % ===*=== Subject 26
        switch cond_idx
            %%%% LF artifacts during first blocs (maybe reref?!) --> stair,
            %%%% train, ...
            case 1
                rejected_indices = [95] ; % (F5)
            case 2
                
            case 3
                rejected_indices = [100] ; % (F5)
            case 4
                
            case 5
                
            case 6
                rejected_indices = [50, 89] ; % (50: T7; 89: F7)
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 27
        % ===*=== Subject 27
        switch cond_idx
            %%%%% LF artifacts, eg during cond3!!
            %%%%% during first blocs: eye blink & mvt after each warm
            %%%%% stimulus!
            %%%%% Closes his eyes sometimes
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                rejected_indices = [16] ; % (all chans)
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 28
        % ===*=== Subject 28
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 29
        % ===*=== Subject 29
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 30
        % ===*=== Subject 30
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                rejected_indices = [13,96] ; % (FC4)
            case 9
                
            case 10
                rejected_indices = [69] ; % (AF4)
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 31
        % ===*=== Subject 31
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                %%%% Artifacts (around ep 87?!) --> mvts
                rejected_indices = [78] ; % (all chans)
            case 4
                rejected_indices = [16,72] ; % F6
            case 5
                
            case 6
                rejected_indices = [1] ; % TP8, FCZ, C2, ...
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 32
        % ===*=== Subject 32
        switch cond_idx
            case 1
                rejected_indices = [17] ; % F5, AF8, ...
            case 2
                rejected_indices = [30, 95] ; % 30: all; 95: C6
            case 3
                
            case 4
                rejected_indices = [16] ; % all chans
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                rejected_indices = [43, 85] ; % FC6
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 33
        % ===*=== Subject 33
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 34
        % ===*=== Subject 34
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                rejected_indices = [16] ; % all chans
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 35
        % ===*=== Subject 35
        switch cond_idx
            case 1
                rejected_indices = [20] ; % FC5
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    case 36
        % ===*=== Subject 36
        switch cond_idx
            case 1
                
            case 2
                
            case 3
                
            case 4
                
            case 5
                
            case 6
                
            case 7
                
            case 8
                
            case 9
                
            case 10
                
            otherwise
                warning(['The condition ', cond_idx, ' does not exist'])
        end
    otherwise
        warning(['The subject ', subj_idx, ' has not been found'])
end


end

