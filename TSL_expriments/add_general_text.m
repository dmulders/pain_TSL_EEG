function add_general_text(w, str_to_disp, xpos, ypos, formatted_text,...
    center_txt, txt_col, backgrd_col, yPositionIsBaseline)
% Draw text on a screen with the Psychtoolbox.
if nargin<9
   yPositionIsBaseline = 0 ;  
end
if nargin<8
    backgrd_col = 0 ; 
end
if nargin<7
    txt_col = 255 ; 
end
if nargin<6
    center_txt = 0 ; 
end

if formatted_text
    if center_txt
        DrawFormattedText(w, str_to_disp, 'center', ypos, txt_col);
    else
        DrawFormattedText(w, str_to_disp, xpos, ypos, txt_col);
    end
else
    Screen('DrawText', w, str_to_disp, xpos, ypos, txt_col, backgrd_col, yPositionIsBaseline) ; 
end

end
