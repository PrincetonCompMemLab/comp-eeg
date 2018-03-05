function [string, RT] = GetEchoStringVert(windowPtr, msg, x, y, textColor, bgColor, useKbCheck, varargin)
% string = GetEchoString(window, msg, x, y, [textColor], [bgColor], [useKbCheck=0], [deviceIndex], [untilTime=inf], [KbCheck args...])
%
%
% Edited version of GetEchoString that echoes the typed text below the
% prompt text, as opposed to next to it. Also displays and converts the
% string to uppercase.
%
% Get a string typed at the keyboard. Entry is terminated by
% <return> or <enter>.
%
% Typed characters are displayed in the window. The delete
% character is handled correctly. Useful for i/o in a Screen window.
%
% If the optional flag 'useKbCheck' is set to 1 then KbCheck is used - with
% potential optional additional 'KbCheck args...' for getting the string
% from the keyboard. Otherwise GetChar is used. 'useKbCheck' == 1 is
% restricted to standard alpha-numeric keys (characters, letters and a few
% special symbols). It can't handle all possible characters and doesn't
% work with non-US keyboard mappings. Its advantage is that it works
% reliably on configurations where GetChar may fail, e.g., on MS-Vista and
% Windows-7.
%
% See also: GetNumber, GetString, GetEchoNumber
%

% 2/4/97    dhb       Wrote GetEchoNumber.
% 2/5/97    dhb       Accept <enter> as well as <cr>.
%           dhb       Allow string return as well.
% 3/3/97    dhb       Updated for new DrawText.
% 3/15/97   dgp       Created GetEchoString based on dhb's GetEchoNumber.
% 3/20/97   dhb       Fixed bug in erase code, it wasn't updated for new
%                       initialization.
% 3/31/97   dhb       More fixes for same bug.
% 2/28/98   dgp       Use GetChar instead of obsolete GetKey. Use SWITCH and LENGTH.
% 3/27/98   dhb       Put an abs around char in switch.
% 12/26/08  yaosiang  Port GetEchoString from PTB-2 to PTB-3.
% 03/20/08  tsh       Added FlushEvents at the start and made bgColor and
%                     textcolor optional
% 10/22/10  mk        Optionally allow to use KbGetChar for keyboard input.

if nargin < 7
    useKbCheck = [];
end

if isempty(useKbCheck)
    useKbCheck = 0;
end

if nargin < 6
    bgColor = [];
end

if nargin < 5
    textColor = [];
end

if ~useKbCheck
    % Flush the keyboard buffer:
    FlushEvents;
end

% Write the message
Screen('DrawText', windowPtr, msg, x, y, textColor, bgColor);
Screen('Flip', windowPtr, 0, 1);

string = '';
RT = 0;
while true
    %Get RT for first typed character
    if useKbCheck && RT == 0
        [char, when] = GetKbChar(varargin{:});
        %Throwaway character so we don't get duplicates
        GetChar;
    else
        char = GetChar;
    end
    if exist('when', 'var')
        RT = when;
        clear when;
    end
    
    if isempty(char)
        string = '';
        return;
    end
    
    switch (abs(char))
        case {13, 3, 10}
            % ctrl-C, enter, or return
            break;
        case 8
            % backspace
            if ~isempty(string)
                string = string(1:length(string)-1);
            end
        otherwise
            string = [string, char]; %#ok<AGROW>
    end
%     fprintf('String = %s\n', string);
    output = sprintf([msg, '\n', upper(string)]);
    %     Screen('DrawText', windowPtr, output, x, y, textColor, bgColor);
    DrawFormattedText(windowPtr, output, x, y, textColor, bgColor);
    Screen('Flip', windowPtr);
%     if first
%         WaitSecs(0.01);
%     end
end
string = upper(string);