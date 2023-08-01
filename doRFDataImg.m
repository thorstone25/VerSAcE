function doRFDataImg(varargin)
% DORFDATAPROC - Activate or deactivate the RFDataProc function
%
% When the button is toggled, this function switches the global variable
% `toggle`'s state. Inputs don't matter.

global TOGGLE_RFDataImg; if isempty(TOGGLE_RFDataImg), TOGGLE_RFDataImg = false; else, TOGGLE_RFDataImg = ~TOGGLE_RFDataImg; end; 
end
