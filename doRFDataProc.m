function doRFDataStore(varargin)
% DORFDATAPROC - Activate or deactivate the RFDataProc function
%
% When the button is toggled, this function switches the global variable
% `toggle`'s state. Inputs don't matter.

global TOGGLE_RFDataProc; if isempty(TOGGLE_RFDataProc), TOGGLE_RFDataProc = false; else, TOGGLE_RFDataProc = ~TOGGLE_RFDataProc; end; 
end
