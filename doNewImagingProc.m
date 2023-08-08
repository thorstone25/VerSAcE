function doNewImagingProc(varargin)
% DONEWIMAGINGPROC - Activate or deactivate the NewImageProc function
%
% When the button is toggled, this function switches the global variable
% `toggle`'s state. Inputs don't matter.

global TOGGLE_newImagingProc; if isempty(TOGGLE_newImagingProc), TOGGLE_newImagingProc = false; else, TOGGLE_newImagingProc = ~TOGGLE_newImagingProc; end; 
end
