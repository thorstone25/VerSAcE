function doImagingProc(varargin)
% DOIMAGINGPROC - Activate or deactivate the NewImageProc function
%
% When the button is toggled, this function switches the global variable
% `toggle`'s state. Inputs don't matter.

global TOGGLE_imagingProc; if isempty(TOGGLE_imagingProc), TOGGLE_imagingProc = false; else, TOGGLE_imagingProc = ~TOGGLE_imagingProc; end; 
end
