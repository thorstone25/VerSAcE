function S = vantageroot(), S = fileparts(which('VSX.m'));
% VANTAGEROOT - Root directory of the Vantage installation.
% 
% S = VANTAGEROOT returns the name of the directory where the Vantage software is
% installed.
% 
% VANTAGEROOT is used to produce platform dependent paths to the various MATLAB and
% toolbox directories.
% 
% Example:
% 
% % Save the setup struct 'vs' within the Vantage directory
% filename = fullfile(vantageroot, 'MatFiles/qups-vsx.mat');
% save(filename, '-struct', 'vs');
% 
% % Run VSX from the Vantage directory
% run VSX;
% 
% See also MATLABROOT