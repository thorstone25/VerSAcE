function vs = update_vstruct()

global vs;

% Properties to be added
allVSXProperties = {'Trans', 'PData', 'Media', 'Resource', 'Event', 'TW', 'TX', ...
                  'TPC', 'TGC', 'Receive', 'Recon', 'ReconInfo', 'Process', ...
                  'SeqControl', 'UI'};

% Load conf data if vs is empty
if isempty(vs)
    vs = load('qups-vsx.mat');
end

% Update with base workspace variables or make an empty struct
field_names = string([fieldnames(vs); allVSXProperties']);
for f = 1:length(field_names)
    if evalin('base', "exist('" + field_names(f) + "', 'var')")
        vs.(field_names(f)) = evalin('base', field_names(f));
    else
        vs.(field_names(f)) = struct(); % Create an empty struct if variable does not exist
    end
end
