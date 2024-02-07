function vs = update_vstruct()

global vs;

% load conf data
if isempty(vs), vs = load('qups-vsx.mat'); end

% update with base workspace variables
for f = string(fieldnames(vs))'
    if evalin('base', "exist('"+f+"', 'var')")
        vs.(f) = evalin('base', f);
    end
end





