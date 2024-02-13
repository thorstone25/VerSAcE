function b = bfQUPS(RData)

% return a b mode image to VSX
persistent vs us chd0 prms;

global TOGGLE_bfQUPS;
global QUPS_BF_PARAMS; 

% on trigger
if isempty(TOGGLE_bfQUPS), TOGGLE_bfQUPS = false; end

% init
if isempty(us) || isempty(chd0)
    load('qups-conf', 'us', 'chd'); % load
    chd0 = chd;
    try chd0 = gpuArray(chd0); end %#ok<TRYNC>    
    
    % load params
    vs = update_vstruct(); % load configuration
    [prms, chd0] = QUPS_RT_update_params(QUPS_BF_PARAMS, vs, chd0); % load params into persistent memory
end

if ~TOGGLE_bfQUPS
    b = zeros(us.scan.size);
else
    
    tic;
    chd = QUPS_RT_load_data(RData, chd0, vs, prms);
    disp("Data loaded in " + toc() + " seconds."); % loading time
    
    % pre-process
    chd = hilbert(filter(chd, prms.D));
    
    % get image
    tic;
    b = DAS(us, chd);
    disp("Beamformed in " + toc() + " seconds."); % bf time
    b = double(gather(abs(b)));
end




