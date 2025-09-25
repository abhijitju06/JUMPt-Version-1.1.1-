%% Run_Main_File.m — JUMPt Setting-2 Robust Pipeline Version 1.1.1
% - Loads params from JUMPt.params / JUMPt_params.m / etc. (any extension)
% - Reads & bins input (robust to mixed-type Excel cells; trims bad rows)
% - Fits corrected half-lives (free-Lys participates in fit)
% - Optionally computes apparent T50 if apparent_T50_calculation == 1
% - Deterministic RNG
% - Author: Dr. Abhijit Dasgupta and Abhisek Bakshi

clear; clc;
fprintf('\n=== JUMPt Setting-2 Robust Pipeline ===\n');

%% Ensure script folder is on path (so binning / solvers are found)
scriptDir = fileparts(mfilename('fullpath'));
if ~isempty(scriptDir) && ~strcmpi(scriptDir, path)
    addpath(scriptDir);
end
addpath(pwd);  % also include current folder (for convenience)

%% Reproducible RNG (single seed used everywhere)
try
    s = RandStream('Threefry','Seed',123456);
    RandStream.setGlobalStream(s);
catch
    rng(123456,'twister');
end

%% -------- Locate & load params (works for .params/.txt/.m) --------
candidates = { ...
    fullfile(pwd, 'JUMPt.params'), ...
    fullfile(pwd, 'JUMPt_params.m'), ...
    fullfile(pwd, 'JUMPt.params.m'), ...
    fullfile(scriptDir, 'JUMPt.params'), ...
    fullfile(scriptDir, 'JUMPt_params.m'), ...
    fullfile(scriptDir, 'JUMPt.params.m') ...
};
paramsFile = '';
for i = 1:numel(candidates)
    if exist(candidates{i}, 'file'), paramsFile = candidates{i}; break; end
end
if isempty(paramsFile)
    d = dir(fullfile(pwd, 'JUMPt.params*'));      if ~isempty(d), paramsFile = fullfile(d(1).folder, d(1).name); end
end
if isempty(paramsFile)
    d = dir(fullfile(scriptDir, 'JUMPt.params*')); if ~isempty(d), paramsFile = fullfile(d(1).folder, d(1).name); end
end
if isempty(paramsFile)
    error('JUMPt.params not found.\nSearched:\n  %s\n  %s\n', pwd, scriptDir);
end
fprintf('Loading params from: %s\n', paramsFile);

% Read as text (robust to non-.m extensions / encodings) and eval into base
rawText = fileread(paramsFile);
if ~isempty(rawText) && rawText(1) == char(65279)  % strip UTF-8 BOM if present
    rawText = rawText(2:end);
end
evalin('base', rawText);

%% -------- Helper to pull values with defaults from params file --------
% (Local functions are defined at end of script: getv, iff)

%% -------- Assemble params struct (defaults) --------
params = struct();
params.input_file                   = getv('input_file', 'Small_Cell_Cere_Inputs.xlsx');
params.output_folder                = getv('output_folder', pwd);
params.bin_size                     = getv('bin_size', 30);                 % proteins per bin
params.number_of_timepoints         = getv('number_of_timepoints', []);     % auto-detect if empty
params.purity_of_SILAC_food         = getv('purity_of_SILAC_food', 99);     % %% heavy
params.apparent_T50_calculation     = getv('apparent_T50_calculation', 0);  % 1 => run Apt_HL_Calculation
% Very wide HL bounds  (extremely small/large HLs allowed)
params.hl_tmin_days                 = getv('hl_tmin_days', 1e-8);
params.hl_tmax_days                 = getv('hl_tmax_days', 1e8);
params.hl_A_min_days                = getv('hl_A_min_days', 1e-8);
params.seedHalfLifeDays             = getv('seedHalfLifeDays', 5);
params.alphaUB                      = getv('alpha_ub', getv('alphaUB', 0.2));
params.rng_seed                     = getv('rng_seed', 123456);
params.disable_physics_filter       = getv('disable_physics_filter', 0);
params.physics_filter_mode          = getv('physics_filter_mode', 'clip');  % 'clip' or 'drop'
params.physics_eps                  = getv('physics_eps', 1e-9);
params.min_observations_per_protein = getv('min_observations_per_protein', 2); % require ≥2 timepoints/protein
params.show_progress = getv('show_progress', 1);  % show waitbar & per-bin messages


% Ensure output folder exists
if ~exist(params.output_folder, 'dir')
    mkdir(params.output_folder);
end

%% -------- Load & bin data from Excel --------
if ~exist(params.input_file, 'file')
    error('Input file "%s" not found.', params.input_file);
end
fprintf('Reading and binning data from: %s\n', params.input_file);
[data, bins] = binning(params);  

% Validate detected time points vs. params.number_of_timepoints
T = numel(data.t);
if isempty(params.number_of_timepoints)
    params.number_of_timepoints = T;
elseif params.number_of_timepoints ~= T
    warning('number_of_timepoints in params (%d) != detected (%d). Using detected value %d.', ...
            params.number_of_timepoints, T, T);
    params.number_of_timepoints = T;
end

theta_FL = 1 - params.purity_of_SILAC_food/100;
fprintf('Parsed %d proteins across %d time points; SILAC food LIGHT fraction = %.4f\n', ...
        size(data.SILAC_data_allTimes,2), T, theta_FL);

%% -------- Fit corrected half-lives (Setting-2) --------
fprintf('\n>>> Fitting corrected half-lives...\n');
calc_half_lives_T(data, params);
fprintf('Completed exporting corrected half-lives.\n');

%% -------- Optional: Apparent T50 --------
if params.apparent_T50_calculation == 1
    if exist('Apt_HL_Calculation.m','file') == 2
        fprintf('\n>>> Computing apparent T50...\n');
        Apt_HL_Calculation(params);
        fprintf('Completed exporting apparent T50.\n');
    else
        warning('apparent_T50_calculation==1 but Apt_HL_Calculation.m not found on path. Skipping.');
    end
end

fprintf('\n*******  JUMPt program is complete *******\n');

%% ===== Local helper functions (script-local; require R2016b+) =====
function v = getv(name, def)
% Fetch a variable from base workspace (set by params file); return default if missing/empty.
    if evalin('base', sprintf('exist(''%s'',''var'')', name))
        tmp = evalin('base', name);
        if ~isempty(tmp)
            v = tmp;
        else
            v = def;
        end
    else
        v = def;
    end
end

function y = iff(c, a, b)
% Inline ternary
    if c, y = a; else, y = b; end
end
