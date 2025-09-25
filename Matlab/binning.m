function [data, bins] = binning(params)
% Robust loader for JUMPt Setting-2 sheets.
% Finds rows by label:
%   - "pulse time" row -> time points
%   - "free Lys" row   -> free-Lys (LIGHT) curve
% Proteins start at (free Lys row + 1).

    % ---------- RNG ----------
    if ~isfield(params,'rng_seed')||isempty(params.rng_seed), params.rng_seed=123456; end
    try, s=RandStream('Threefry','Seed',params.rng_seed); RandStream.setGlobalStream(s);
    catch, rng(params.rng_seed,'twister'); end

    % ---------- defaults ----------
    if ~isfield(params,'bin_size')||isempty(params.bin_size), params.bin_size=30; end
    if ~isfield(params,'physics_filter_mode')||isempty(params.physics_filter_mode), params.physics_filter_mode='clip'; end
    if ~isfield(params,'disable_physics_filter'), params.disable_physics_filter=0; end
    if ~isfield(params,'physics_eps')||isempty(params.physics_eps), params.physics_eps=1e-9; end
    if ~isfield(params,'purity_of_SILAC_food')||isempty(params.purity_of_SILAC_food), params.purity_of_SILAC_food=99; end
    if ~isfield(params,'min_observations_per_protein')||isempty(params.min_observations_per_protein)
        params.min_observations_per_protein = 2;
    end
    if ~isfield(params,'input_file') || ~exist(params.input_file,'file')
        error('binning:input','Input file "%s" not found.', string(getfield(params,'input_file','(unset)')));
    end

    % ---------- read Excel ----------
    raw = readcell(params.input_file);

    nRows = size(raw,1); nCols = size(raw,2);
    col1  = strings(nRows,1);
    for i=1:nRows
        try, col1(i) = string(raw{i,1}); catch, col1(i) = ""; end
    end
    col1lower = lower(strtrim(col1));

    % Find labeled rows (robust to wording/case/spaces)
    pulseRow = find(contains(col1lower, "pulse") & contains(col1lower, "time"), 1, 'first');
    freeRow  = find(contains(col1lower, "free")  & contains(col1lower, "lys"),  1, 'first');

    if isempty(pulseRow)
        error('binning:labels','Could not find a "pulse time" row (e.g., "pulse time (days)").');
    end
    if isempty(freeRow)
        error('binning:labels','Could not find a "free Lys" row (e.g., "free Lys (pSILAC ratio)").');
    end

    % ---------- detect time columns using those two rows ----------
    isNumPulse = false(1, nCols);
    isNumFree  = false(1, nCols);
    valsPulse  = nan(1, nCols);
    for j=1:nCols
        vP = toNumScalar(safeGet(raw, pulseRow, j));
        vF = toNumScalar(safeGet(raw, freeRow , j));
        valsPulse(j)  = vP;
        isNumPulse(j) = ~isnan(vP);
        isNumFree(j)  = ~isnan(vF);
    end
    timeMask = isNumPulse & isNumFree;
    % Never treat the first two columns (labels) as time columns
    timeMask(1:min(2,nCols)) = false;

    if ~any(timeMask)
        % fallback to params.number_of_timepoints if provided
        if ~isfield(params,'number_of_timepoints') || isempty(params.number_of_timepoints) || params.number_of_timepoints<2
            error('binning:times','Could not detect time columns from the labeled rows.');
        end
        T = params.number_of_timepoints;
        timeCols = 3:(2+T);
    else
        timeCols = find(timeMask);
    end

    % Build time vector (from the "pulse time" row)
    t_vec = valsPulse(timeCols).';              % T×1
    % Free-Lys (from the "free Lys" row)
    T = numel(timeCols);
    theta_AL = nan(T,1);
    for k=1:T, theta_AL(k) = toNumScalar(safeGet(raw, freeRow, timeCols(k))); end

    % ---------- proteins start at (freeRow + 1) ----------
    data_start_row = freeRow + 1;
    if data_start_row > nRows
        error('binning:rows','No protein rows found after the "free Lys" row.');
    end

    Ppept = raw(data_start_row:end, 1);
    Pgene = raw(data_start_row:end, 2);
    ProtY = raw(data_start_row:end, timeCols);     % N x T cells

    % Coerce to strings (keep empty as "")
    Ppept = cellstr(string(Ppept));
    Pgene = cellstr(string(Pgene));

    % Coerce protein matrix to double (NaNs allowed)
    Nraw = size(ProtY,1);
    Ymat = nan(Nraw, T);
    for i=1:Nraw
        for k=1:T
            Ymat(i,k) = toNumScalar(ProtY{i,k});
        end
    end

    % ---------- optional physics filter ----------
    if ~params.disable_physics_filter
        switch lower(params.physics_filter_mode)
            case 'clip'
                for k=1:T
                    if isnan(theta_AL(k)), continue; end
                    mask = ~isnan(Ymat(:,k));
                    bad  = Ymat(mask,k) < (theta_AL(k) - params.physics_eps);
                    if any(bad)
                        tmp = Ymat(mask,k);
                        tmp(bad) = theta_AL(k) - params.physics_eps;
                        Ymat(mask,k) = tmp;
                    end
                end
            case 'drop'
                drop_mask = false(Nraw,1);
                for i=1:Nraw
                    for k=1:T
                        if ~isnan(Ymat(i,k)) && ~isnan(theta_AL(k)) && Ymat(i,k) < (theta_AL(k) - params.physics_eps)
                            drop_mask(i) = true; break;
                        end
                    end
                end
                if any(drop_mask), Ymat(drop_mask,:) = NaN; end
        end
    end

    % ---------- trim non-protein / bad rows ----------
    n_obs  = sum(~isnan(Ymat), 2);
    minObs = max(1, params.min_observations_per_protein);
    isBlank = @(s) (strlength(strtrim(string(s))) == 0);
    keep = (n_obs >= minObs) & ~( isBlank(Ppept) & isBlank(Pgene) );

    Ppept = Ppept(keep);
    Pgene = Pgene(keep);
    Ymat  = Ymat(keep,:);
    N     = size(Ymat,1);

    % ---------- build outputs ----------
    ProtInfo = table(Ppept, Pgene, 'VariableNames', {'Peptide','Gene'});

    if isfinite(params.bin_size) && params.bin_size>0
        edges0 = 0:params.bin_size:N; if edges0(end)~=N, edges0=[edges0 N]; end
    else
        edges0 = [0 N];
    end

    theta_FL = 1 - params.purity_of_SILAC_food/100;

    data = struct();
    data.t                   = t_vec(:);
    data.LysRatio            = theta_AL(:);
    data.SILAC_data_allTimes = Ymat.';        % T x N
    data.ProtInfo_allTimes   = ProtInfo;
    data.protBins_allTimes   = edges0;
    data.theta_FL            = theta_FL;

    B = numel(edges0)-1;
    bins = repmat(struct('t',[], 'theta_AL',[], 'theta_PiL',[], 'thetaFL',[], ...
                         'peptide',[], 'gene',[]), 1, B);
    for b=1:B
        s0=edges0(b)+1; e0=edges0(b+1); idx=s0:e0;
        bins(b).t        = data.t;
        bins(b).theta_AL = data.LysRatio.';
        bins(b).theta_PiL= Ymat(idx,:);
        bins(b).thetaFL  = theta_FL;
        bins(b).peptide  = Ppept(idx);
        bins(b).gene     = Pgene(idx);
    end
end

% ---- helpers ----
function v = safeGet(C,i,j)
    if i<=size(C,1) && j<=size(C,2), v = C{i,j}; else, v = []; end
end

function d = toNumScalar(x)
% Single scalar double (NaN if not numeric); extracts first numeric substring from text.
    d = NaN; if isempty(x), return; end
    if iscell(x), if numel(x)>=1, d = toNumScalar(x{1}); end; return; end
    if isnumeric(x) || islogical(x), d = double(x(1)); return; end
    if isdatetime(x), d = double(days(x(1) - datetime(0,1,1))); return; end
    if isduration(x), d = double(days(x(1))); return; end
    if isstring(x) || ischar(x)
        xs = char(x);
        tok = regexp(xs,'[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?','match','once');
        if ~isempty(tok)
            y = str2double(tok); if ~isnan(y), d = y; return; end
        end
        y = str2double(xs); if ~isnan(y), d = y; end
    end
end
