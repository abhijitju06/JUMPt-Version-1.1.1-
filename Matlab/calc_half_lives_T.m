function [] = calc_half_lives_T(data, params)
% Setting-2 solver (robust & deterministic)
% - 20 random starts (rng default)
% - alpha >= 0 (optional UB from params.alphaUB)
% - choose minimum TOTAL SSE solution  (MultiStart path)
% - manual path prefers best free-Lys SSE, then total SSE
% - positive HL CIs via log-gamma space
% - progress bar 

    % ---------- options / UI ----------
    showProgress = true;
    if isfield(params,'show_progress') && ~isempty(params.show_progress)
        showProgress = logical(params.show_progress);
    end

    baseOpts = optimoptions(@lsqnonlin, ...
        'Display','off', 'FiniteDifferenceType','central', ...
        'MaxFunctionEvaluations',1e5, ...
        'StepTolerance',1e-12, 'FunctionTolerance',1e-10);

    % ---------- output file ----------
    [~, base, ~] = fileparts(params.input_file);
    outfile = fullfile(params.output_folder, ...
        sprintf('results_Corrected_T50_%s.xlsx', base));

    SILAC_food_impurity = (100 - params.purity_of_SILAC_food) / 100;

    % ---------- bounds (lower) ----------
    if isfield(params,'hl_A_min_days') && ~isempty(params.hl_A_min_days)
        hl_min_A = max(eps, params.hl_A_min_days);
    else
        hl_min_A = 1e-8;
    end
    if isfield(params,'hl_tmin_days') && ~isempty(params.hl_tmin_days)
        hl_min_P = max(eps, params.hl_tmin_days);
    else
        hl_min_P = 1e-8;
    end
    a_min = 0;   % alpha >= 0

    % ---------- bin structure ----------
    edges0  = data.protBins_allTimes;
    B       = numel(edges0) - 1;
    all_rows = table();

    % headers for time-point columns
    T        = numel(data.t);
    tpNames  = arrayfun(@(k) sprintf('time_point_%d',k), 1:T, 'uni', false);

    % ---------- progress bar  ----------
    h = [];
    if showProgress && usejava('desktop')
        try
            h = waitbar(0, 'Fitting corrected half-lives (bin 0/0)...', ...
                        'Name', 'JUMPt', ...
                        'CreateCancelBtn','setappdata(gcbf,''canceling'',1)');
            setappdata(h,'canceling',0);
        catch
            h = [];
        end
    end
    wb_cleanup = onCleanup(@() safe_close_waitbar(h));  

    % ---------- main loop over bins ----------
    rng('default');               % deterministic across machines/versions
    haveMS = exist('MultiStart','class') == 8;

    for b = 1:B
        if ~isempty(h) && isvalid(h)
            if getappdata(h,'canceling'), error('User cancelled.'); end
            waitbar((b-1)/B, h, sprintf('Fitting corrected half-lives (bin %d/%d)...', b, B));
            drawnow limitrate
        end

        s0  = edges0(b)   + 1;
        e0  = edges0(b+1);
        idx = s0:e0;

        % ----- build per-bin observation matrix: Lys + proteins -----
        Yi   = data.SILAC_data_allTimes(:, idx);     % T×M (proteins)
        K    = size(Yi, 2) + 1;                       % 1+M (free-Lys + proteins)
        LysP = [data.LysRatio(:), Yi];                % T×K

        % ----- initial starts (positive) -----
        nStarts = 20;
        X0 = rand(K+1, nStarts);  % [HL_A, HL_1..HL_M, alpha]  all >0

        % ----- bounds -----
        LB = [max(eps,hl_min_A); repmat(max(eps,hl_min_P), K-1, 1); a_min];

        % per-bin upper bound for alpha (optional)
        if isfield(params,'alphaUB') && ~isempty(params.alphaUB)
            UBbin = [Inf(1, K), params.alphaUB];    % K HLs + alpha
        else
            UBbin = [];                              % no UB on alpha
        end

        % ----- objective -----
        obj = @(p) GO_lsqnonlin(p, data.t, LysP, SILAC_food_impurity);

        % ----- solve (MultiStart or manual) -----
        if haveMS
            if showProgress
                fprintf('  Bin %d/%d: %d proteins — MultiStart(%d starts)\n', ...
                        b, B, numel(idx), nStarts);
            end
            ms   = MultiStart('UseParallel','never','Display','off');
            cust = CustomStartPointSet(X0.');   % rows = start points
            prob = createOptimProblem('lsqnonlin', 'x0', X0(:,1), ...
                    'objective', obj, 'lb', LB, 'ub', UBbin, 'options', baseOpts);
            [bestX, ~] = run(ms, prob, cust);

            % refine once to get J for CI
            [bestX,~,bestRes,~,~,~,bestJ] = lsqnonlin(obj, bestX, LB, UBbin, baseOpts);

        else
            if showProgress
                fprintf('  Bin %d/%d: %d proteins — %d starts\n', ...
                        b, B, numel(idx), nStarts);
            end

            % ---- prefer best free-Lys SSE, tie-break by total SSE
            bestX=[]; bestRes=[]; bestJ=[]; bestSSE_T = inf; bestSSE_F = inf;
            for s = 1:nStarts
                if ~isempty(h) && isvalid(h)
                    if getappdata(h,'canceling'), error('User cancelled.'); end
                    waitbar((b-1 + min(s/nStarts,0.99))/B, h, ...
                        sprintf('Fitting bin %d/%d (start %d/%d)...', b, B, s, nStarts));
                    drawnow limitrate
                end

                [x_s,~,res_s,~,~,~,J_s] = lsqnonlin(obj, X0(:,s), LB, UBbin, baseOpts);

                % total SSE uses exactly the residual vector that GO_lsqnonlin returns
                sseT = sum(res_s.^2, 'omitnan');

                % free-Lys SSE: if no NaNs, residuals form T×K and col 1 is free-Lys
                Tloc = size(LysP,1);
                if numel(res_s) == Tloc * K
                    resM = reshape(res_s, Tloc, K);
                    sseF = nansum(resM(:,1).^2);
                else
                    % fallback: approximate by simulating and measuring free-Lys residuals
                    sseF = free_sse_from_sim(x_s, data.t, LysP(:,1), SILAC_food_impurity);
                end

                if (sseF < bestSSE_F) || (sseF == bestSSE_F && sseT < bestSSE_T)
                    bestSSE_F = sseF; bestSSE_T = sseT;
                    bestX = x_s; bestRes = res_s; bestJ = J_s;
                end
            end
        end

        % ----- CIs in LOG-GAMMA space (guarantees positive HL bounds) -----
        ln2       = log(2);
        gammaBest = ln2 ./ bestX(1:K);          % gamma = ln2 / HL
        thetaBest = log(gammaBest);             % theta = log(gamma) > -Inf

        if ~isempty(bestJ)
            % bestJ columns 1..K are dr/dHL_i. Convert to dr/dtheta_i:
            dHL_dgamma = -(ln2) ./ (gammaBest.^2);   % K×1
            Jtheta = bestJ(:,1:K);
            for i = 1:K
                Jgamma_col  = bestJ(:,i) * dHL_dgamma(i);   % dr/dgamma_i
                Jtheta(:,i) = Jgamma_col * gammaBest(i);    % dr/dtheta_i
            end
            % Build Jacobian for beta = [theta ; alpha]
            Jbeta = [Jtheta, bestJ(:,K+1)];

            beta = [thetaBest; bestX(end)];
            try
                Ci_beta = nlparci(beta, bestRes, 'jacobian', Jbeta, 'alpha', 0.05);
            catch
                Ci_beta = nan(numel(beta), 2);
            end

            % Back-transform: theta -> gamma -> HL
            thetaL = Ci_beta(1:K,1);  thetaU = Ci_beta(1:K,2);
            gL = exp(thetaL);         gU = exp(thetaU);      % strictly positive
            HL_low  = ln2 ./ gU;                           % lower HL bound
            HL_high = ln2 ./ gL;                           % upper HL bound

            CiHL = nan(numel(bestX),2);
            CiHL(1:K,1) = HL_low;
            CiHL(1:K,2) = HL_high;
        else
            CiHL = nan(K+1,2);
        end

        % ----- residual SSE per protein (for reporting) -----
        y0    = ones(K,1);
        p_ode = [gammaBest; bestX(end)];
        [~, Ysim] = ode15s(@(tt,yy) PT_ODE_Ratio(tt,yy,p_ode,SILAC_food_impurity), data.t, y0);
        diffP   = LysP(:,2:end) - Ysim(:,2:end);
        protSSE = nansum(diffP.^2, 1).';

        % ----- assemble output table for this bin -----
        HL_i      = bestX(2:K);
        CI_low_i  = CiHL(2:K,1);
        CI_high_i = CiHL(2:K,2);

        % >>> HL display / extreme flag (NEW) <<<
        HL_hi_thr = 100;    % default thresholds
        HL_lo_thr = 0.5;
        if isfield(params,'HL_display_upper') && ~isempty(params.HL_display_upper)
            HL_hi_thr = params.HL_display_upper;
        end
        if isfield(params,'HL_display_lower') && ~isempty(params.HL_display_lower)
            HL_lo_thr = params.HL_display_lower;
        end
        HL_disp = strings(numel(HL_i),1);
        HL_extreme = false(numel(HL_i),1);
        for ii = 1:numel(HL_i)
            if HL_i(ii) > HL_hi_thr
                HL_disp(ii) = sprintf('>%g', HL_hi_thr);
                HL_extreme(ii) = true;
            elseif HL_i(ii) < HL_lo_thr
                HL_disp(ii) = sprintf('<%g', HL_lo_thr);
                HL_extreme(ii) = true;
            else
                HL_disp(ii) = num2str(HL_i(ii), '%.6g'); % nicely formatted numeric
            end
        end
        % <<< end HL display block >>>

        tpVals = num2cell(Yi.', 2);   % one cell per protein, value is 1×T
        T_bin  = table( ...
            data.ProtInfo_allTimes.Peptide(idx), ...
            data.ProtInfo_allTimes.Gene(idx), ...
            'VariableNames', {'Proteins/Peptides','Gene/Protein-Gene'} ...
        );
        for k = 1:T
            T_bin.(tpNames{k}) = cellfun(@(row) row(k), tpVals);
        end
        T_bin.('HalfLife_in_days')        = HL_i(:);                 % numeric
        T_bin.('HalfLife_display')        = cellstr(HL_disp);        % ">100" / "<0.5" / number
        T_bin.('HL_is_extreme')           = HL_extreme(:);           % logical flag
        T_bin.('HL_CI_low')               = CI_low_i(:);
        T_bin.('HL_CI_high')              = CI_high_i(:);
        T_bin.('Confidence_Interval')     = CI_high_i(:) - CI_low_i(:);
        T_bin.('residual_error')          = protSSE;

        all_rows = [all_rows; T_bin]; %#ok<AGROW>

        if showProgress
            fprintf('  Bin %d/%d done: %d proteins.\n', b, B, numel(idx));
        end
    end

    % progress bar closes via onCleanup; calling explicitly is harmless
    safe_close_waitbar(h);

    % ---------- write outputs ----------
    if showProgress && usejava('desktop')
        try
            h2 = waitbar(1, 'Writing results...', 'Name', 'JUMPt'); 
        catch
        end
    end
    try
        % params sheet may fail if params has non-scalar fields; ignore
        writetable(struct2table(params), outfile, 'Sheet', 'parameter_file');
    catch
    end
    writetable(all_rows, outfile, 'Sheet', 'results');
    fprintf('Wrote %d proteins to "%s" (sheet: results).\n', height(all_rows), outfile);
end

% ---------- helpers ----------
function sseF = free_sse_from_sim(param, tspan, Lys_obs, impurity)
% Fallback: compute free-Lys SSE directly (unweighted) if needed.
    t = tspan(:);
    % infer K from parameter length: [HL_A, HL_1..HL_M, alpha]
    K = numel(param) - 1;
    HL    = max(param(1:K), realmin('double'));
    alpha = param(end);
    gamma = log(2) ./ HL(:);
    y0    = ones(K,1);
    p_ode = [gamma; alpha];
    [~, Y] = ode15s(@(tt,yy) PT_ODE_Ratio(tt,yy,p_ode,impurity), t, y0);
    rF = Lys_obs(:) - Y(:,1);
    sseF = nansum(rF.^2);
end

function safe_close_waitbar(h)
    try
        if ~isempty(h) && isvalid(h), delete(h); end
    catch
    end
end
