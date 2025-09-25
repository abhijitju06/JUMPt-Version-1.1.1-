function r = GO_lsqnonlin(param, tspan, Lys_P_ratio, SILAC_food_impurity)
% GO_lsqnonlin
% Residual function for Setting-2 global least-squares fit.
% We simulate the coupled ODE (free-Lys + proteins) and return a
% vector of residuals between observed and simulated ratios.
%
% --- Parameterization ----------------------------------------------------
% param = [HL_A, HL_1 .. HL_M, alpha]
%   HL_A     : half-life (days) of free-Lys pool (index 1)
%   HL_1..M  : half-lives (days) of proteins in the current bin
%   alpha    : coupling scalar that links protein kinetics to free-Lys
%
% --- Inputs --------------------------------------------------------------
% tspan              : column or row vector of time points (days)
% Lys_P_ratio        : T x (1+M) matrix of observed ratios
%                      [:,1]   -> free-Lys
%                      [:,2:]  -> protein ratios
% SILAC_food_impurity: scalar in [0,1], e.g., 0.01 for 99% purity
%
% --- Output --------------------------------------------------------------
% r : vectorized residuals used by lsqnonlin / MultiStart
%
% --- Notes ---------------------------------------------------------------
% 1) This version emphasizes the free-Lys trajectory using a weight
%    roughly equal to the number of proteins in the bin.
% 2) Determinism of initial starts should be handled in the caller
%    (e.g., rng(..), CustomStartPointSet, MultiStart('UseParallel','never')).

    % ---- Shapes & parsing ------------------------------------------------
    t = tspan(:);                        % ensure column
    [~, K] = size(Lys_P_ratio);          % K = 1 (free-Lys) + M (proteins)
    M = K - 1;

    % parameters -> rates
    HL    = param(1:K);
    alpha = param(end);

    % Guard against zero/negative HL from optimizer wander
    % (smallest positive normal double keeps behavior but avoids Inf)
    HL = max(HL, realmin('double'));

    gamma = log(2) ./ HL(:);             % Kx1 rates (free-Lys first)

    % ---- ODE simulate ----------------------------------------------------
    % State ordering y = [free-Lys ; protein_1 ; ... ; protein_M]
    y0    = ones(K,1);                   % initial condition
    p_ode = [gamma; alpha];              % pack parameters for PT_ODE_Ratio

    % Use default tolerances here; caller may set odeset if needed
    [~, Ysim] = ode15s(@(tt,yy) PT_ODE_Ratio(tt, yy, p_ode, SILAC_food_impurity), t, y0);

    % ---- Residual assembly w/ legacy-style weighting --------------------
    % Pipeline effectively privileged free-Lys. We emulate that by
    % giving free-Lys a weight ~= number of proteins in the bin.
    wF = M/5;              % try M/5 .. 2*M to fine-tune if ever needed
    wP = 1.0;

    % residual matrix (T x K)
    R = Lys_P_ratio - Ysim(:,1:K);

    % apply weights
    if K >= 1
        R(:,1)     = wF * R(:,1);        % free-Lys residuals
    end
    if K >= 2
        R(:,2:end) = wP * R(:,2:end);    % protein residuals
    end

    % vectorize and drop NaNs (missing values are ignored)
    r = R(~isnan(R));
end
