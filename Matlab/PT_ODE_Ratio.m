function dydt = PT_ODE_Ratio(~, y, param, SILAC_food_impurity)
% param = [gamma_A; gamma_1..gamma_M; alpha]; y(1)=free-Lys, y(2..)=proteins
    y = y(:);
    y1 = y(1);
    yprot = y(2:end);                 % M×1

    p = param(:);
    gamma_A   = p(1);
    gamma_vec = p(2:end-1);           % M×1
    alpha     = p(end);

    dy1 = gamma_A*(SILAC_food_impurity - y1) + alpha * dot(gamma_vec, (yprot - y1));
    dyP = gamma_vec .* (y1 - yprot);

    dydt = [dy1; dyP];                % column
end
