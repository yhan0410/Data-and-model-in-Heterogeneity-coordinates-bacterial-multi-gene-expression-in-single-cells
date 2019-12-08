function P_PCC_a = correlation_approximation_protein(varianceMatrix, parameters)
%   Calculate correlation in proteins using linear approximation
    [~, ~, ~, J] = TL_solver(parameters);
    MT_1_var = varianceMatrix(1,1);
    MT_2_var = varianceMatrix(2,2);
    RibT_var = varianceMatrix(3,3);
    M1_M2_cov = varianceMatrix(1,2);
    M1_RibT_cov = varianceMatrix(1,3);
    M2_RibT_cov = varianceMatrix(2,3);
    P_1_var_a = J(1,1)^2 * MT_1_var + J(1,2)^2 * MT_2_var + J(1,3)^2 * RibT_var ...
        + 2 * J(1,1) * J(1,2) * M1_M2_cov + 2 * J(1,1) * J(1,3) * M1_RibT_cov + 2 * J(1,2) * J(1,3) * M2_RibT_cov;
    P_2_var_a = J(2,1)^2 * MT_1_var + J(2,2)^2 * MT_2_var + J(2,3)^2 * RibT_var ...
        + 2 * J(2,1) * J(2,2) * M1_M2_cov + 2 * J(2,1) * J(2,3) * M1_RibT_cov + 2 * J(2,2) * J(2,3) * M2_RibT_cov;
    P_cov_a = J(1,1) * J(2,1) * MT_1_var + J(1,2) * J(2,2) * MT_2_var + J(1,3) * J(2,3) * RibT_var ...
        + (J(1,1) * J(2,2) + J(2,1) * J(1,2)) * M1_M2_cov + (J(1,1) * J(2,3) + J(2,1) * J(1,3)) * M1_RibT_cov + (J(2,3) * J(1,2) + J(1,3) * J(2,2)) * M2_RibT_cov;
    P_PCC_a = P_cov_a ./ sqrt(P_1_var_a .* P_2_var_a);
end