function M_PCC_a = correlation_approximation_mRNA(varianceMatrix, parameters)
%   Calculate correlation in mRNAs using linear approximation
    [~, ~, ~, J] = TX_solver(parameters);
    DT_1_var = varianceMatrix(1,1);
    DT_2_var = varianceMatrix(2,2);
    RNAPT_var = varianceMatrix(3,3);
    D1_D2_cov = varianceMatrix(1,2);
    D1_RNAPT_cov = varianceMatrix(1,3);
    D2_RNAPT_cov = varianceMatrix(2,3);
    M_1_var_a = J(1,1)^2 * DT_1_var + J(1,2)^2 * DT_2_var + J(1,3)^2 * RNAPT_var ...
        + 2 * J(1,1) * J(1,2) * D1_D2_cov + 2 * J(1,1) * J(1,3) * D1_RNAPT_cov + 2 * J(1,2) * J(1,3) * D2_RNAPT_cov;
    M_2_var_a = J(2,1)^2 * DT_1_var + J(2,2)^2 * DT_2_var + J(2,3)^2 * RNAPT_var ...
        + 2 * J(2,1) * J(2,2) * D1_D2_cov + 2 * J(2,1) * J(2,3) * D1_RNAPT_cov + 2 * J(2,2) * J(2,3) * D2_RNAPT_cov;
    M_cov_a = J(1,1) * J(2,1) * DT_1_var + J(1,2) * J(2,2) * DT_2_var + J(1,3) * J(2,3) * RNAPT_var ...
        + (J(1,1) * J(2,2) + J(2,1) * J(1,2)) * D1_D2_cov + (J(1,1) * J(2,3) + J(2,1) * J(1,3)) * D1_RNAPT_cov + (J(2,3) * J(1,2) + J(1,3) * J(2,2)) * D2_RNAPT_cov;
    M_PCC_a = M_cov_a ./ sqrt(M_1_var_a .* M_2_var_a);
end