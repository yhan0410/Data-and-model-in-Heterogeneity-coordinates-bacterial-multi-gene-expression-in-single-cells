function [RNAPF, M_1, M_2, J]=TX_solver(parameters)
%	TX_solver calculates output variables RNAPF, M_1, M_2, and Jacobian matrix at steady state
%	Input: kinetic parameters and steady state input variables RNAPT, DT_1, DT_2, and DT_3
    alphap_1 = parameters.alphap_1;
    alpham_1 = parameters.alpham_1;
    alphap_2 = parameters.alphap_2;
    alpham_2 = parameters.alpham_2;
    alphap_3 = parameters.alphap_3;
    alpham_3 = parameters.alpham_3;
    gamma_1 = parameters.gamma_1;
    gamma_2 = parameters.gamma_2;
    m_1 = parameters.m_1;
    m_2 = parameters.m_2;
    m_3 = parameters.m_3;    
    alpha_1 = alpham_1/alphap_1;
    alpha_2 = alpham_2/alphap_2;
    alpha_3 = alpham_3/alphap_3;
    RNAPT = parameters.RNAPT;
    DT_1 = parameters.DT_1;
    DT_2 = parameters.DT_2;
    DT_3 = parameters.DT_3;
%   Core functions
    TXfunc = @(RNAPF, RNAPT, DT_1, DT_2, DT_3, alpha_1, alpha_2, alpha_3, m_1, m_2, m_3)...
        RNAPF * (1 + m_1*DT_1/(alpha_1+RNAPF) + m_2*DT_2/(alpha_2+RNAPF) + m_3*DT_3/(alpha_3+RNAPF)) - RNAPT;
    func = @(RNAPF) TXfunc(RNAPF, RNAPT, DT_1, DT_2, DT_3, alpha_1, alpha_2, alpha_3, m_1, m_2, m_3);
%   Solve steady-state RNAPF first, and then M_1 and M_2
    RNAPF = fzero(func,[0, RNAPT]);
    M_1 = m_1 * alpham_1 / gamma_1 * DT_1 * RNAPF / (alpha_1 + RNAPF);
    M_2 = m_2 * alpham_2 / gamma_2 * DT_2 * RNAPF / (alpha_2 + RNAPF);
%     M_3 = m_3 * alpham_3 / gamma_3 * DT_3 * RNAPF / (alpha_3 + RNAPF);
%   Jacobian matrix
    dRNAPFdDT_1 = (-m_1 * RNAPF / (alpha_1 + RNAPF)) / ...
        (1 + m_1 * alpha_1 * DT_1 / (alpha_1 + RNAPF)^2 ...
           + m_2 * alpha_2 * DT_2 / (alpha_2 + RNAPF)^2 ...
           + m_3 * alpha_3 * DT_3 / (alpha_3 + RNAPF)^2);
    dRNAPFdDT_2 = (-m_2 * RNAPF / (alpha_2 + RNAPF)) / ...
        (1 + m_1 * alpha_1 * DT_1 / (alpha_1 + RNAPF)^2 ...
           + m_2 * alpha_2 * DT_2 / (alpha_2 + RNAPF)^2 ...
           + m_3 * alpha_3 * DT_3 / (alpha_3 + RNAPF)^2);
    dM_1dDT_1 = m_1 * alpham_1 / gamma_1 * ... 
        (alpha_1 * dRNAPFdDT_1 * DT_1 / (alpha_1 + RNAPF)^2 + RNAPF / (alpha_1 + RNAPF));
    dM_2dDT_1 = m_2 * alpham_2 / gamma_2 * ...
        alpha_2 * dRNAPFdDT_1 * DT_2 / (alpha_2 + RNAPF)^2;
    dM_1dDT_2 = m_1 * alpham_1 / gamma_1 * ... 
        alpha_1 * dRNAPFdDT_2 * DT_1 / (alpha_1 + RNAPF)^2;
    dM_2dDT_2 = m_2 * alpham_2 / gamma_2 * ...
        (alpha_2 * dRNAPFdDT_2 * DT_2 / (alpha_2 + RNAPF)^2 + RNAPF / (alpha_2 + RNAPF));
    
    dRNAPFdRNAPT = 1 / ...
        (1 + m_1 * alpha_1 * DT_1 / (alpha_1 + RNAPF)^2 ...
           + m_2 * alpha_2 * DT_2 / (alpha_2 + RNAPF)^2 ...
           + m_3 * alpha_3 * DT_3 / (alpha_3 + RNAPF)^2);
    dM_1dRNAPT = m_1 * gamma_1 / alpham_1 * ... 
        alpha_1 * dRNAPFdRNAPT * DT_1 / (alpha_1 + RNAPF)^2;
    dM_2dRNAPT = m_2 * gamma_2 / alpham_2 * ... 
        alpha_2 * dRNAPFdRNAPT * DT_2 / (alpha_2 + RNAPF)^2;
    J = [dM_1dDT_1, dM_1dDT_2, dM_1dRNAPT; dM_2dDT_1, dM_2dDT_2, dM_2dRNAPT];
end