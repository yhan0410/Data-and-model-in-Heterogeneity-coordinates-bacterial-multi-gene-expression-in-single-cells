function [RibF, P_1, P_2, J] = TL_solver(parameters)
%	TL_solver calculates output variables RibF, P_1, P_2, and Jacobian matrix at steady state
%	Input: kinetic parameters and steady state input variables RibT, MT_1, MT_2
    betap_1 = parameters.betap_1; %Translation initiation rate for heterologous gene (default: 2e-3). 
    betam_1 = parameters.betam_1; %Translation enlongation rate for heterologous gene (default 20/1500). Assume 1500 AA and 20 AA/sec 
    betap_2 = parameters.betap_2; %Translation initiation rate for endogenous gene (default: 2e-4). 
    betam_2 = parameters.betam_2; %Translation enlongation rate for endogenous gene (default 20/360). Assume 360 AA per protein (average value in E. coli)
    beta_1 = betam_1 / betap_1;
    beta_2 = betam_2 / betap_2;
    n_1 = 30;
    n_2 = 2;
    lambda = parameters.lambda;% log(2)/5400
    RibT = parameters.RibT;
    MT_1 = parameters.MT_1;
    MT_2 = parameters.MT_2;
%   Core functions
    TLfunc = @(RibF, RibT, MT_1, MT_2, beta_1, beta_2, n_1, n_2) ...
        RibF * (1 + n_1 * MT_1 / (beta_1 + RibF) + n_2 * MT_2 /(beta_2 + RibF)) - RibT;        
    func = @(RibF) TLfunc(RibF, RibT, MT_1, MT_2, beta_1, beta_2, n_1, n_2);
%   Solve steady-state RibF first, and then P_1 and P_2
    RibF = fzero(func, [0,RibT]);
    P_1 = n_1 * betam_1 / lambda * MT_1 * RibF / (beta_1 + RibF);
    P_2 = n_2 * betam_2 / lambda * MT_2 * RibF / (beta_2 + RibF);
%   Jacobian matrix
    dRibFdMT_1 = (-n_1 * RibF / (beta_1 + RibF)) / ...
        (1 + n_1 * beta_1 * MT_1 / (beta_1 + RibF)^2 + n_2 * beta_2 * MT_2 / (beta_2 + RibF)^2);
    dP_1dMT_1 = n_1 * betam_1 / lambda * ... 
        (beta_1 * dRibFdMT_1 * MT_1 / (beta_1 + RibF)^2 + RibF / (beta_1 + RibF));
    dP_2dMT_1 = n_2 * betam_2 /  lambda * ...
        beta_2 * dRibFdMT_1 * MT_2 / (beta_2 + RibF)^2;
    dRibFdMT_2 = (-n_2 * RibF / (beta_2 + RibF)) / ...
        (1 + n_1 * beta_1 * MT_1 / (beta_1 + RibF)^2 + n_2 * beta_2 * MT_2 / (beta_2 + RibF)^2);
    dP_1dMT_2 = n_1 * betam_1 / lambda * ... 
        beta_1 * dRibFdMT_2 * MT_1 / (beta_1 + RibF)^2;
    dP_2dMT_2 = n_2 * betam_2 / lambda * ...
        (beta_2 * dRibFdMT_2 * MT_2 / (beta_2 + RibF)^2 + RibF / (beta_2 + RibF));
    dRibFdRibT = 1 / ...
        (1 + n_1 * beta_1 * MT_1 / (beta_1 + RibF)^2 + n_2 * beta_2 * MT_2 / (beta_2 + RibF)^2);
    dP_1dRibT = n_1 * betam_1 / lambda * ... 
        beta_1 * dRibFdRibT * MT_1 / (beta_1 + RibF)^2;
    dP_2dRibT = n_2 * betam_2  / lambda * ... 
        beta_2 * dRibFdRibT * MT_2 / (beta_2 + RibF)^2;
    J = [dP_1dMT_1, dP_1dMT_2, dP_1dRibT; dP_2dMT_1, dP_2dMT_2, dP_2dRibT];
end