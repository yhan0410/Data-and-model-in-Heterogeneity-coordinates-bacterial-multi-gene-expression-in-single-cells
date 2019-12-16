%--------------------------------------------------------------------------
% Author:  Yichao Han
% Created: 2nd Jul 2018
% Updated(1.0): 12th Nov 2018
% Updated(2.0): 30th Nov 2018
% Updated(3.0): 3rd Dec 2019
%--------------------------------------------------------------------------
% The model is used to numerically simulate single-cell level gene
% expression resource competition
%--------------------------------------------------------------------------
% Variables in transcription process
% DT_1, heterologous total promoters
% DT_2, endogenous total promoters
% DF_1, heterologous empty promoters
% DF_2, endogenous empty promoters
% DC_1, heterologous RNAP-bound promoters
% DC_2, endogenous RNAP-bound promoters
% RNAPT, total RNAP
% RNAPF, free RNAP
% M_1, heterologous mRNAs 
% M_2, endogenous mRNAs 
%--------------------------------------------------------------------------
% Variables in translation process
% MT_1, heterologous total mRNAs 
% MT_2, endogenous total mRNAs 
% MF_1, heterologous free mRNAs
% MF_2, endogenous free mRNAs
% MC_1, heterologous mRNA-Ribosome complex
% MC_2, endogenous mRNA-Ribosome complex
% RibF, free ribosome
% RibT, total ribosome
% P_1, heterologous proteins
% P_2, endogenous proteins 
%--------------------------------------------------------------------------
% mRNA dynamics
% Mass balance (i=1,2)
% RNAPT = RNAPF + DC_1 + DC_2;
% DT_i  = DC_1 + DF_i;
%
% ODEs
% RNAPF/dt= alpham_1 * DC_1 - alphap_1 * RNAPF * DF_1 + alpham_2 * DC_2 - alphap_2 * RNAPF * DF_2;
% DC_1/dt = alphap_1 * RNAPF * DF_1 - alpham_1 * DC_1;
% DC_2/dt = alphap_2 * RNAPF * DF_2 - alpham_2 * DC_2;
% M_1/dt  = alpham_1 * DC_1 - lambda * M_1;
% M_2/dt  = alpham_2 * DC_2 - lambda * M_2;
% 
% Protein dynamics
% Mass balance (i=1,2)
% RibT = RibF + n_1 * MC_1 + n_2 * MC_2;
% MT_i = MC_i + MF_i;
%
% ODEs
% RibF/dt = n_1 * (betam_1 * MC_1 - betap_1 * RibF * MF_1 ) + n_2 * (betam_2 * MC_2 - betap_2 * RibF * MF_2);
% MC_1/dt = betap_1 * RibF * MF_1 - betam_1 * MC_1;
% MC_2/dt = betap_2 * RibF * MF_2 - betam_2 * MC_2;
% P_1/dt  = betam_1 * MC_1 - lambda * P_1;
% P_2/dt  = betam_2 * MC_2 - lambda * P_2;
%% Set up the default parameters
clear;
default_parameters = parameter_table();
%% Transcriptional resource competition model
%	Analyze competition
parameters = default_parameters;
RibT_Range=7000:2000:15000;
MT_1_Range=0:1:250;
num_RibT_Range = length(RibT_Range);
num_MT_1_Range = length(MT_1_Range);
%   Initialize the outputs of steady state
P_1 = zeros(num_RibT_Range,num_MT_1_Range);
P_2 = zeros(num_RibT_Range,num_MT_1_Range);
RibF = zeros(num_RibT_Range,num_MT_1_Range);
for j = 1:num_RibT_Range
    parameters.RibT = RibT_Range(j);
    for i = 1:num_MT_1_Range
        parameters.MT_1 = MT_1_Range(i);
        [RibF(j,i),P_1(j,i),P_2(j,i)] = TL_solver(parameters);
    end
end
%--------------------------------------------------------------------------
%	Analyze correlation via linear approximation
parameters = default_parameters;
MT_1_cv_range = linspace(0, 0.5, 101);
RibT_cv_range = linspace(0, 0.5, 101);
MT_1_mean_range = linspace(0, 500, 101);
RibT_mean_range = linspace(0, 20000, 101);
P_PCC_a = zeros(101,101);
% parameters.MT_2_cv2 = 0;
% parameters.corr_M1_M2 = 0;
% for i = 1:length(MT_1_cv_range)
%     for j = 1:length(MT_1_mean_range)
for i = 1:length(RibT_cv_range)
    for j = 1:length(RibT_mean_range)
        %parameters.MT_1 = MT_1_mean_range(j);
        parameters.RibT = RibT_mean_range(j);
        MT_1_var = parameters.MT_1^2 * parameters.MT_1_cv2;
        %MT_1_var = parameters.MT_1^2 * MT_1_cv_range(i);
        %RibT_var = parameters.RibT^2 * parameters.RibT_cv2;%
        RibT_var = parameters.RibT^2 * RibT_cv_range(i);
        MT_2_var = parameters.MT_2^2 * parameters.MT_2_cv2;
        Cov_M1_M2 = sqrt(MT_1_var * MT_2_var) * parameters.corr_M1_M2;
        Cov_M1_RibT = sqrt(MT_1_var * RibT_var) * parameters.corr_M1_RibT;
        Cov_M2_RibT = sqrt(MT_2_var * RibT_var) * parameters.corr_M2_RibT;
        varianceMatrix = [MT_1_var, Cov_M1_M2, Cov_M1_RibT;
                          Cov_M1_M2 ,MT_2_var,Cov_M2_RibT ;
                          Cov_M1_RibT,Cov_M2_RibT,RibT_var];
        P_PCC_a(i,j) = correlation_approximation_protein(varianceMatrix, parameters);
    end
end
%% Transcriptional resource competition model
%	Analyze competition
parameters = default_parameters;
RNAPT_Range = 12000:-2000:4000;
DT_1_Range = 0:1:100;
num_DT_1_Range = length(DT_1_Range);
num_RNAPT_Range = length(RNAPT_Range);
%   Initialize the outputs of steady state
M_1 = zeros(num_RNAPT_Range,num_DT_1_Range);
M_2 = zeros(num_RNAPT_Range,num_DT_1_Range);
RNAPF = zeros(num_RNAPT_Range,num_DT_1_Range);
%   Calculate steady states under different initial conditions
for j = 1:length(RNAPT_Range)
    parameters.RNAPT = RNAPT_Range(j);
    for i = 1:length(DT_1_Range)
        parameters.DT_1 = DT_1_Range(i);
        [RNAPF(j,i), M_1(j,i), M_2(j,i)] = TX_solver(parameters);
    end
end
%%
%--------------------------------------------------------------------------
%	Analyze correlation via linear approximation
parameters = default_parameters;
DT_1_mean_range = linspace(0, 200, 101);
RNAPT_cv2_range = linspace(0, 0.5, 101);
DT_1_cv2_range = linspace(0, 0.5, 101);
M_PCC_a = zeros(101,101);
% for i = 1:length(DT_1_mean_range)
%     DT_1 = DT_1_mean_range(i);
%     RNAPT_var = RNAPT^2 * default_parameters.RNAPT_cv2;
for i = 1:length(RNAPT_cv2_range)
    RNAPT_var = parameters.RNAPT^2 * RNAPT_cv2_range(i);
    for j = 1:length(DT_1_cv2_range)
        DT_1_var = parameters.DT_1^2 * DT_1_cv2_range(j);
        [RNAPF,M_1,M_2,J]  = TX_solver(parameters);
        RNAPT_var = parameters.RNAPT^2 * parameters.RNAPT_cv2;
        DT_2_var = parameters.DT_2^2 * parameters.DT_2_cv2;
        Cov_D1_D2 = sqrt(DT_1_var * DT_2_var) * parameters.corr_M1_M2;
        Cov_D1_RNAPT = sqrt(DT_1_var * RNAPT_var) * parameters.corr_M1_RibT;
        Cov_D2_RNAPT = sqrt(DT_2_var * RNAPT_var) * parameters.corr_M2_RibT;
        varianceMatrix = [DT_1_var, Cov_D1_D2, Cov_D1_RNAPT;
                          Cov_D1_D2 ,DT_2_var, Cov_D2_RNAPT;
                          Cov_D1_RNAPT,Cov_D2_RNAPT,RNAPT_var];
        %M_PCC_a: approximated Pearson correlation coefficient between M1 and M2
        M_PCC_a(i,j) = correlation_approximation_mRNA(varianceMatrix, parameters);
    end
end
