function parameters = parameter_table()
    parameters = struct;
%   Parameters used in translational resource competition model
    parameters.RibT = 10000; 
    parameters.MT_1 = 300;
    parameters.MT_2 = 4000;
    parameters.betap_1 = 2e-3; %Translation initiation rate for heterologous gene (default: 2e-3). 
    parameters.betam_1 = 20/1500;%Translation enlongation rate for heterologous gene (default 20/1500). Assume 1500 AA and 20 AA/sec 
    parameters.betap_2 = 2e-4; %Translation initiation rate for endogenous gene (default: 2e-4). 
    parameters.betam_2 = 20/360;%Translation enlongation rate for endogenous gene (default 20/360). Assume 360 AA per protein (average value in E. coli)
    parameters.beta_1 = parameters.betam_1 / parameters.betap_1;
    parameters.beta_2 = parameters.betam_2 / parameters.betap_2;
    parameters.n_1 = 30;
    parameters.n_2 = 2;
    parameters.lambda = log(2)/5400;
%   Parameters used in transcriptional resource competition model 
    parameters.RNAPT = 8000;
    parameters.DT_1 = 50;
    parameters.DT_2 = 700;
    parameters.DT_3 = 100;    
    parameters.alphap_1 = 1/10;
    parameters.alpham_1 = 60/5000;
    parameters.alphap_2 = 1/10;
    parameters.alpham_2 = 60/1000;
    parameters.alphap_3 = 1/5;
    parameters.alpham_3 = 90/5000;
    parameters.gamma_1 = 1/120;
    parameters.gamma_2 = 1/60;
    parameters.gamma_3 = log(2)/5400;
    parameters.m_1 = 8;
    parameters.m_2 = 1.6;
    parameters.m_3 = 11;
%   cv2: squared coefficient of variation
    parameters.RNAPT_cv2 = 0.1;
    parameters.DT_1_cv2 = 0.1;
    parameters.DT_2_cv2 = 0.1;
    parameters.RibT_cv2 = 0.1;
    parameters.MT_1_cv2 = 0.1;
    parameters.MT_2_cv2 = 0.1;
%   corr_A_B: Pearson correlation coefficient between A and B
    parameters.corr_M1_M2 = 0.1;
    parameters.corr_M1_RibT = 0;
    parameters.corr_M2_RibT = 0;
    parameters.corr_D1_D2 = 0.1;
    parameters.corr_D1_RNAPT = 0;
    parameters.corr_D2_RNAPT = 0;
end

