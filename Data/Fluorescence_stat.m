function [PCC, Num_Cell, Mean, Std, CV2, Summary] = Fluorescence_stat(data)
    [PCC,pVal] = corr(data(:,1), data(:,2));
    PCC_ci = bootci(1000, @corr, data(:,1), data(:,2));
    Num_Cell = length(data(:,1));
    Mean = mean(data);
    Mean_ci = bootci(1000, @mean, data);
    Std = std(data);
    Std_ci = bootci(1000, @std, data);
    CV2 = (Std ./ Mean).^2;
    Noise_cal = @(x) (std(x) / mean(x))^2;
    CV2_ci(:,1) = bootci(1000, Noise_cal, data(:,1));
    CV2_ci(:,2) = bootci(1000, Noise_cal, data(:,2));
    %Summary = [Mean, Std,CV2, PCC];
    Summary = [Num_Cell, Mean(1), Mean_ci(:,1)', Mean(2), Mean_ci(:,2)', ...
        Std(1), Std_ci(:,1)', Std(2), Std_ci(:,2)', ...
        CV2(1), CV2_ci(:,1)', CV2(2), CV2_ci(:,2)', ...
        PCC, PCC_ci',pVal];
end