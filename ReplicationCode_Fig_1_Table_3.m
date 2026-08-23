%% Replication File to the paper "Robust Forecast Evaluation Under Extreme Shocks"
% Authors: Fabrizio Iacone and Andrea Viselli
% Description: Real GDP growth forecast evaluation (SPF vs. naive benchmark)
% -------------------------------------------------------------------------

clc; clear; close all;

%% 1. Setup and Data Loading
% Load the dataset containing SPF forecast error statistics
% Variables in .mat file: Realiz1 (observed values), SPFfor_Step1 (SPF forecasts)
load('SPF_Error_Statistics.mat');

% Generate quarterly date vector from 1968:Q4 to 2025:Q4
% Note: Using 1st of the month for the respective quarters
dates = datetime(1968, 11, 1) : calmonths(3) : datetime(2025, 11, 1); 
dates = dates';

%% 2. Data Processing and Trimming (2000:Q1 - 2024:Q4)
% Compute forecast errors:
% 1: SPF median forecast
err_spf = Realiz1 - SPFfor_Step1;

% 2: Naive benchmark (predicting 0 growth)
err_naive = Realiz1 - 0;

% Trim the sample to the period of interest: 2000:Q1 to 2024:Q4 (100 observations)
% Index 126 corresponds to 2000:Q1 in the generated date vector
t0 = 126; 
t1 = t0 + 100 - 1; 

err_spf = err_spf(t0:t1); 
err_naive = err_naive(t0:t1);
dates = dates(t0:t1);

nobs = size(err_spf, 1);

% Define the indices for normal and Covid-19 periods
% Covid-19 observations are 2020:Q1 to 2020:Q3 (observations 81 to 83)
covid_idx = 81:83;
normal_idx = [1:80, 84:nobs];

%% 3. Figure 1: Forecast Errors Plot
% Plotting the quarterly SPF one-step-ahead forecast errors vs naive benchmark
figure('Name', 'Figure 1: Forecast Errors');
plot(dates, err_spf, '-k', 'LineWidth', 1.2); hold on;
plot(dates, err_naive, ':k', 'LineWidth', 1.2);
title('Quarterly growth of Real GDP forecast errors (2000 - 2024)');
ylabel('Forecast errors');
legend('SPF', 'Naive benchmark', 'Location', 'northwest');
legend boxoff
grid on; hold off;

%% 4. Table 3: Summary Statistics
% Compute Average and MSE for the whole sample and excluding the pandemic

% Whole sample statistics
avg_spf_all = sum(err_spf) / nobs;
avg_naive_all = sum(err_naive) / nobs;
mse_spf_all = (err_spf' * err_spf) / nobs;
mse_naive_all = (err_naive' * err_naive) / nobs;

% Excluding the pandemic (observations 81-83)
nobs_no_covid = nobs - length(covid_idx);
avg_spf_nc = sum(err_spf(normal_idx)) / nobs_no_covid;
avg_naive_nc = sum(err_naive(normal_idx)) / nobs_no_covid;
mse_spf_nc = (err_spf(normal_idx)' * err_spf(normal_idx)) / nobs_no_covid;
mse_naive_nc = (err_naive(normal_idx)' * err_naive(normal_idx)) / nobs_no_covid;

%% 5. Table 3: Ordinary Diebold-Mariano (DM) Test
% Compute the standard DM test over the full sample (including Covid-19)
loss_diff = err_spf.^2 - err_naive.^2; 

% Bandwidth: set to floor(T^(1/2)) as specified in Section 4
bw = floor(nobs^(1/2));  	
b = bw / nobs;							

% Compute Fixed-Smoothing critical values (Coroneo and Iacone, 2020)
cv_fixed_smoothing = 1.96 + 2.9694*b + 0.4160*(b^2) - 0.5324*(b^3);							

% Long Run Variance (LRV) estimation using Bartlett kernel
demeaned_loss_diff = loss_diff - mean(loss_diff);							
long_run_var = sum(demeaned_loss_diff.^2) / nobs;							
for jcov = 1:bw							
    gamma_j = sum(demeaned_loss_diff(1+jcov:nobs) .* demeaned_loss_diff(1:nobs-jcov)) / nobs;
    long_run_var = long_run_var + 2 * ((bw - jcov) / bw) * gamma_j;							
end	

% Ordinary DM test statistic computation
dm_stat_ordinary = abs(sqrt(nobs) * mean(loss_diff) / sqrt(long_run_var));							

%% 6. Table 3: Winsorized Diebold-Mariano (DM*) Test
% Winsorizing the loss function for the Covid-19 period (observations 81-83)
% using the infimum and supremum of the normal sample.

loss_diff_wins = loss_diff; % Create a copy for winsorization
min_loss_normal = min(loss_diff(normal_idx));
max_loss_normal = max(loss_diff(normal_idx));

% Apply winsorization to Covid-19 quarters
for i = covid_idx
    if loss_diff_wins(i) < min_loss_normal 
        loss_diff_wins(i) = min_loss_normal; 
    end
    if loss_diff_wins(i) > max_loss_normal 
        loss_diff_wins(i) = max_loss_normal; 
    end
end

% Recompute LRV with the winsorized loss differential							
demeaned_loss_wins = loss_diff_wins - mean(loss_diff_wins);							
long_run_var_wins = sum(demeaned_loss_wins.^2) / nobs;							
for jcov = 1:bw							
    gamma_j_wins = sum(demeaned_loss_wins(1+jcov:nobs) .* demeaned_loss_wins(1:nobs-jcov)) / nobs;
    long_run_var_wins = long_run_var_wins + 2 * ((bw - jcov) / bw) * gamma_j_wins;							
end	

% Winsorized DM test statistic (DM*) computation							
dm_stat_wins = abs(sqrt(nobs) * mean(loss_diff_wins) / sqrt(long_run_var_wins));							

% Display the combined Table 3
fprintf('===============================================================================\n');
fprintf('Table 3: Forecast error summary and Diebold-Mariano test\n');
fprintf('-------------------------------------------------------------------------------\n');
fprintf('%-22s %-10s %10s %10s %10s %10s\n', ...
    'Sample', 'Forecast', 'AFE', 'MSE', 'DM', 'DM*');
fprintf('-------------------------------------------------------------------------------\n');
fprintf('%-22s %-10s %10.4f %10.4f %10.3f %10.3f\n', ...
    'Whole sample', 'SPF', avg_spf_all, mse_spf_all, dm_stat_ordinary, dm_stat_wins);
fprintf('%-22s %-10s %10.4f %10.4f %10s %10s\n', ...
    '', 'Naive', avg_naive_all, mse_naive_all, '', '');
fprintf('%-22s %-10s %10.4f %10.4f %10s %10s\n', ...
    'Excluding pandemic', 'SPF', avg_spf_nc, mse_spf_nc, '', '');
fprintf('%-22s %-10s %10.4f %10.4f %10s %10s\n', ...
    '', 'Naive', avg_naive_nc, mse_naive_nc, '', '');
fprintf('-------------------------------------------------------------------------------\n');
fprintf(['Notes: The table reports average forecast errors (AFE) and mean squared errors for the\n' ...
    '       SPF Real GDP nowcasts and naive benchmark of zero growth. The pandemic period\n' ...
    '       refers to 2020Q1-2020Q3. DM* denotes the winsorised Diebold-Mariano statistic.\n' ...
    '       The 5%% critical value is %.3f (see Coroneo and Iacone, 2020).\n'], ...
    cv_fixed_smoothing);
fprintf('===============================================================================\n\n');

%% 7. Monte Carlo Calibration Parameters (Section 3 of Paper)
% Fit an AR(1) model to the loss differential to calibrate the Monte Carlo simulation
fprintf('--- Monte Carlo Calibration (AR(1) fit to pre-Covid data) ---\n');

% Standardizing variance for the first 80 observations
% Dividing by sqrt(68) aligns the variance with the MC design where sigma^2 = 1.
mdl = arima(1, 0, 0); 
EstMdl = estimate(mdl, loss_diff(1:80) / sqrt(68), 'Display', 'off');

% Print empirical parameter equivalents
disp('Estimated AR(1) Coefficient (rho):');
disp(EstMdl.AR{1});

disp('Mean parameter estimate equivalent:');
fprintf('Pre-Covid Mean (1:80): %f\n', mean(loss_diff(1:80)) / sqrt(68));
fprintf('Covid Mean (81:83):    %f\n', mean(loss_diff(81:83)) / sqrt(68));
