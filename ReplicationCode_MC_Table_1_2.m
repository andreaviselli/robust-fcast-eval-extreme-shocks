%% Replication File to the paper "Robust Forecast Evaluation Under Extreme Shocks"
% Authors: Fabrizio Iacone and Andrea Viselli
% Description: Monte Carlo study (Section 3)
% -------------------------------------------------------------------------

clc; clear; close all;

% -------------------------------------------------------------------------
% Monte Carlo design
% -------------------------------------------------------------------------

seed    = 0;
nobs    = 100;
nrounds = 10000;

% Resetting the seed before each experiment
reset_rng_each_design = true;

rng(seed);

% -------------------------------------------------------------------------
% DGP specifications
% -------------------------------------------------------------------------
% The simulated loss differential is:
%
%   LL_t = u_t + kappa + theta * 1{t is in Covid-shock period},
%
% where u_t follows an AR(1) process.
%
% The Covid-shock period starts at observation 81 and lasts for m periods.
% The current parameter choices correspond to those in the paper.

sigmasq = 1;       % Innovation variance; specification from the SPF setting

% -------------------------------------------------------------------------
% Parameter choices used in the Monte Carlo experiment
% -------------------------------------------------------------------------

% Size experiment: no average loss differential and no Covid shock.
phi_grid_size = [0, 0.25, 0.5];
m_grid_size   = [1, 3];

kappa_size = 0;
theta_size = 0;

% Power / Covid-shock experiments.
phi_power = 0.0;

m_grid_power = [1, 3];

design_grid = [
     0.0   -40
     0.0   -80
    -0.4     0
    -0.8     0
    -0.4   -40
    -0.8   -40
    -0.4   -80
    -0.8   -80
     0.4   -40
     0.8   -40
];

% -------------------------------------------------------------------------
% Global quantities used in the long-run variance and critical value
% -------------------------------------------------------------------------

bw = floor(nobs^(1/2));      % Bandwidth used in the LRV estimate
b  = bw / nobs;              % Fixed-b ratio

DMcv = 1.96 ...
     + 2.9694 * b ...
     + 0.4160 * (b^2) ...
     - 0.5324 * (b^3);

% -------------------------------------------------------------------------
% Storage for Monte Carlo results
% -------------------------------------------------------------------------

% Size results:
% columns are m, phi, DM, DM*
size_results = zeros(numel(m_grid_size) * numel(phi_grid_size), 4);

% Power / Covid-shock results:
% columns are m, kappa, theta, DM, DM*
power_results = zeros(numel(m_grid_power) * size(design_grid, 1), 5);

% -------------------------------------------------------------------------
% Run size experiment
% -------------------------------------------------------------------------

row = 0;

for im = 1:numel(m_grid_size)

    m = m_grid_size(im);

    for iphi = 1:numel(phi_grid_size)

        phi = phi_grid_size(iphi);

        [rejDM, rejDMstar] = run_mc_spec( ...
            nobs, nrounds, sigmasq, phi, kappa_size, theta_size, m, ...
            bw, DMcv, seed, reset_rng_each_design);

        row = row + 1;

        size_results(row, :) = [m, phi, rejDM, rejDMstar];

    end

end

% -------------------------------------------------------------------------
% Run power / Covid-shock experiment
% -------------------------------------------------------------------------

row = 0;

for im = 1:numel(m_grid_power)

    m = m_grid_power(im);

    for idesign = 1:size(design_grid, 1)

        kappa = design_grid(idesign, 1);
        theta = design_grid(idesign, 2);

        [rejDM, rejDMstar] = run_mc_spec( ...
            nobs, nrounds, sigmasq, phi_power, kappa, theta, m, ...
            bw, DMcv, seed, reset_rng_each_design);

        row = row + 1;

        power_results(row, :) = [m, kappa, theta, rejDM, rejDMstar];

    end

end

% -------------------------------------------------------------------------
% Display Monte Carlo settings and empirical size/power
% -------------------------------------------------------------------------

lineLong  = repmat('=', 1, 78);
lineShort = repmat('-', 1, 78);

fprintf('\n');
fprintf('%s\n', lineLong);
fprintf(' SPF Forecast Evaluation: Monte Carlo Results\n');
fprintf('%s\n\n', lineLong);

fprintf('Monte Carlo settings\n');
fprintf('%s\n', lineShort);
fprintf(' nobs                  = %d\n', nobs);
fprintf(' nrounds               = %d\n', nrounds);
fprintf(' seed                  = %d\n', seed);
fprintf(' sigmasq               = %.4f\n', sigmasq);
fprintf(' bandwidth             = %d\n', bw);
fprintf(' fixed-b ratio         = %.4f\n', b);
fprintf(' DM critical value     = %.4f\n', DMcv);
fprintf('\n');

fprintf('Parameter interpretation\n');
fprintf('%s\n', lineShort);
fprintf(' kappa       = average loss differential\n');
fprintf(' theta    = additional loss differential during the Covid-shock period\n');
fprintf(' m        = length of the Covid-shock period\n');
fprintf(' phi      = AR(1) coefficient of the loss-differential component\n');
fprintf('\n');

% -------------------------------------------------------------------------
% Table 1: empirical size
% -------------------------------------------------------------------------

fprintf('%s\n', lineLong);
fprintf(' Table 1. Empirical size: kappa = 0, theta = 0\n');
fprintf('%s\n', lineLong);

fprintf(' %-8s %10s %14s %14s\n', 'm', 'phi', 'DM', 'DM*');
fprintf('%s\n', lineShort);

for i = 1:size(size_results, 1)

    fprintf(' %-8.0f %10.2f %14.4f %14.4f\n', ...
        size_results(i, 1), ...
        size_results(i, 2), ...
        size_results(i, 3), ...
        size_results(i, 4));

end

fprintf('\n');

% -------------------------------------------------------------------------
% Table 2: empirical size/power under alternative designs
% -------------------------------------------------------------------------

fprintf('%s\n', lineLong);
fprintf(' Table 2. Empirical size/power: phi = %.2f\n', phi_power);
fprintf('%s\n', lineLong);

for im = 1:numel(m_grid_power)

    m = m_grid_power(im);

    fprintf('\n');
    fprintf(' Covid-shock length: m = %d\n', m);
    fprintf('%s\n', lineShort);
    fprintf(' %-8s %10s %10s %14s %14s\n', 'm', 'kappa', 'theta', 'DM', 'DM*');
    fprintf('%s\n', lineShort);

    rows_m = power_results(power_results(:, 1) == m, :);

    for i = 1:size(rows_m, 1)

        fprintf(' %-8.0f %10.2f %10.0f %14.4f %14.4f\n', ...
            rows_m(i, 1), ...
            rows_m(i, 2), ...
            rows_m(i, 3), ...
            rows_m(i, 4), ...
            rows_m(i, 5));

    end

end

fprintf('\n');
fprintf('%s\n', lineLong);
fprintf(' End \n');
fprintf('%s\n\n', lineLong);
