%% Replication File to the paper "Robust Forecast Evaluation Under Extreme Shocks"
% Authors: Fabrizio Iacone and Andrea Viselli
% Description: Monte Carlo study (Section 3)
% -------------------------------------------------------------------------

clc; clear; close all;

% -------------------------------------------------------------------------
% Monte Carlo settings
% -------------------------------------------------------------------------

seed    = 0;
nobs    = 100;
nrounds = 10000;
sigmasq = 1;

% Empirically calibrated design reported at the end of Section 3 of the paper.
% The shock period starts at observation 81 and lasts for m observations.
kappa = -0.8;
theta = -80;
phi   = 0.25;
m     = 3;

% Reset the random-number generator for run_mc_spec.
reset_rng = true;

% Bartlett bandwidth and fixed-smoothing critical value.
bw = floor(nobs^(1/2));
b  = bw / nobs;

DMcv = 1.96 ...
     + 2.9694 * b ...
     + 0.4160 * (b^2) ...
     - 0.5324 * (b^3);

% -------------------------------------------------------------------------
% Run the calibrated experiment
% -------------------------------------------------------------------------

[rejDM, rejDMstar] = run_mc_spec( ...
    nobs, nrounds, sigmasq, phi, kappa, theta, m, ...
    bw, DMcv, seed, reset_rng);

% -------------------------------------------------------------------------
% Display results
% -------------------------------------------------------------------------

lineLong  = repmat('=', 1, 72);
lineShort = repmat('-', 1, 72);

fprintf('\n%s\n', lineLong);
fprintf(' Empirically Calibrated Monte Carlo Experiment\n');
fprintf('%s\n', lineLong);
fprintf(' nobs                  = %d\n', nobs);
fprintf(' nrounds               = %d\n', nrounds);
fprintf(' seed                  = %d\n', seed);
fprintf(' sigmasq               = %.4f\n', sigmasq);
fprintf(' kappa                 = %.2f\n', kappa);
fprintf(' theta                 = %.0f\n', theta);
fprintf(' phi                   = %.2f\n', phi);
fprintf(' m                     = %d\n', m);
fprintf(' bandwidth             = %d\n', bw);
fprintf(' DM critical value     = %.4f\n', DMcv);
fprintf('%s\n', lineShort);
fprintf(' %-18s %14s\n', 'Test', 'Rejection rate');
fprintf('%s\n', lineShort);
fprintf(' %-18s %14.4f\n', 'DM', rejDM);
fprintf(' %-18s %14.4f\n', 'DM*', rejDMstar);
fprintf('%s\n\n', lineLong);
