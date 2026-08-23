function [rejDM, rejDMstar] = run_spf_mc( ...
    nobs, nrounds, sigmasq, phi, kappa, theta, m, ...
    bw, DMcv, seed, reset_rng_each_design)

    % ---------------------------------------------------------------------
    % Reset random seed if requested
    % ---------------------------------------------------------------------

    if reset_rng_each_design
        rng(seed);
    end

    % ---------------------------------------------------------------------
    % Storage for Monte Carlo rejection indicators
    % ---------------------------------------------------------------------

    pass1 = zeros(nrounds, 1);   % Rejection indicator for the regular DM test
    pass2 = zeros(nrounds, 1);   % Rejection indicator for the winsorized DM* test

    % ---------------------------------------------------------------------
    % Pre-allocation for the simulated DGP
    % ---------------------------------------------------------------------

    eps = zeros(nobs, 1);        % Innovations
    ud  = zeros(nobs, 1);        % AR(1) loss-differential component
    LL  = zeros(nobs, 1);        % Loss differential
    dLL = zeros(nobs, 1);        % Demeaned loss differential

    % ---------------------------------------------------------------------
    % Monte Carlo loop
    % ---------------------------------------------------------------------

    for index = 1:nrounds

        % -----------------------------------------------------------------
        % DGP
        % -----------------------------------------------------------------
        % Simulate an AR(1) process for the loss-differential component.

        eps = randn(nobs, 1) * sqrt(sigmasq);

        ud(1) = sqrt(1 / (1 - phi^2)) * eps(1);

        for t = 2:nobs
            ud(t) = phi * ud(t - 1) + eps(t);
        end

        % Add the average loss differential.

        LL = ud + kappa;

        % Add the Covid-period loss differential.
        % The Covid-shock period starts at t = 81 and lasts for m observations.

        for t = 81:(81 + m - 1)
            LL(t) = LL(t) + theta;
        end

        % -----------------------------------------------------------------
        % Regular DM test
        % -----------------------------------------------------------------

        % Long-run variance estimate.

        dLL = LL - mean(LL);

        lrv = sum(dLL .* dLL) / nobs;

        for jcov = 1:bw
            lrv = lrv ...
                + 2 * ((bw - jcov) / bw) ...
                * sum(dLL((1 + jcov):nobs) .* dLL(1:(nobs - jcov))) / nobs;
        end

        % DM test statistic.

        ts = sqrt(nobs) * mean(LL) / sqrt(lrv);

        if abs(ts) > DMcv
            pass1(index) = 1;
        end

        % -----------------------------------------------------------------
        % Winsorizing the loss function
        % -----------------------------------------------------------------
        % The Covid-shock observations are winsorized using the minimum and
        % maximum loss differentials observed outside the Covid-shock period.
        %
        % The non-Covid periods are:
        %   - observations 1 to 80
        %   - observations 81 + m to 100

        minLL = min([min(LL(1:80)) min(LL((81 + m):100))]);
        maxLL = max([max(LL(1:80)) max(LL((81 + m):100))]);

        for cvd = 1:m

            if LL(80 + cvd) < minLL
                LL(80 + cvd) = minLL;
            end

            if LL(80 + cvd) > maxLL
                LL(80 + cvd) = maxLL;
            end

        end

        % -----------------------------------------------------------------
        % Long-run variance estimate after winsorization
        % -----------------------------------------------------------------

        dLL = LL - mean(LL);

        lrv = sum(dLL .* dLL) / nobs;

        for jcov = 1:bw
            lrv = lrv ...
                + 2 * ((bw - jcov) / bw) ...
                * sum(dLL((1 + jcov):nobs) .* dLL(1:(nobs - jcov))) / nobs;
        end

        % -----------------------------------------------------------------
        % DM* test statistic
        % -----------------------------------------------------------------

        tsDM = abs(sqrt(nobs) * mean(LL) / sqrt(lrv));

        if abs(tsDM) > DMcv
            pass2(index) = 1;
        end

    end

    % ---------------------------------------------------------------------
    % Empirical rejection frequencies
    % ---------------------------------------------------------------------

    rejDM     = mean(pass1);
    rejDMstar = mean(pass2);

end