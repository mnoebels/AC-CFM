function settings = get_default_settings()
    settings = struct;
    settings.verbose = 0;
    settings.mpopt = mpoption('verbose', 0, 'model', 'AC', 'out.all', 0);
    settings.max_iterations = 100;
    settings.enforce_q = 1;
    settings.uvls_per_step = 0.05;
    settings.uvls_max_steps = 5;
    settings.dP_limit = 0.15;
    settings.dSlack_limit = 0.01;
    settings.P_overhead = 0.1;
    settings.Q_tolerance = 0.1;
    settings.DC_fallback = 1;
end

