function result = accfm_branch_scenarios(network, scenarios, settings)
% running AC-CFM for a set of initial contingencies
%   accfm_branch_scenarios(network, scenarios, settings) runs the AC-CFM in
%   network with for the initial contingencies given in scenarios.
%   result_cascade = accfm_branch_scenarios( ___ ) returns a struct
%   including the individual results for every scenario.

    % use default settings if not defined
    if ~exist('settings', 'var')
        settings = get_default_settings();
    end

    %% initialise variables
    number_of_scenarios = length(scenarios);

    result.scenarios = scenarios;
    contingency = zeros(number_of_scenarios, 1);
    lost_load_final = zeros(number_of_scenarios, 1);
    
    tripped_buses_in_generation = zeros(number_of_scenarios, settings.max_iterations);
    tripped_buses_in_scenario = zeros(number_of_scenarios, size(network.bus, 1), 'logical');
    
    tripped_lines_in_generation = zeros(number_of_scenarios, settings.max_iterations);
    tripped_lines_in_scenario = zeros(number_of_scenarios, size(network.branch, 1), 'logical');
    line_criticality = zeros(size(network.branch, 1), 1);
    
    ufls_buses = zeros(number_of_scenarios, size(network.bus, 1), 'logical');
    uvls_buses = zeros(number_of_scenarios, size(network.bus, 1), 'logical');
    tripped_gens_in_generation = zeros(number_of_scenarios, settings.max_iterations);
    tripped_gens_in_scenario = zeros(number_of_scenarios, size(network.gen, 1), 'logical');
    
    lost_load = zeros(number_of_scenarios, settings.max_iterations);
    ls_ufls = zeros(number_of_scenarios, 1);
    ls_uvls = zeros(number_of_scenarios, 1);
    ls_vcls = zeros(number_of_scenarios, 1);
    ls_opf = zeros(number_of_scenarios, 1);
    ls_tripped = zeros(number_of_scenarios, 1);
    load_at_time = cell(number_of_scenarios, 1);
    
    computational = zeros(number_of_scenarios, 2); 
    vcls = zeros(number_of_scenarios, 6);
    
    result.network = network;
    result.settings = settings;
    
    if settings.keep_networks_after_cascade
        networks_after_cascade = cell(number_of_scenarios, 1);
    end
    
    %% the loop
    % output progress if not running on cluster
    startTime = tic;
    if ~isdeployed
        fprintf('\t Completion: ');
        showTimeToCompletion;
        p = parfor_progress( number_of_scenarios );
    end
    
    % use parallel computing toolbox
    parfor i = 1:number_of_scenarios
        %output progress if not running on cluster
        if ~isdeployed
            p = parfor_progress;
            showTimeToCompletion( p/100, [], [], startTime );
        else
            fprintf('Scenario %d', i);
        end
        
        % AC-CFM only produces an error if iteration limit is reached
        try
            % run the model
            result_cascade = accfm(network, struct('branches', scenarios{i}), settings);

            % get results and store them in the result arrays
            contingency(i) = length(scenarios{i});
            tripped_lines_in_generation(i, :) = sum(result_cascade.branch_tripped);
            tripped_lines_in_scenario(i, :) = sum(result_cascade.branch_tripped, 2);
            line_criticality = line_criticality + sum(result_cascade.branch_tripped, 2);
            
            tripped_buses_in_generation(i, :) = sum(result_cascade.bus_tripped);
            tripped_buses_in_scenario(i, :) = sum(result_cascade.bus_tripped, 2);
            
            ufls_buses_i = zeros(size(result_cascade.bus, 1), 1);
            ufls_buses_i(sum(result_cascade.bus_ufls, 2) > 0) = 1;
            ufls_buses(i, :) = ufls_buses_i;
            
            uvls_buses_i = zeros(size(result_cascade.bus, 1), 1);
            uvls_buses_i(sum(result_cascade.bus_uvls, 2) > 0) = 1;
            uvls_buses(i, :) = uvls_buses_i;
            
            tripped_gens_in_generation(i, :) = sum(result_cascade.gen_tripped);
            tripped_gens_in_scenario(i, :) = sum(result_cascade.gen_tripped, 2);
            
            ls_i = 1 - result_cascade.load / sum(network.bus(:, 3));
            lost_load(i, :) = ls_i;
            
            ls_ufls(i) = result_cascade.ls_ufls;
            ls_uvls(i) = result_cascade.ls_uvls;
            ls_vcls(i) = result_cascade.ls_vcls;
            ls_opf(i) = result_cascade.ls_opf;
            ls_tripped(i) = result_cascade.ls_tripped;
            lost_load_final(i) = ls_i(end);
            
            rng(100000 + i);
            load_at_time{i} = get_load_vs_time(result_cascade);
            
            vc_edges = find(strcmp(result_cascade.G.Edges.Type, 'VC'));
            vc_nodes = findnode(result_cascade.G, result_cascade.G.Edges.EndNodes(vc_edges, 2));
            vc_island_sizes = result_cascade.G.Nodes.Buses(vc_nodes);

            opf_edges = find(strcmp(result_cascade.G.Edges.Type, 'OPF'));
            opf_nodes = findnode(result_cascade.G, result_cascade.G.Edges.EndNodes(opf_edges, 2));
            opf_island_sizes = result_cascade.G.Nodes.Buses(opf_nodes);

            dc_edges = find(strcmp(result_cascade.G.Edges.Type, 'DC'));
            dc_lines = result_cascade.G.Edges.Weight(dc_edges);
            
            computational(i, :) = [result_cascade.pf_count result_cascade.elapsed];
            vcls(i, :) = [length(vc_edges) mean(vc_island_sizes) length(opf_edges) mean(opf_island_sizes) length(dc_edges) mean(dc_lines)];

            if settings.keep_networks_after_cascade
                networks_after_cascade{i} = struct('version', result_cascade.version, 'baseMVA', result_cascade.baseMVA, 'bus', result_cascade.bus, 'branch', result_cascade.branch, 'gen', result_cascade.gen, 'gencost', result_cascade.gencost);
            end
            
            if isdeployed
                fprintf(' (%.2f%%)\n', ls_i(end) * 100);
            end
        catch
            % iteration limit reached
            
            lost_load_final(i) = -1;
        end
    end
    
    toc(startTime);
    
    result.contingency = contingency;
    result.tripped_lines_in_scenario = tripped_lines_in_scenario;
    result.tripped_lines_in_generation = tripped_lines_in_generation;
    result.line_criticality = line_criticality;
    
    result.tripped_buses_in_scenario = tripped_buses_in_scenario;
    result.tripped_buses_in_generation = tripped_buses_in_generation;
    result.ufls_buses = ufls_buses;
    result.uvls_buses = uvls_buses;
    result.tripped_gens_in_scenario = tripped_gens_in_scenario;
    result.tripped_gens_in_generation = tripped_gens_in_generation;
    result.lost_load = lost_load;
    result.ls_ufls = ls_ufls;
    result.ls_uvls = ls_uvls;
    result.ls_vcls = ls_vcls;
    result.ls_opf = ls_opf;
    result.ls_tripped = ls_tripped;
    result.lost_load_final = lost_load_final;
    result.load_at_time = load_at_time;
    result.computational = computational;
    result.vcls = vcls;
    
    if settings.keep_networks_after_cascade
        result.networks_after_cascade = networks_after_cascade;
    end
end