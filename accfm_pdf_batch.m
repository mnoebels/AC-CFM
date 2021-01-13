function result = accfm_pdf_batch(network, distribName, distribParams, number_of_scenarios, settings, output_file)
% batch processing for the AC Cascading Fault Model following a probability
% distribution
%   accfm_pdf_batch(network, distribName, distribParams, number_of_scenarios, output_file) 
%   runs the AC-CFM in network, which can either be a matpower case struct
%   or a filename of a matpower case file for number_of_scenarios scenarios
%   with the initial contingency size following a probability distribution
%   specified by distribName and defined by parameters distribParams. Call 
%   randraw for a list of available distributions.
%   The results are stored in the file specified in output_file.
%   You may want to run rng('shuffle') before calling this function.

%   AC-CFM
%   Copyright (c) 2020, Matthias Noebels
%   This file is part of AC-CFM.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

    % load the network if it is given as a file
    if (isstring(network) || ischar(network)) && isfile(network)
        variableInfo = who('-file', network);
        if ismember('network', variableInfo)
            load(network, 'network');
        else
            [network, info] = loadcase(network);
            if ~info
                error("Invalid network case file");
            end
        end
    elseif iscell(network)
        number_of_scenarios = length(network);
    elseif (isstring(network) || ischar(network)) && ~isfile(network)
        error("Network file not found");
    end

    % convert non-numeric input parameters
    if ~isnumeric(distribParams)
        distribParams = str2num(distribParams);
    end
    if ~isnumeric(number_of_scenarios)
        number_of_scenarios = str2num(number_of_scenarios);
    end
    
    % model settings
    if ~exist('settings', 'var') || ~isstruct(settings)
        settings = get_default_settings;
    end
	
    % create the initial contingecies
	contingencies = round(randraw(distribName, distribParams, number_of_scenarios));

	scenarios = cell(number_of_scenarios, 1);
    for i = 1:number_of_scenarios
        if iscell(network)
            n_branch = size(network{i}.branch, 1);
        else
            n_branch = size(network.branch, 1);
        end
        
        if contingencies{i} > n_branch
            contingencies{i} = n_branch;
        end
        
		scenarios{i} = randperm(n_branch, contingencies(i));
    end
    
    if ~iscell(network)
        % run an initial PF to see if the case is valid
        network = runpf(network, settings.mpopt);
        if ~network.success
            error("Initial PF of network failed, please check network case");
        end

        % run AC-CFM
        result = accfm_branch_scenarios(network, scenarios, settings);
    else
        % output progress if not running on cluster
        startTime = tic;
        if ~isdeployed
            parfor_progress( number_of_scenarios );
        end
        
        result(number_of_scenarios) = struct('version', [], 'baseMVA', [], 'bus', [], 'gen', [], 'branch', [], 'gencost', [], 'success', [], 'branch_tripped', [], 'bus_tripped', [], 'bus_uvls', [], 'bus_ufls', [], 'gen_tripped', [], 'load', [], 'pf_count', [], 'ls_total', [], 'ls_ufls', [], 'ls_uvls', [], 'ls_vcls', [], 'ls_opf', [], 'ls_tripped', [], 'elapsed', []);
    
        parfor i = 1:number_of_scenarios
            % output progress if not running on cluster
            if ~isdeployed
                parfor_progress;
            else
                fprintf('Scenario %d', i);
            end
            
            network_i = runpf(network{i}, settings.mpopt);
            if ~network_i.success
                error("Initial PF of network failed, please check network case");
            end
            
            try
                % run AC-CFM
                result_i = accfm(network_i, struct('branches', scenarios{i}), settings);
                result(i).version = result_i.version;
                result(i).baseMVA = result_i.baseMVA;
                result(i).bus = result_i.bus;
                result(i).branch = result_i.branch;
                result(i).gen = result_i.gen;
                result(i).gencost = result_i.gencost;
                result(i).success = result_i.success;
                result(i).branch_tripped = result_i.branch_tripped;
                result(i).bus_tripped = result_i.bus_tripped;
                result(i).gen_tripped = result_i.gen_tripped;
                result(i).bus_uvls = result_i.bus_uvls;
                result(i).bus_ufls = result_i.bus_ufls;
                result(i).load = result_i.load;
                result(i).pf_count = result_i.pf_count;
                result(i).ls_total = result_i.ls_total;
                result(i).ls_ufls = result_i.ls_ufls;
                result(i).ls_uvls = result_i.ls_uvls;
                result(i).ls_vcls = result_i.ls_vcls;
                result(i).ls_opf = result_i.ls_opf;
                result(i).ls_tripped = result_i.ls_tripped;
                result(i).elapsed = result_i.elapsed;
                
                result(i).contingency = scenarios{i};
            catch
                result(i).success = 0;
            end
        end
    end

    if exist('output_file', 'var')
        % save result struct
        save(output_file, 'result');
    end