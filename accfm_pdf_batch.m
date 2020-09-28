function result = accfm_pdf_batch(network, pdf, alpha, number_of_scenarios, settings, output_file)
% batch processing for the AC Cascading Fault Model following a probability
% distribution
%   accfm_pdf_batch(network, pdf, alpha, number_of_scenarios, output_file) 
%   runs the AC-CFM in network, which can either be a matpower case struct
%   or a filename of a matpower case file for number_of_scenarios scenarios
%   with the initial contingency size following a probability distribution
%   specified in pdf and defined by parameter alpha. Currently, only 'zipf'
%   can be used as a pdf. The results are stored in the file specified in
%   output_file.

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
    if ~isnumeric(alpha)
        alpha = str2double(alpha);
    end
    if ~isnumeric(number_of_scenarios)
        number_of_scenarios = str2double(number_of_scenarios);
    end
    
    % model settings
    if ~exist('settings', 'var') || ~isstruct(settings)
        settings = get_default_settings;
    end
	
    % create the initial contingecies
	rng('shuffle');
	if strcmp(pdf, 'zipf')
		contingencies = zipfrnd(alpha, [number_of_scenarios 1]);
	end

	scenarios = cell(number_of_scenarios, 1);
    for i = 1:number_of_scenarios
        if iscell(network)
            n_branch = size(network{i}.branch, 1);
        else
            n_branch = size(network.branch, 1);
        end
        
		scenarios{i} = round(rand(contingencies(i), 1) * n_branch);
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
            fprintf('\t Completion: ');
            showTimeToCompletion;
            p = parfor_progress( number_of_scenarios );
        end
    
        result(number_of_scenarios) = struct('version', [], 'baseMVA', [], 'bus', [], 'gen', [], 'branch', [], 'gencost', [], 'gentype', [], 'genfuel', [], 'bus_name', [], 'order', [], 'et', [], 'success', [], 'iterations', [], 'bus_id', [], 'gen_id', [], 'branch_id', [], 'branch_tripped', [], 'bus_tripped', [], 'bus_uvls', [], 'bus_ufls', [], 'gen_tripped', [], 'load', [], 'generation_before', [], 'pf_count', [], 'G', [], 'ls_total', [], 'ls_ufls', [], 'ls_uvls', [], 'ls_vcls', [], 'ls_opf', [], 'ls_tripped', [], 'elapsed', []);
        
        parfor i = 1:number_of_scenarios
            %output progress if not running on cluster
            if ~isdeployed
                p = parfor_progress;
                showTimeToCompletion( p/100, [], [], startTime );
            else
                fprintf('Scenario %d', i);
            end
            
            network_i = runpf(network{i}, settings.mpopt);
            if ~network_i.success
                error("Initial PF of network failed, please check network case");
            end
            
            % run AC-CFM
            result(i) = accfm(network_i, struct('branches', scenarios{i}), settings);
        end
    end

    if exist('output_file', 'var')
        % save result struct
        save(output_file, 'result');
    end