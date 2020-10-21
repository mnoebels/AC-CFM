function network = disp2load(network)
%DISP2LOAD converts dispatchable loads to fixed loads.
%   It is the inverse of the load2disp function supplied with Matpower.

    define_constants;

    % get all dispatchable loads from the gen matrix
    loads = network.gen(isload(network.gen), :);
    
    % convert to non-dispatchable loads
    network.bus(ismember(network.bus(:, BUS_I), loads(:, GEN_BUS)), [PD QD]) = network.bus(ismember(network.bus(:, BUS_I), loads(:, GEN_BUS)), [PD QD]) - loads(:, [PG QG]);
    
    % remove from gen and gencost matrix
    network.gencost(isload(network.gen), :) = [];
    network.gen(isload(network.gen), :) = [];
end

