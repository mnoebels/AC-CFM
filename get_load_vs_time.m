function load_at_time = get_load_vs_time(result)
% GET_LOAD_VS_TIME converts cascade generations into times
%   Conversion is based on a probability distribution of outage time
%   differences.

%   AC-CFM
%   Copyright (c) 2020, Matthias Noebels
%   This file is part of AC-CFM.
%   Covered by the 3-clause BSD License (see LICENSE file for details).

    [load_at_time, t] = apply_recursion(digraph(result.G.Edges, result.G.Nodes), 'event', 1, zeros(1000, 3));
    
    initial = findnode(result.G, 'root');
    event = findnode(result.G, 'event');
    
    load_at_time = [result.G.Nodes.Load(initial) result.G.Nodes.Generators(initial) result.G.Nodes.Lines(initial); result.G.Nodes.Load(event) result.G.Nodes.Generators(event) result.G.Nodes.Lines(event); load_at_time(1:t, :)];
    load_at_time = [([0 0 1:t]).' load_at_time];
end

function [load_at_time, t] = apply_recursion(G, node, t, load_at_time)
    
    [eid, nid] = outedges(G, node);
    
    nid = unique(nid);
    
    if length(nid) > 1
        load_at_time(t:end, :) = 0;
        t_this = 0;
        for i = 1:length(nid)
            [load_at_time_isl, t_isl] = apply_recursion(G, nid(i), t, load_at_time);
            load_at_time(t:end, :) = load_at_time(t:end, :) + load_at_time_isl(t:end, :);
            t_this = max([t_this, t_isl]);
        end
        t = t_this;
    elseif length(nid) == 1
        if any(strcmp(G.Edges.Type(eid), 'OL'))
            delay = round(exprnd(19.5797));
            load_at_time(t:(t + delay - 1), 1) = G.Nodes.Load(findnode(G, node));
            load_at_time(t:(t + delay - 1), 2) = G.Nodes.Generators(findnode(G, node));
            load_at_time(t:(t + delay - 1), 3) = G.Nodes.Lines(findnode(G, node));
            t = t + delay;
        end
        
        [load_at_time, t] = apply_recursion(G, nid(1), t, load_at_time);
    else
        load_at_time(t:end, 1) = G.Nodes.Load(findnode(G, node));
        load_at_time(t:end, 2) = G.Nodes.Generators(findnode(G, node));
        load_at_time(t:end, 3) = G.Nodes.Lines(findnode(G, node));
    end

end