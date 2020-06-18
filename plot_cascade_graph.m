function H = plot_cascade_graph(result)
% PLOT_GRAPH plots the cascade graph returned by the AC-CFM
%   plot_cascade_graph(result) plots the cascade graph in result.
%   H = plot_cascade_graph( ___ ) returns the GraphPlot object.

    % group isolated nodes for better visibility
    ISLEdges = find(ismember(result.G.Edges.Type, 'ISL'));
    ISLEdges = ISLEdges(result.G.Nodes.Buses(findnode(result.G, result.G.Edges.EndNodes(ISLEdges, 2))) == 1 & strcmp(result.G.Nodes.Type(findnode(result.G, result.G.Edges.EndNodes(ISLEdges, 2))), 'failure'));
    ISLStartNodes = findnode(result.G, result.G.Edges.EndNodes(ISLEdges, 1));
    ISLEndNodes = findnode(result.G, result.G.Edges.EndNodes(ISLEdges, 2));

    [ISLCombinedStartNodes, ISLCombinedFirstEdge] = unique(ISLStartNodes);
    ISLCombinedCounts = accumarray(ISLStartNodes, 1);

    result.G.Edges.Type(ISLEdges(ISLCombinedFirstEdge)) = {'(ISL)'};
    result.G.Nodes.Buses(findnode(result.G, result.G.Edges.EndNodes(ISLEdges(ISLCombinedFirstEdge), 2))) = ISLCombinedCounts(ISLCombinedStartNodes);

    result.G = rmnode(result.G, setdiff(ISLEndNodes, findnode(result.G, result.G.Edges.EndNodes(ISLEdges(ISLCombinedFirstEdge), 2))));

    ColouredEdges = ismember(result.G.Edges.Type, 'UFLS') | ismember(result.G.Edges.Type, 'UVLS') | ismember(result.G.Edges.Type, 'VC') | ismember(result.G.Edges.Type, 'XL');
    WidthEdges = ~ismember(result.G.Edges.Type, 'UFLS') & ~ismember(result.G.Edges.Type, 'UVLS') & ~ismember(result.G.Edges.Type, 'VC') & ~ismember(result.G.Edges.Type, 'XL');

    ECData = zeros(numedges(result.G), 1);
    ECData(ColouredEdges) = result.G.Edges.Weight(ColouredEdges) ./ result.G.Edges.Base(ColouredEdges) * 100;

    LWidth = repmat(2, numedges(result.G), 1);
    LWidth(WidthEdges) = result.G.Edges.Weight(WidthEdges);

    figure;
    H = plot(result.G, 'NodeLabel', result.G.Nodes.Buses, 'EdgeLabel', result.G.Edges.Type, 'LineWidth', LWidth, 'EdgeCData', ECData, 'Layout', 'layered');
    highlight(H, find(ismember(result.G.Nodes.Type, 'root')), 'NodeColor', 'magenta');
    highlight(H, find(ismember(result.G.Nodes.Type, 'success')), 'NodeColor', 'green');
    highlight(H, find(ismember(result.G.Nodes.Type, 'failure')), 'NodeColor', 'red');
    colorbar;
    caxis([0 100]);
end

