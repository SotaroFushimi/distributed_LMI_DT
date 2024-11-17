function G = generate_delaunaygraph(n)
    % Function to generate the graph Laplacian matrix using a Delaunay triangulation
    % Input: n - Number of agents
    % Output: L - Laplacian matrix

    % Generate random positions for agents
    x = rand(n,1);
    y = rand(n,1);
    points = [x, y];
    
    % Generate Delaunay triangulation
    DT = delaunayTriangulation(points);
    
    % Create adjacency matrix
    Adj = zeros(n, n);
    edges = DT.edges; % Get edge list from Delaunay triangulation
    for k = 1:size(edges, 1)
        i = edges(k, 1);
        j = edges(k, 2);
        Adj(i, j) = 1;
        Adj(j, i) = 1;
    end
    
    % Calculate the Laplacian matrix
    D = diag(sum(Adj, 2)); % Degree matrix
    L = D - Adj   ;        % Laplacian matrix
    
    % Display the resulting graph
    G = graph(Adj);
    figure;
    plot(G, 'XData', x, 'YData', y);
    title('Delaunay Graph and Laplacian Matrix');
end
