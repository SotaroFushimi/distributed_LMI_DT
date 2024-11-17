function G = generate_diskgraph(n)

    % arguments
    %     prob = 0.3 % probability to generate edges
    % end
    
    Adj = zeros(n,n); % adjacency matrix
    
    r = 0.2;
    x = rand(n,1);
    y = rand(n,1);
    for i = 1:n
        for j = 1:n
            if norm([x(i,1)-x(j,1),y(i,1)-y(j,1)]) <= r && i>j
                Adj(i, j) = 1;
                Adj(j,i) = 1;
            end
        end
    end
    % diag(laplacian(graph(Adj)))

    G = graph(Adj);