function G = generate_cliquegraph(n, l)
    % Function to generate a graph with n agents divided into l cliques of random sizes
    % Input: n - Total number of agents
    %        l - Number of cliques
    % Output: G - Graph object representing the generated graph

    % Initialize adjacency matrix for n agents
    Adj = zeros(n, n);

    % Generate random clique sizes that sum up to n
    cliqueSizes = generate_random_partition(n, l);

    % Initialize agent indices and store each clique's indices
    agentIndex = 1;
    cliqueIndices = cell(l, 1);
    
    % Create cliques and populate the adjacency matrix
    for i = 1:l
        cliqueSize = cliqueSizes(i);
        agentsInClique = agentIndex:(agentIndex + cliqueSize - 1);
        cliqueIndices{i} = agentsInClique;

        % Fully connect agents within each clique
        for p = agentsInClique
            for q = agentsInClique
                if p ~= q
                    Adj(p, q) = 1;
                    Adj(q, p) = 1;
                end
            end
        end

        % Move to the next set of agents
        agentIndex = agentIndex + cliqueSize;
    end

    % Connect each consecutive clique with a single edge
    for i = 1:l-1
        % Connect the last node of the current clique to the first node of the next clique
        currentCliqueLast = cliqueIndices{i}(end);
        nextCliqueFirst = cliqueIndices{i+1}(1);
        Adj(currentCliqueLast, nextCliqueFirst) = 1;
        Adj(nextCliqueFirst, currentCliqueLast) = 1;
    end

    % Calculate the Laplacian matrix
    D = diag(sum(Adj, 2)); % Degree matrix
    L = D - Adj;           % Laplacian matrix
    
    % Create graph object and plot
    G = graph(Adj);
    % figure;
    % plot(G);
    % title(['Random Clique Graph with ', num2str(l), ' Cliques and ', num2str(n), ' Nodes']);
end

function cliqueSizes = generate_random_partition(total, parts)
    % Helper function to generate a random partition of 'total' into 'parts' integers
    % that sum up to 'total'
    % Example: generate_random_partition(10, 3) might return [3, 2, 5]
    
    % Generate 'parts - 1' random cut points, then sort them
    cuts = sort([0, randperm(total - 1, parts - 1), total]);
    
    % Calculate sizes based on the cut points
    cliqueSizes = diff(cuts);
end
