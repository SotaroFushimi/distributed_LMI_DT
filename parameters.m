%% define parameters
%%%%%%%%% --------------- start ---------------
params.n = 12; %% the number of nodes 偶数とする
params.l = 3;
params.G = generate_cliquegraph(params.n,params.l);

% define A
eig_A_theta     = 2*pi* rand(params.n/2,1);
eig_A_theta     = [eig_A_theta;-eig_A_theta];
eig_A_r         = 1+5*rand(params.n/2,1);
eig_A_r         = [eig_A_r;eig_A_r];
eig_A           = complex(eig_A_r.*cos(eig_A_theta), eig_A_r.*sin(eig_A_theta));
tmp = rand(params.n,params.n);
params.A        = tmp*diag(eig_A)\tmp;

% define other system matrices
params.B = rand(params.n,params.n);
params.C = eye(params.n);
params.D = eye(params.n);
params.Bw = eye(params.n);
params.Dw = eye(params.n);

%%%%%%%%% --------------  end  ---------------
fprintf('Rank deficiency of the contrability matrix - n:%8.2e \n', rank(ctrb(params.A,params.B))-params.n);


% params.solver_chosens = 'sedumi';
params.solver_chosens = 'sdpt3';