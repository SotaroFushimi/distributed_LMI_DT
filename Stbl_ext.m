function [K_opt,P_opt,eig_max,X_opt] = Stbl_ext(params,flag)

    n = params.n;
    G = params.G;

    A = params.A;
    B = params.B;
    C = params.C;
    D = params.D;

    Bw = params.Bw;
    Dw = params.Dw; 

    %%  -------- solve LMI -------- 
    ops = sdpsettings('solver',params.solver_chosens);
    ops.verbose = flag;

    Z = [];
    Q = sdpvar(n);
    X = diag(sdpvar(n,1));
    % gamma = 1.0981;
    
    L = laplacian(G);
    for i = 1:n
        tmp = [];
        for j = 1:n
            if L(i,j) ~= 0
                tmp = [tmp sdpvar(1,1)];
            else
                tmp = [tmp 0];
            end
        end
        Z = [Z;tmp];
    end
    
    % Objective and constraints
    
    LMI = [];

    LMI_1 = [X+X'-Q,(A*X+B*Z)';
             A*X+B*Z,Q];
    
    LMI = [LMI_1 >= 0, Q >= 0];

    optimize(LMI,[],ops)
    %%  --------  end LMI  --------


    %% results 
    Z_opt = value(Z);
    Q_opt = value(Q);
    X_opt = value(X);
    K_opt = Z_opt*inv(X_opt);
    P_opt = Q_opt;
    eig_max=max( abs(eig( A + B*K_opt )));

    fprintf('-------------------------------------------\n');
    fprintf('--------- Stabilizer Extended LMI ---------\n')
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' min of Xs eigval               : %8.2e \n', min(eig(X_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt)); 
    fprintf(' max of A+BKs eigval (abs)      : %8.2e \n', eig_max);
    fprintf('-------------------------------------------\n');

end