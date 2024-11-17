function [gamma_opt,K_opt,P_opt] = Hinfty_centralized(params,flag)

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

    Q = sdpvar(n,n,'symmetric');
    Z = sdpvar(n,n,'full');
    gamma= sdpvar(1,1);
    
    % Objective and constraints
    
    LMI = [];
    
    LMI_1 = [-Q, A*Q+B*Z, Bw, zeros(n);
             (A*Q+B*Z)', -Q, zeros(n), (C*Q+D*Z)';
             Bw', zeros(n), -gamma * eye(n), Dw';
             zeros(n), C*Q + D*Z, Dw, -gamma * eye(n)];
    
    LMI = [LMI_1 <= 0, Q >= 0, gamma>=0];

    result = optimize(LMI,gamma,ops)
    %%  --------  end LMI  --------
    gamma_opt = value(gamma);
    if result.problem == 1
        gamma_opt = 0;
    end

    Z_opt = value(Z);
    Q_opt = value(Q);
    K_opt = Z_opt*inv(Q_opt);
    P_opt = inv(Q_opt);
    

    fprintf('-------------------------------------------\n');
    fprintf('----------- Centralized control -----------\n')
    fprintf(' gamma_opt                      : %8.3e \n', gamma_opt);
    fprintf(' min of Ps eigval               : %8.2e \n', min(eig(P_opt)));
    fprintf(' condition number of P          : %8.2e \n', max(eig(P_opt))/min(eig(P_opt)));
    fprintf(' Norm of K                      : %8.2e \n', norm(K_opt));
    fprintf(' max of A+BKs eigval (abs)      : %8.2e \n', max( abs(eig( A + B*K_opt )) ));
    fprintf('-------------------------------------------\n');

end