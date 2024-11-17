clear all;
num_sim = 30;

gamma_result_hist = zeros(num_sim,5);
stab_result_hist = zeros(num_sim,5);
sim_results(num_sim) = struct();

for iii = 1:num_sim
    iii
    parameters;
    Hinfty_test;
    gamma_result_hist(iii,:) = gamma_result;
    stab_result_hist(iii,:) = stab_result;

    sim_results(iii).gamma = gamma_result;
    sim_results(iii).stab = stab_result;
    sim_results(iii).A = params.A;
    sim_results(iii).B = params.B;
    sim_results(iii).C = params.C;
    sim_results(iii).D = params.D;
    sim_results(iii).Bw = params.Bw;
    sim_results(iii).Dw = params.Dw;
    sim_results(iii).G = params.G;
    sim_results(iii).P_diag = P_opt_diag;
    sim_results(iii).K_diag = K_opt_diag;
    sim_results(iii).P_ext = P_opt_ext;
    sim_results(iii).K_ext = K_opt_ext;
    sim_results(iii).X_ext = X_opt_ext;
    sim_results(iii).P_S = P_opt_S;
    sim_results(iii).K_S = K_opt_S;
    sim_results(iii).P_proposed = P_opt_proposed;
    sim_results(iii).K_proposed = K_opt_proposed;
    sim_results(iii).X_proposed = X_opt_proposed;
    sim_results(iii).P_cen = P_opt_cen;
    sim_results(iii).K_cen = K_opt_cen;
end

%% change failed gamma values for plotting
for ll = 1:num_sim
    for lll = 1:4
        if stab_result_hist(ll,lll)>=1
            gamma_result_hist(ll,lll) = 10^6;
        else
            gamma_result_hist(ll,lll) = gamma_result_hist(ll,lll)/gamma_result_hist(ll,5);
        end
        if gamma_result_hist(ll,lll)<0.999
            gamma_result_hist(ll,lll) = 10^6;
        end
    end
end

%% 
clf;
gamma_result_hist(gamma_result_hist>10^5) = 10^5;
ms = 10;
ran = 1:num_sim;
fig = figure;
semilogy(ran,gamma_result_hist(:,1),'LineStyle','none','MarkerSize',20,'Marker','pentagram','MarkerFaceColor','k');
hold on

semilogy(ran,gamma_result_hist(:,2),'LineStyle','none','MarkerSize',15,'Marker','square','MarkerFaceColor','g');
semilogy(ran,gamma_result_hist(:,3),'LineStyle','none','MarkerSize',15,'Marker','square','MarkerFaceColor','b');
semilogy(ran,gamma_result_hist(:,4),'LineStyle','none','MarkerSize',13,'Marker','o','MarkerFaceColor','r');

hold off
ax = gca;
ax.XAxis.FontSize = 14;
ax.YAxis.FontSize = 14;

xlabel('Sample number','FontSize',18)
ylabel('$\gamma_*/\gamma_\mathrm{*,cen}$','Interpreter', 'latex','FontSize',22)

xlim([0 num_sim+1])
ylim([1 3])

legend('block-diagonal','extended LMI','clique-wise','proposed', 'FontSize', 14)

saveas(fig, 'data/Hinfty/Hinfty_result_new.png')
save('data/Hinfty/Hinfty_result_new.mat', 'sim_results');
