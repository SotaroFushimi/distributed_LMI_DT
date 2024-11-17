% Load the saved data
load('data/Hinfty/Hinfty_result_N10l3.mat', 'sim_results');

% Extract the number of simulations
num_sim = length(sim_results);

% Initialize arrays to hold gamma and stability results
gamma_result_hist = zeros(num_sim, 5);
stab_result_hist = zeros(num_sim, 5);

% Populate the arrays from the saved data
for iii = 1:num_sim
    gamma_result_hist(iii, :) = sim_results(iii).gamma;
    stab_result_hist(iii, :) = sim_results(iii).stab;
end

% Adjust gamma values for failed stability cases
for ll = 1:num_sim
    for lll = 1:4
        if stab_result_hist(ll, lll) >= 1
            gamma_result_hist(ll, lll) = 10^6;
        else
            gamma_result_hist(ll, lll) = gamma_result_hist(ll, lll) / gamma_result_hist(ll, 5);
        end
        if gamma_result_hist(ll, lll) < 0.999
            gamma_result_hist(ll, lll) = 10^6;
        end
    end
end

% Plotting
clf;
gamma_result_hist(gamma_result_hist > 10^5) = 10^5; % Cap large values for visibility
ran = 1:num_sim;
fig = figure;
ms = 10;

% Plot each method's gamma results with different markers and colors
semilogy(ran, gamma_result_hist(:, 1), 'LineStyle', 'none', 'MarkerSize', 20, 'Marker', 'pentagram', 'MarkerFaceColor', 'k');
hold on
semilogy(ran, gamma_result_hist(:, 2), 'LineStyle', 'none', 'MarkerSize', 15, 'Marker', 'square', 'MarkerFaceColor', 'g');
semilogy(ran, gamma_result_hist(:, 3), 'LineStyle', 'none', 'MarkerSize', 15, 'Marker', 'square', 'MarkerFaceColor', 'b');
semilogy(ran, gamma_result_hist(:, 4), 'LineStyle', 'none', 'MarkerSize', 13, 'Marker', 'o', 'MarkerFaceColor', 'r');
hold off

% Customize axes and labels
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;

xlabel('Sample number', 'FontSize', 20)
ylabel('$\gamma_*/\gamma_\mathrm{cen}$', 'Interpreter', 'latex', 'FontSize', 30)

xlim([0 num_sim + 1])
ylim([1 7])

legend('block-diagonal', 'extended LMI', 'clique-wise', 'proposed', 'FontSize', 14)

% Save the recreated figure
% saveas(fig, 'data/Hinfty/Hinfty_script_N10l3_2.png');
