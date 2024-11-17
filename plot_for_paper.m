% Load the first dataset
load('data/Hinfty/Hinfty_result_N10l3.mat', 'sim_results');
sim_results_1 = sim_results;

% Load the second dataset
load('data/Hinfty/Hinfty_result_N40l5.mat', 'sim_results');
sim_results_2 = sim_results;

num_sim = 30;

% Initialize arrays to hold gamma and stability results for both datasets
gamma_result_hist_1 = zeros(num_sim, 5);
stab_result_hist_1 = zeros(num_sim, 5);
gamma_result_hist_2 = zeros(num_sim, 5);
stab_result_hist_2 = zeros(num_sim, 5);

% Populate the arrays for the datasets
for iii = 1:num_sim
    gamma_result_hist_1(iii, :) = sim_results_1(iii).gamma;
    stab_result_hist_1(iii, :) = sim_results_1(iii).stab;
    gamma_result_hist_2(iii, :) = sim_results_2(iii).gamma;
    stab_result_hist_2(iii, :) = sim_results_2(iii).stab;
end

% Adjust gamma values for failed stability cases in the first dataset
for ll = 1:num_sim
    for lll = 1:4
        if stab_result_hist_1(ll, lll) >= 1
            gamma_result_hist_1(ll, lll) = 10^6;
        else
            gamma_result_hist_1(ll, lll) = gamma_result_hist_1(ll, lll) / gamma_result_hist_1(ll, 5);
        end
        if gamma_result_hist_1(ll, lll) < 0.999
            gamma_result_hist_1(ll, lll) = 10^6;
        end
    end
end

% Adjust gamma values for failed stability cases in the second dataset
for ll = 1:num_sim
    for lll = 1:4
        if stab_result_hist_2(ll, lll) >= 1
            gamma_result_hist_2(ll, lll) = 10^6;
        else
            gamma_result_hist_2(ll, lll) = gamma_result_hist_2(ll, lll) / gamma_result_hist_2(ll, 5);
        end
        if gamma_result_hist_2(ll, lll) < 0.999
            gamma_result_hist_2(ll, lll) = 10^6;
        end
    end
end

% Cap large gamma values for better visualization
gamma_result_hist_1(gamma_result_hist_1 > 10^5) = 10^5;
gamma_result_hist_2(gamma_result_hist_2 > 10^5) = 10^5;

% Plotting
fig = figure;
width = 800;
height = 450;
set(fig, 'Position', [100, 100, width, height]); % [左端, 下端, 横幅, 縦幅]

% Plot first dataset
ran_1 = 1:num_sim;
semilogy(ran_1, gamma_result_hist_1(:, 1), 'LineStyle', 'none', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', 'k');
hold on;
semilogy(ran_1, gamma_result_hist_1(:, 2), 'LineStyle', 'none', 'MarkerSize', 15, 'Marker', 'square', 'MarkerFaceColor', 'g');
semilogy(ran_1, gamma_result_hist_1(:, 3), 'LineStyle', 'none', 'MarkerSize', 15, 'Marker', 'square', 'MarkerFaceColor', 'b');
semilogy(ran_1, gamma_result_hist_1(:, 4), 'LineStyle', 'none', 'MarkerSize', 13, 'Marker', 'o', 'MarkerFaceColor', 'r');

% Offset for second dataset
offset = num_sim + 1;
ran_2 = (1:num_sim) + offset;

% Plot second dataset
semilogy(ran_2, gamma_result_hist_2(:, 1), 'LineStyle', 'none', 'MarkerSize', 16, 'Marker', 'pentagram', 'MarkerFaceColor', 'k');
semilogy(ran_2, gamma_result_hist_2(:, 2), 'LineStyle', 'none', 'MarkerSize', 15, 'Marker', 'square', 'MarkerFaceColor', 'g');
semilogy(ran_2, gamma_result_hist_2(:, 3), 'LineStyle', 'none', 'MarkerSize', 15, 'Marker', 'square', 'MarkerFaceColor', 'b');
semilogy(ran_2, gamma_result_hist_2(:, 4), 'LineStyle', 'none', 'MarkerSize', 13, 'Marker', 'o', 'MarkerFaceColor', 'r');

% Add separation line between datasets
sep_line_x = (num_sim + offset + 1) / 2;
plot([sep_line_x, sep_line_x], [1, 7], '--k', 'LineWidth', 1.5);
hold off;

% Customize axes and labels
ax = gca;
ax.XAxis.FontSize = 16;
ax.YAxis.FontSize = 16;

xlabel('Sample number', 'FontSize', 20);
ylabel('$\gamma_*/\gamma_\mathrm{cen}$', 'Interpreter', 'latex', 'FontSize', 30);

% Set axis limits
xlim([0 num_sim + offset + 1]);
ylim([1 2.7]);
xticks([1 10 20 30 32 41 51 61]);
xticklabels({"1", "10", "20", "30", "1", "10", "20", "30"});
% yticks([1 1.1 1.2 1.3 1.4 1.5 1.6 2.7]);
% yticklabels({"1", "1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "2.7"});

% Add legend
legend('Block-diagonal', 'Extended LMI', 'Clique-wise', 'Proposed', 'Location', 'northwest','FontSize', 16);

% Add dataset labels
text(num_sim / 2, 2.8, '(N, l) = (10, 3)', 'FontSize', 18, 'HorizontalAlignment', 'center');
text(offset + num_sim / 2, 2.8, '(N, l) = (40, 5)', 'FontSize', 18, 'HorizontalAlignment', 'center');

% Save the figure as PNG
saveas(fig, 'data/Hinfty/Combined_Hinfty_N10l3_N40l5_final.png');
