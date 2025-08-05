
    % Time vector for plotting (from 0 to m*T)
    t_plot = (0:m) * T;

    % Plot y1, y2, y3, y4 over time
    figure;
    plot(t_plot, y1, 'LineWidth', 1.5); hold on;
    plot(t_plot, y2, 'LineWidth', 1.5);
    plot(t_plot, y3, 'LineWidth', 1.5);
    plot(t_plot, y4, 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Output y_i');
    legend('y1','y2','y3','y4');
    title('Outputs y1, y2, y3, y4 over Time');
    grid on;