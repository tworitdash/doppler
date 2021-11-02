function [] = plott(x, y, xl, yl, tl, lw)

    figure;
    plot(x, y, 'LineWidth', lw); 
    grid on;
    
%     axis equal tight; 

    xlim([min(x) max(x)]);
    ylim([min(y) max(y)]);

%     T = [-x(end):5:x(end)];
%     L = {};
%     for m = 1:length(T)
%         if T(m) == 0
%             L{m} = '0';
%         else
%             L{m} = num2str(T(m), '%3.1f');
%         end
%     end

%     set(gca, 'XTick', T, 'XTickLabel', L);
%     set(gca, 'YTick', T, 'YTickLabel', L);

    set(gca, 'FontSize', 18, 'LineWidth', 4);

%     title(tl, 'Interpreter', 'latex', 'FontSize', 18);
%     xlabel(xl, 'Interpreter', 'latex', 'FontSize', 18);
%     ylabel(yl, 'Interpreter', 'latex', 'FontSize', 18);
%     zlabel(zl, 'Interpreter', 'latex', 'FontSize', 18);
    title(tl,  'FontSize', 18);
    xlabel(xl, 'FontSize', 18);
    ylabel(yl, 'FontSize', 18);
    
end