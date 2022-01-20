function [] = plott2(x, y, xl, yl, tl, lw, dtext, c, m)
    plot(x, y, m, 'Color', c, 'LineWidth', lw, 'DisplayName', dtext); 
    grid on;
    
%     axis equal tight; 

    xlim([min(x) max(x)]);
%     ylim([min(y) max(y)]);

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

    set(gca, 'FontSize', 20, 'LineWidth', 4);
    box on;
    

%     title(tl, 'Interpreter', 'latex', 'FontSize', 18);
%     xlabel(xl, 'Interpreter', 'latex', 'FontSize', 18);
%     ylabel(yl, 'Interpreter', 'latex', 'FontSize', 18);
%     zlabel(zl, 'Interpreter', 'latex', 'FontSize', 18);
    title(tl,  'FontSize', 20);
    xlabel(xl, 'FontSize', 20);
    ylabel(yl, 'FontSize', 20);
    lgd = legend;
    lgd.FontSize = 16;
    lgd.FontWeight = 'bold';
    % for x
%     set(gca,'XTickLabel',[])
%     set(gca,'XTick',[])
    set(gca,'xcolor', c) 
    % for y
%     set(gca,'YTickLabel',[])
%     set(gca,'YTick',[])
    set(gca,'ycolor', c)
    
    
    
end