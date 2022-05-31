function [] = surplot_pcolor(x, y, z, xl, yl, zl, tl)

    figure;
    h = pcolor(x, y, z); colorbar; colormap('jet');  shading flat; 
%     h = contour(x, y, z); colorbar; colormap('jet');  
    
    
    dt = datatip(h, 0, 0);
    
    h.ZData = h.CData;
    
    
%     % Generate an invisible datatip to ensure that DataTipTemplate is generated
%     dt = datatip(h,h.XData(1),h.YData(1),'Visible','off');
%     % Replace datatip row labeled Z with CData
%     idx = strcmpi('Z',{h.DataTipTemplate.DataTipRows.Label});
%     newrow = dataTipTextRow(,h.CData);
%     h.DataTipTemplate.DataTipRows(idx).Value = str2double(newrow.Value);
%     % Remove invisible datatip
%     delete(dt)
    
    
  
%     caxis([-7.5 7.5]);
    
%     axis equal tight; 

    xlim([min(min(x)) max(max(x))]);
    ylim([min(min(y)) max(max(y))]);

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
    box on;
%     title(tl, 'Interpreter', 'latex', 'FontSize', 18);
%     xlabel(xl, 'Interpreter', 'latex', 'FontSize', 18);
%     ylabel(yl, 'Interpreter', 'latex', 'FontSize', 18);
%     zlabel(zl, 'Interpreter', 'latex', 'FontSize', 18);
    title(tl,  'FontSize', 18);
    xlabel(xl, 'FontSize', 18);
    ylabel(yl, 'FontSize', 18);
    zlabel(zl, 'FontSize', 18);
    
end