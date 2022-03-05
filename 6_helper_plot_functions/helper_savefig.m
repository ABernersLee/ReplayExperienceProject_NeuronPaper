function helper_savefig(label)

% saveas(gcf,[label '.fig'],'fig')
% saveas(gcf,[label '.tif'],'tiff')
% saveas(gcf,[label '.eps'],'epsc')

saveas(gcf,[label '.fig'],'fig')
saveas(gcf,[label '.tif'],'tiff')


illustrator_ready(gcf)

% saveas(gcf,[label '.eps'],'epsc')
% printeps(1,[label '.eps'])

try
    export_fig(gcf,label,'-transparent','-eps','-painters') %,'-font_space', '', '-regexprep', 'Arial')
    
catch
    disp(['Problem, skipped'])
    
    set(gcf,'renderer','Painters')
    saveas(gcf,[label '.eps'],'epsc')
end

