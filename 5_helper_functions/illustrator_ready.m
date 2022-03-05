function illustrator_ready(f)

%%%%% to make it easier to bring into illustrator

% set(f,'units','normalized','outerposition',[0 0 2/3 1])
set(f,'paperpositionmode','auto');

ax = 1:size(f.Children,1);
for iax = 1:size(ax,2)
    if isa(f.Children(iax),'matlab.graphics.axis.Axes')
        set(f.Children(iax),'FontSize',12)
        set(f.Children(iax),'FontName','Arial')    
        set(f.Children(iax),'FontWeight','bold')
    elseif isprop(f.Children(iax),'FontName')        
        set(f.Children(iax),'FontName','Arial')
        set(f.Children(iax),'FontWeight','bold')
    end
end



% set(f,'renderer','Painters')

end