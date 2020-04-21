
for i = 1:length(g_num)
    
    tag = g_num([3,5,1,7],i);
    x = g_coord(1,tag);
    y = g_coord(2,tag);
    X = [x(1),x(2),x(4),x(3)];
    Y = [y(1),y(2),y(4),y(3)];
    
    fill(X,Y,'w','LineWidth',1)
    text((x(1)+x(4))/2,(y(1)+y(4))/2,num2str(i),'Color','red','FontSize',14)
    hold on
end

%
for i = 1:length(g_coord)
    coord = g_coord(:,i);
    text(coord(1),coord(2),num2str(i),'Color','blue','FontSize',14)
    hold on
end
%}

set(gca,'xtick',[],'xticklabel',[])
set(gca,'ytick',[],'yticklabel',[])