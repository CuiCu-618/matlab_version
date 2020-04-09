nels = size(g_num,2);
% faces with tip 1 and 8
F = [1,2,3,4;1,5,8,4;1,2,6,5;7,6,5,8;7,3,4,8;7,6,2,3];

for iel = 1:nels
    num = g_num(:,iel);
    coord = g_coord(:,num);
    patch('Faces',F,'Vertices',coord','FaceColor','none','LineWidth',1,'EdgeColor','k');
    hold on
end
xlabel('x')
ylabel('y')
zlabel('z')
title("undeformed")

figure
for iel = 1:nels
    num = g_num(:,iel);
    coord = g_coord(:,num);
    coord = coord + loads(nf(:,num));
    patch('Faces',F,'Vertices',coord','FaceColor','none','LineWidth',1,'EdgeColor','k');
    hold on
end
xlabel('x')
ylabel('y')
zlabel('z')
title("deformed")