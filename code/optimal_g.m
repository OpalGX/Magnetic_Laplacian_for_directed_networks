close all
clear all
clc

addpath functions/
addpath tensor_toolbox/

tic
n = 10; %nodes at each cluster
g_array = [1/6 1/3 1/2 2/3 5/6];
k_array = [6]
p_in = 0;
p_out = 1;
objective = zeros(length(k_array), length(g_array));
i = 1;
for k = k_array
    j = 1;
    [G, A] = ring(k, n, p_in, p_out);
    p = plot(G, 'layout', 'force', 'MarkerSize', 15, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5],'NodeLabel',{});
    title(strcat(num2str(k),' cycle'));
    saveas(p,strcat('Plot/',num2str(k),'_ring.jpg'));

    img = imagesc(A); %plot color map of original matrix
    colormap(parula);
    colorbar;
    title(strcat(num2str(k),' ring adjacency matrix'));
    saveas(img,strcat('Plot/',num2str(k),'_ring_colormap_b.jpg'));

    for g = g_array
        [G, D, p] = meigenmaps(G,g);
        p = plot(G, 'XData', cos(G.Nodes.phase0) ,'YData',sin(G.Nodes.phase0), 'MarkerSize', 15, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5], 'NodeLabel',{});
        highlight(p,1:10,'NodeColor',[0 0.4470 0.7410]);
        title(strcat(num2str(k), ' cycle, g= ', num2str(g,2)));
        xlabel('phase(\phi_0)');
        ylabel('phase(\phi_1)');
        saveas(p,strcat('Plot/',num2str(k),'_ring_meigenmap_g=_',num2str(g, 2),'.jpg'));
            
        %colormap
        [out,idx] = sort(G.Nodes.phase0);
        M = adjacency(G);
        M_ordered = M(idx,idx);
        img = imagesc(M_ordered); %plot color map of original matrix
        colormap(parula);
        colorbar;
        title(strcat(num2str(k),' ring re-ordered g = ', num2str(g,2)));
        saveas(img,strcat('Plot/',num2str(k),'_ring_colormap_a_g=_',num2str(g, 2),'.fig'));
        
        objective(i,j) = D;
        j = j+1;
    end
    i = i+1;
end


img = imagesc(real(objective)); %plot color map of original matrix
colorbar;
title('Objective function for k-cycle');
set(gca,'Ytick',1:4,'YTickLabel',{'3'; '4'; '5'; '6'})
set(gca,'Xtick',1:6,'XTickLabel',{'1/6'; '1/5'; '1/4'; '1/3'; '1/2'; '1'})
xlabel('g');
ylabel('k');
saveas(img,strcat('Plot/objective_function.jpg'));
        

   
figure
x = g_array;

p = plot(x,real(objective(1,:)),x,real(objective(2,:)),x,real(objective(3,:)),x,real(objective(4,:)));
xlabel('g')
ylabel('Objective function')
legend('3 ring','4 ring','5 ring','6 ring')
names = {'CRHS'; 'ELLY'; 'LGWD'; 'ECFS'; 'THMS'};
%set(gca,'xtick',g_array,'xticklabel', num2str(g_array))
saveas(p,'Plot/Optimal_g_k_ring.fig');

toc 


