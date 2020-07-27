close all
clear all
clc

addpath functions/
addpath tensor_toolbox/

tic
n = 10; %nodes at each cluster
g_array = linspace(0, 1, 13);
k_array = [2 3 4 5 6];
p_out = 1;
p_in = 0;
objective = zeros(length(k_array), length(g_array));
bcubed = zeros(length(k_array), length(g_array));
precision = zeros(length(k_array), length(g_array));
recall = zeros(length(k_array), length(g_array));
bc = zeros(max(k_array),1);
rc = zeros(max(k_array), 1);
pr = zeros(max(k_array), 1);

i = 1;
for k = k_array
    j = 1;
    [G, A] = ring(k, n, p_in, p_out);
    G = max_connected_subgraph(G);
    %p = plot(G, 'layout', 'force', 'MarkerSize', 15, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5],'NodeLabel',G.Nodes.L);  
    %title(strcat(num2str(k),' cycle'));
    %saveas(p,strcat('Plot/',num2str(k),'_ring.jpg'));

    %img = imagesc(A); %plot color map of original matrix
    %colormap(parula);
    %colorbar;
    %title(strcat(num2str(k),' ring adjacency matrix'));
    %saveas(img,strcat('Plot/',num2str(k),'_ring_colormap_b.jpg'));
    
    for g = g_array
        [G, D, p] = meigenmaps(G,g);
        %p = plot(G, 'XData', cos(G.Nodes.phase0) ,'YData', sin(G.Nodes.phase0), 'MarkerSize', 15, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5], 'NodeLabel',{});
        %title(strcat(num2str(k), ' cycle, g= ', num2str(g,2)));
        %ax = gca;
        %ax.XAxisLocation = 'origin';
        %ax.YAxisLocation = 'origin';
        %xticks([-1  0  1])
        %xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
        %yticks([-1 0 1])
        %xlabel('$\sin{\theta}$','Interpreter','latex');
        %ylabel('$\cos{\theta}$','Interpreter','latex');
        %saveas(p,strcat('Plot/',num2str(k),'_ring_meigenmap_g=_',num2str(g, 2),'.jpg'));
            
        %colormap
        [out,idx] = sort(G.Nodes.phase0);
        M = adjacency(G);
        M_ordered = M(idx,idx);
        %img = imagesc(M_ordered); %plot color map of original matrix
        %colormap(parula);
        %colorbar;
        %title(strcat(num2str(k),' ring re-ordered g = ', num2str(g,2)));
        %saveas(img,strcat('Plot/',num2str(k),'_ring_colormap_a_g=_',num2str(g, 2),'.fig'));
        
        
        %k-mean cluster
        X = [cos(G.Nodes.phase0), sin(G.Nodes.phase0)];
        n_cluster = k;
        [cluster,centroid] = kmeans(X, n_cluster, 'MaxIter',1000, 'Replicates',1000, 'Start', 'uniform');
        G.Nodes.C = cluster;    
        %calculate precision and recall
        [BCubed,f_precision,f_recall] = Calculate_Cluster_BCubed_precision(G.Nodes.C,G.Nodes.L);


      
        objective(i,j) = D;
        precision(i,j) = f_precision;

        %precision(i,j) = calculate_clustering_precision(M_ordered, k, n);
        
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
        

   
Fig = figure
x = g_array;
p = plot(x,real(objective(1,:)),x,real(objective(2,:)),x,real(objective(3,:)),x,real(objective(4,:)),x,real(objective(5,:)));
ylim([0 1.2]);
xlabel('g')
ylabel('Objective function')
title(strcat('objective, po=', num2str(p_out,2), ', pi=', num2str(p_in,2)));
legend('2 cycle','3 cycle','4 cycle','5 cycle', '6 cycle')
%set(gca,'xtick',g_array,'xticklabel', num2str(g_array))
saveas(Fig,'Plot/Optimal_g_k_ring.fig');


%Fig = figure
%x = g_array;
%p = plot(x,bcubed(1,:),x,bcubed(2,:),x,bcubed(3,:),x,bcubed(4,:),x,bcubed(5,:));
%xlabel('g')
%ylabel('Bcubed')
%title(strcat('bcubed, po=', num2str(p_out,2), ', pi=', num2str(p_in,2)));
%legend('2 cycle','3 cycle','4 cycle','5 cycle', '6 cycle')
%saveas(Fig,'Plot/bcubed_g_k_ring.fig');

%Fig = figure
%x = g_array;
%p = plot(x,recall(1,:),x,recall(2,:),x,recall(3,:),x,recall(4,:),x,recall(5,:));
%xlabel('g')
%ylabel('Recall')
%title(strcat('recall, po=', num2str(p_out,2), ', pi=', num2str(p_in,2)));
%legend('2 cycle','3 cycle','4 cycle','5 cycle', '6 cycle')
%saveas(Fig,'Plot/recall_g_k_ring.jpg');

Fig = figure
x = g_array;
p = plot(x,precision(1,:), '-.o',x,precision(2,:), '--^',x,precision(3,:), ...
    ':s',x,precision(4,:), '-.p',x,precision(5,:), '--h', 'LineWidth',2, 'MarkerSize',10);
ylim([0 1.2]);
xlabel('g')
ylabel('Precision')
%title(strcat('precision, po=', num2str(p_out,2), ', pi=', num2str(p_in,2)));
legend('2 cycle','3 cycle','4 cycle','5 cycle', '6 cycle')
saveas(Fig,'Plot/precision_g_k_ring.jpg');

%save('processed_data/precision','precision')

toc 


