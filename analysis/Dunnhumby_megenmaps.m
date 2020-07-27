close all
clear all
clc

addpath functions/
addpath tensor_toolbox/

tic
load('processed_data/corr')
load('processed_data/condP')


%visualize shopping data
n = length(condP);
W = condP - eye(n);
G0 = digraph(W);
G0.Nodes.Name = commodity_name;

%visualize the eigenmap of largest connected component
G = max_connected_subgraph(G0);
imagesc(adjacency(G)); %plot color map of original matrix
colormap(parula);
colorbar;

[G, D, p] = meigenmaps(G,0.25);

title('Dunnhumby data meigenmap');
saveas(p,'Plot/Dunnhumby data meigenmap.fig');


%visualize the eigenmap of normalized Laplacian connected component
[G, p] = normLaplacian(G0);
title('Dunnhumby norm Laplacian');
saveas(p,'Plot/Dunnhumby_norm_Laplacian.fig');
[out,idx] = sort(G.Nodes.Eig2);
M = adjacency(G);
M_nLap = M(idx,idx);
imagesc(M_nLap); %plot color map of matrix re-ordered by normal Laplacian
colormap(parula);
colorbar;


threshold = 0.5;
g_array = 0.25;
objective = zeros(length(threshold), length(g_array));
col = 1;l
for g = g_array
    row = 1;
    for i = 1:length(threshold)
        W = (condP > threshold(i))-eye(n);
        G0 = digraph(W);
        G0.Nodes.Name = commodity_name;
        P = plot(G0,'NodeLabel',G0.Nodes.Name);
        %colormap
        img = imagesc(W); %plot color map of original matrix
        colormap(parula);
        colorbar;
        saveas(img,strcat('Plot/Dunnhumby_colormap_b_',num2str(threshold(i)),'.fig'));

        %visualize the eigenmap of largest connected component
        G = max_connected_subgraph(G0);
        [G, D, p] = meigenmaps(G,g);
        [n_nodes, ~]=size(G.Nodes.phase0);
        p = plot(G, 'XData', 0.3*(rand([n_nodes 1])-0.5)+cos(G.Nodes.phase0) ,'YData',0.2*(rand([n_nodes 1])-0.5)+sin(G.Nodes.phase0),'NodeColor',[1 0.5 0.1250], 'MarkerSize', 10, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5], 'NodeFontSize',10,'ArrowSize',10);
        ax = gca;
        ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';
        axis([-0.6 2 -1.4 0.3]);
        xticks([-0.5  0  1])
        %xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
        %yticks([-1 0 1])
        ax.FontSize = 16;
        xlabel('$\cos{\theta}$','Interpreter','latex', 'Fontsize', 18);
        ylabel('$\sin{\theta}$','Interpreter','latex', 'Fontsize', 18);
        %title(strcat('Dunnhumby meigenmap threshold =  ',num2str(threshold(i)), ' g= ', num2str(g,2)));
        saveas(p,strcat('Plot/Dunnhumby_meigenmap_',num2str(threshold(i)),'_g=_',num2str(g, 2),'.jpg'));

        %colormap
        [out,idx] = sort(G.Nodes.phase0);
        M = adjacency(G);
        M_nLap = M(idx,idx);
        img = imagesc(M_nLap); %plot color map of original matrix
        colormap(parula);
        colorbar;
        saveas(img,strcat('Plot/Dunnhumby_colormap_a_',num2str(threshold(i)),'_g=_',num2str(g, 2),'.fig'));
        objective(row,col) = D;
        row = row+1;
    end
    col = col+1;
end


x = g_array;
p = plot(x,real(objective(1,:)));
title('Smallest Eigenvalue')
xlabel('g')
ylabel('Eigenvalue')
legend('conP > 0.4')
saveas(p,strcat('Plot/Dunnhumby_eigenvalues.fig'));
