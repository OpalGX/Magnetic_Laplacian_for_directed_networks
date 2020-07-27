close all
clear all
clc

addpath functions/
addpath tensor_toolbox/

tic

load('datasets/melted_cormat')

G0 = digraph(M, names); 

G = max_connected_subgraph(G0);
A = adjacency(G);
p = plot(G,'Layout','force');
saveas(p,'Plot/dominance_network.fig');

p = imagesc(A); 
colormap(parula);
colorbar;
ax = gca;
ax.FontSize = 14;
set(gca,'Ytick',1:length(names),'YTickLabel',1:length(names))
saveas(p,'Plot/dominance_colormap.jpg');


G = max_connected_subgraph(G0);
[G, D, p, V, D] = meigenmaps(G,0.1);
[n, ~]=size(G.Nodes.phase0);
p = plot(G, 'XData', 0.1*(rand([n 1])-0.5)+G.Nodes.phase0 ,'YData',1*(rand([n 1])-0.5)+G.Nodes.phase1,'NodeColor',[1 0.5 0.1250], 'MarkerSize', 15, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5], 'NodeFontSize',12,'ArrowSize',10);
%p = plot(G, 'XData', 0.4*(rand([n 1])-0.5)+cos(G.Nodes.phase0) ,'YData',0.4*(rand([n 1])-0.5)+sin(G.Nodes.phase0),'NodeColor',[1 0.5 0.1250], 'MarkerSize', 10, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5], 'NodeFontSize',10,'ArrowSize',8);

axis([-0.2 2 -1 7]);
ax = gca;
ax.FontSize = 12;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xlabel('phase(\phi_0)', 'Fontsize',16);
ylabel('phase(\phi_1)', 'Fontsize',16);
saveas(p,'Plot/dominance_network_eigenmaps.jpg');

[out,idx] = sort(G.Nodes.phase0);
A_reordered = A(idx,idx);
p = imagesc(A_reordered); %plot color map of matrix re-ordered by normal Laplacian
colormap(parula);
colorbar;
ax = gca;
ax.FontSize = 16;
saveas(p,'Plot/dominance_colormap_reordered.jpg');


toc