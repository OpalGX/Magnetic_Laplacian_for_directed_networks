close all
clear all
clc

addpath functions/
addpath tensor_toolbox/

tic
load('processed_data/corr')
load('processed_data/condP')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3-ring example

%network
[G, A] = ring(3, 10, 0.5, 0.5);

p = imagesc(adjacency(G)); %plot color map of original matrix
colormap(parula);
colorbar;
ax = gca;
ax.FontSize = 16;
saveas(p,'Plot/3_ring_adjacency.jpg');

x = cat(2,randperm(10),randperm(10)+15, randperm(10)+30);
y = cat(2,randperm(10),randperm(10)+30, randperm(10));
p = plot(G,'XData', x ,'YData',y, 'MarkerSize', 15, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5], 'NodeLabel',{});
highlight(p,1:10,'NodeColor',[0 0.4470 0.7410]);
highlight(p,11:20,'NodeColor',[0.8500 0.3250 0.0980]);
highlight(p,21:30,'NodeColor',[0.9290 0.6940 0.1250]);
xticks([]);
yticks([]);
ax = gca;
ax.FontSize = 16;
saveas(p,'Plot/3_ring_network.jpg');

%original adjacency matrix

%magnetic laplacian
[G, D, p, V, D] = meigenmaps(G,1/3);
p = plot(G, 'XData', cos(G.Nodes.phase0) ,'YData', sin(G.Nodes.phase0), 'MarkerSize', 15, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5], 'NodeLabel',{});
highlight(p,1:10,'NodeColor',[0 0.4470 0.7410]);
highlight(p,11:20,'NodeColor',[0.8500 0.3250 0.0980]);
highlight(p,21:30,'NodeColor',[0.9290 0.6940 0.1250]);
ax = gca;
ax.XAxisLocation = 'origin';
ax.YAxisLocation = 'origin';
xticks([-1  0  1])
%xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})
yticks([-1 0 1])
%xticklabels({'-3\pi','-2\pi','-\pi','0','\pi','2\pi','3\pi'})

ax.FontSize = 16; 
xlabel('$\cos{\theta}$','Interpreter','latex', 'FontSize',18);
ylabel('$\sin{\theta}$','Interpreter','latex', 'FontSize', 18);
saveas(p,'Plot/3_ring_meigenmap.jpg');

%normal laplacian
W = adjacency(G);
Ws = (W+W.')/2;
deg = outdegree(G); D = diag(deg); % Degree matrix.?
L = D - Ws; % Magnetic Laplacian.
Ln = inv(sqrtm(D))*L*inv(sqrtm(D));% Normalized Laplacian.
[V,D] = eig(Ln); % Eigenvectors.

G.Nodes.Eig0 = V(:,1);
G.Nodes.Eig1 = V(:,2);

p = plot(G,'XData',G.Nodes.Eig0,'YData', G.Nodes.Eig1, 'MarkerSize', 15, 'LineWidth', 1, 'EdgeColor', [0.5 0.5 0.5], 'NodeLabel',{});
highlight(p,1:10,'NodeColor',[0 0.4470 0.7410]);
highlight(p,11:20,'NodeColor',[0.8500 0.3250 0.0980]);
highlight(p,21:30,'NodeColor',[0.9290 0.6940 0.1250]);
xlabel('Eig_1', 'FontSize',18);
ylabel('Eig_2', 'FontSize',18);
ax = gca;
ax.FontSize = 16; 
saveas(p,'Plot/3_ring_norm_Laplacian.jpg');

