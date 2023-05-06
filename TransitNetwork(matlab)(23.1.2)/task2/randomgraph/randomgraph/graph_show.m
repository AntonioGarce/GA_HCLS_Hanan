%%
clear
close all 
clc

%% Load matricies 
load('KRAKOW_data_final.mat')
weight = load('cracow.mat').cracow;
weight=lt(weight,1000);    %make weight matrix to detrmine the links between nodes

%% Graph Construction 
G=graph(weight,'omitselfloops');
hh= plot((G),'Markersize',4);
disp(height(G.Nodes))
numnodes = height(G.Nodes);
nodeColor = zeros(numnodes,3);
nodeMarker = zeros(1,numnodes);
for k = 1:numnodes
    switch(G.degree(k))
        case 1
            nodeColor(k,:) = [1 0 0];
        case 2
            nodeColor(k,:)  = [0 1 0];
        case 3
            nodeColor(k,:)  = [0 0 1];
        otherwise
            nodeColor(k,:)  = [1 1 0];
    end
end
hh.NodeColor = nodeColor;
hh.NodeLabelMode = 'auto';
hh.Marker = 'o';
hh.MarkerSize = 6;

for k=1:height(G.Edges)
    
end



