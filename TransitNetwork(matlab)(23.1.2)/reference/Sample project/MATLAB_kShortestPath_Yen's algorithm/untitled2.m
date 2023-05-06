clc,close all,clear 

totalroutes = 30;
minroutelen = 12;
maxroutelen = 25;
%% read excel file
[num1,txt1,~] = xlsread('KRAKOW - FULL DATA.xlsx','TRAFFIC-DEMAND');
[num2,txt2,~] = xlsread('KRAKOW - FULL DATA.xlsx','GRAPH');
[num3,txt3,~] = xlsread('KRAKOW - FULL DATA.xlsx','CurrentSolution');

srcnodes = num1(:,1);
dstnodes = num1(:,2);
demands = num1(:,3);
%% get nodes and make matrix
nodes_demand = unique([srcnodes;dstnodes]);
numnodes_demand = length(nodes_demand);

srcnodes = num2(:,1);
dstnodes = num2(:,2);
demands = num2(:,3);
%% get nodes and make matrix
nodes_time = unique([srcnodes;dstnodes]);
numnodes_time = length(nodes_time);

nodes_cur = unique(num3(:,1));
ind = find(num3(:,1)==216);
disp(ind);

index_cur = find(isnan(nodes_cur)==1);
nodes_cur = nodes_cur(1:index_cur(1)-1);
numnodes_cur = length(nodes_cur);

for i=1:length(nodes_cur)
    index = find(nodes_time==nodes_cur(i));
%     disp(index);
    if isempty(index)
        disp(nodes_cur(i))
    end
end

% for i=1:length(nodes_time)
%     if nodes_time(i)==216
%         disp("Aaa");
%     end
% end

