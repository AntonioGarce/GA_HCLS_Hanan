clc,clear,close all 

gen_num = 150;
population=30;
crossover_prob = 0.5;
mutation_prob = 0.3;
Meme_population = 4;
HCLS_MaxIter = 22;
totalroutes = 22;
minroutelen = 10;
maxroutelen = 35;
%% read excel file
[num1,txt1,~] = xlsread('Cracow.xlsx','TRAFFIC-DEMAND');
[num2,txt2,~] = xlsread('Cracow.xlsx','Connections');

srcnodes = num2(:,1);
dstnodes = num2(:,2);
demands = num1(:,3);
%% get nodes and make matrix
nodes = unique([srcnodes;dstnodes]);
% nodenames = cell(1,length(nodes));
numnodes = length(nodes);
fulldemandMatrix = zeros(numnodes,numnodes);

%% set elements of demand matrix.
srcnodes = num1(:,1);
dstnodes = num1(:,2);
for i=1:length(srcnodes)
    indsrc = find(nodes==srcnodes(i));
    inddst = find(nodes==dstnodes(i));
    if(~isempty(indsrc) && ~isempty(inddst))
        fulldemandMatrix(indsrc,inddst) = fulldemandMatrix(indsrc,inddst)+demands(i);
        fulldemandMatrix(inddst,indsrc) = fulldemandMatrix(indsrc,inddst);
    end
end
%% set elements of timetravel matrix
fulltravelTimeMatrix = zeros(numnodes,numnodes);

srcnodes = num2(:,1);
dstnodes = num2(:,2);
time = num2(:,3);

for i=1:length(srcnodes)
    indsrc = find(nodes==srcnodes(i));
    inddst = find(nodes==dstnodes(i));
    fulltravelTimeMatrix(indsrc,inddst) = time(i);
    fulltravelTimeMatrix(inddst,indsrc) = time(i);
end
fullGraph = graph(fulltravelTimeMatrix);
figure(1);
hh = plot(fullGraph,'NodeLabel',unique(fullGraph.Edges.EndNodes),'layout', 'force');
numnodes = height(fullGraph.Nodes);
nodeColor = zeros(numnodes,3);
nodeMarker = zeros(1,numnodes);
for k = 1:numnodes
    switch(fullGraph.degree(k))
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

old_solution = readSolution('bestRoute1.xlsx','Sheet1');
new_solution = readSolution('KRAKOW - FULL DATA.xlsx', 'CurrentSolution');

num_solution = length(old_solution);
similarRoute_i = zeros(num_solution,1);
num_CommNode = zeros(num_solution,1);
for k=1:num_solution
    maxComm = 0;
    for j = k+1:num_solution
        if ~any(similarRoute_i==j)
            comm = length(intersect(old_solution{k},new_solution{j}));
            if maxComm<=comm
                maxComm = comm;
                similarRoute_i(k,1) = j;
                num_CommNode(k,1) = comm;
            end
        end
    end
end

[num_CommNode, sort_ind] = sort(num_CommNode,'descend');
similarRoute_i = similarRoute_i(sort_ind);

% for k=1:num_solution
%     T={['Transition' num2str(k)]};
%     xlsxwrite('transition.xlsx', )
% end
% new_solution = readSolution()
disp('end');
