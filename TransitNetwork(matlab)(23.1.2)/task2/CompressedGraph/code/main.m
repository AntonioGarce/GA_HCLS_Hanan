clc,close all,clear 

gen_num = 100;
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
% [num1,txt1,~] = xlsread('KRAKOW - FULL DATA.xlsx','TRAFFIC-DEMAND');
% [num2,txt2,~] = xlsread('KRAKOW - FULL DATA.xlsx','GRAPH');

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
hh = plot(fullGraph,'NodeLabel',unique(fullGraph.Edges.EndNodes));
for i = 1:numnodes
    for j = i+1:numnodes
        if(fulltravelTimeMatrix(i,j) == 0)
            fulltravelTimeMatrix(i,j) = Inf;
            fulltravelTimeMatrix(j,i) = Inf;
        end
    end
end
originTravelTimeMatrix = fulltravelTimeMatrix;


%% Get Compressed nodes
numCompressNodes = 0;
compressednode = {};       % {[1],[2,3,4],[5,6],[7],...}
for nd = 1:numnodes         % index of nodes
    [ind_neigh] = find(~isinf(fulltravelTimeMatrix(nd,:)));
    order_node = getOrderOfNode(nd,fulltravelTimeMatrix);
    if order_node > 2       % if order of node greater than 2
        numCompressNodes = numCompressNodes + 1;
        compressednode{numCompressNodes} = nd;
    elseif order_node == 2
        for ind_n = 1:3
            order_neigh = getOrderOfNode(ind_neigh(ind_n),fulltravelTimeMatrix);
            if order_neigh==1
                numCompressNodes = numCompressNodes + 1;
                compressednode{numCompressNodes} = [ind_neigh(ind_n), ...
                    getCompressedNode(ind_neigh(ind_n), nd, fulltravelTimeMatrix)];
                fulltravelTimeMatrix(compressednode{numCompressNodes},:) = inf;
                break;
            elseif order_neigh>2
                numCompressNodes = numCompressNodes + 1;
                compressednode{numCompressNodes} = getCompressedNode(ind_neigh(ind_n), nd, fulltravelTimeMatrix);
                fulltravelTimeMatrix(compressednode{numCompressNodes},:) = inf;
                break;
            end
        end
    elseif getOrderOfNode(nd,originTravelTimeMatrix) == 1
        nei = ind_neigh(find(ind_neigh~=nd));
        order_neigh = getOrderOfNode(nei,originTravelTimeMatrix);
        if order_neigh~=2 && order_neigh>0
            numCompressNodes = numCompressNodes + 1;
            compressednode{numCompressNodes} = nd;
        end
    end
end

numCompressNodes = length(compressednode);       %number of nodes in a compressed graph

%% get travel and demand matrix for compressed graph
compTravelTimeMatrix = zeros(numCompressNodes); %travel time matrix for compreseed graph
compDemandMatrix = zeros(numCompressNodes);     %demand matrix for compressed graph
compLenMatrix = zeros(numCompressNodes,1);        %length matrix for compressed graph(length of compressed nodes)

for src_ind = 1:numCompressNodes    %src_ind: index of source node(compressed graph)
    srcnode = compressednode{src_ind};
    compLenMatrix(src_ind) = length(srcnode);   %srouce node
    for dst_ind = src_ind+1:numCompressNodes    %dst_ind: index of source node(compressed graph)
        dstnode = compressednode{dst_ind};
        compDemandMatrix(src_ind,dst_ind) = sum(sum(fulldemandMatrix(srcnode,dstnode)));
        compDemandMatrix(dst_ind,src_ind) = compDemandMatrix(src_ind,dst_ind);
        if length(srcnode)~=1 && length(dstnode)~=1     %if both nodes are compressed nodes
            compTravelTimeMatrix(src_ind,dst_ind) = inf;
            compTravelTimeMatrix(dst_ind,src_ind) = inf;
        else
            compTravelTimeMatrix(src_ind,dst_ind) = min([originTravelTimeMatrix(srcnode(1),dstnode(1)), ...
                originTravelTimeMatrix(srcnode(end),dstnode(1)),originTravelTimeMatrix(srcnode(1),dstnode(end)) ...
                originTravelTimeMatrix(srcnode(end),dstnode(end))]);
            compTravelTimeMatrix(dst_ind,src_ind) = compTravelTimeMatrix(src_ind,dst_ind);            
        end 
    end
end

%% get self travel time for compressed n ode
for node_ind = 1:numCompressNodes
    if compLenMatrix(node_ind)~=1
        srcTime = 0;
        srcnode = compressednode{node_ind};
        for k=1:length(srcnode)-1
            srcTime = srcTime + originTravelTimeMatrix(srcnode(k),srcnode(k+1));
        end
        compTravelTimeMatrix(node_ind,node_ind) = srcTime;
    end
end

%% construct the graph and show it
%compWeight: weight matrix for graph (it is same as compTravelMatrix but
%inf cells in compTravelMatrix are changed zero cells.
compWeight = compTravelTimeMatrix;
for k = 1:numCompressNodes
    for kk = 1:numCompressNodes
        if compTravelTimeMatrix(k,kk) == inf
            compWeight(k,kk) = 0;
            compWeight(kk,k) = 0;
        end
    end
end

compressedGraph = graph(compWeight);
figure(2);
plot(compressedGraph);
title("Compressed Graph");
saveas(gcf,'compressed Graph.fig');
f_s = 0;
%% generate initial population

originCompTravelTimeMatrix = compTravelTimeMatrix;
comptravelTimeMatrix = originCompTravelTimeMatrix;

originDemandMatrix = compDemandMatrix;
DemandMatrix = originDemandMatrix;

rou_i = 0;

% endNodes contains nodes which don't have degree of two 
endNodes = find(sum(~isinf(compTravelTimeMatrix))~=3);
% get all shortest paths which has propal length and their lengths for all OD pairs of end nodes.
for k=1:length(endNodes)
    for kk=k+1:length(endNodes)
        shortest_path = shortestpath(compressedGraph,endNodes(k),endNodes(kk));
        lenpath = sum(compLenMatrix(shortest_path));
        if lenpath<=maxroutelen && lenpath>=minroutelen
            rou_i = rou_i + 1;
            shortest_paths{rou_i} = shortest_path;
            shortest_lengths(rou_i) = lenpath;
        end
    end
end

% sort the shortest paths in descend order by lengths
[shortest_lengths,sortind] = sort(shortest_lengths,'descend');
shortest_paths = shortest_paths(sortind);
origin_shortest_paths = shortest_paths;
% numPaths: number of shortest paths
numPaths = length(shortest_paths);
tic
for gen_iter=1:population
    shortest_paths = origin_shortest_paths;
    % select random shortest path as a solution route
    sol = {};
    sol{1} = shortest_paths{randi(numPaths)};
    % generate remain routes
    for iter=2:totalroutes
        %get all nodes which are contained in the selected routes
        totalnodes = [];  
        for s = 1:length(sol)
            totalnodes = [totalnodes sol{s}];
        end
        containedNode = unique(totalnodes);
        disp(length(containedNode));
        %% get the routes
        maxNewNodes = 0;
        maxLen=0;
        maxLen1 = 0;
        selnode1 = 0; % selected node index which have a common point with existing routes
        % for the first part, we select the paths which have maximum new nodes
        if iter<=8
            for k=1:numPaths
                if isempty(shortest_paths{k})
                    continue;
                end
                if length(intersect(containedNode,shortest_paths{k}))==1 && length(shortest_paths{k})>maxLen1
                    selnode1 = k;
                    maxLen1 = length(shortest_paths{k});
                end
                numnewnodes = length(shortest_paths{k})-length(intersect(containedNode,shortest_paths{k})); 
                if maxNewNodes<numnewnodes
                    maxNewNodes = numnewnodes;
                    selnode = k;
                end
            end
            if selnode1 == 0
                sol{iter} = shortest_paths{selnode};
                shortest_paths{selnode} = [];
            else
                sol{iter} = shortest_paths{selnode1};
                shortest_paths{selnode1} = [];
            end
        % second part: find the routes which contain the remain nodes
        elseif iter>8
            newPath = [];
            for k=1:numCompressNodes
                lenpath =0;
                if isempty(find(containedNode==k,1))
                    while(lenpath<minroutelen || lenpath>maxroutelen)
                        newPathL = shortestpath(compressedGraph,endNodes(randi(length(endNodes))),k);
                        newPathR = shortestpath(compressedGraph,k,endNodes(randi(length(endNodes))));
                        newPath = [newPathL,newPathR(2:end)];
                        lenpath = sum(compLenMatrix(newPath));
                    end
                    sol{iter} = newPath;
                    break;
                end
            end
            if isempty(newPath)
                while(length(newPath)<minroutelen || length(newPath)>maxroutelen)
                    randpaths = kShortestPath(comptravelTimeMatrix,endNodes(randi(length(endNodes))),endNodes(randi(length(endNodes))),2);
                    if length(randpaths)<2
                        continue;
                    end
                    newPath = randpaths{2};
                    lenpath = sum(compLenMatrix(newPath));
                end
                sol{iter} = newPath;
            end
        end
    end
    generation{gen_iter} = sol;
end
toc
%% test
fullSolution = {};
for k=1:length(sol)
    full_route = [];
    comp_route = sol{k};
    for kk =1:length(comp_route)
        full_route = [full_route, compressednode{comp_route(kk)}];
    end
    fullSolution{k} = full_route;
end
ffval_full  = ff_fun(fullSolution,fulldemandMatrix,fulltravelTimeMatrix,0.5,0.3,0.2,0.5,0.3,0.2,5,10);

for k=1:length(fullSolution)
    writematrix(fullSolution{k},"result.xlsx",'Sheet','solution','WriteMode','append');
end

%% 
% for k=1:population
%     generation{k} = sol;
% end

objMatrix = zeros(population,1);
FFMatrix = {};
for pop_i = 1:population
    solution = generation{pop_i};
    ffval = ff_fun(solution,compDemandMatrix,originCompTravelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
   
    fullSolution = {};
    for k=1:length(solution)
        full_route = [];
        comp_route = solution{k};
        for kk =1:length(comp_route)
            full_route = [full_route, compressednode{comp_route(kk)}];
        end
        fullSolution{k} = full_route;
    end
    ffval_full  = ff_fun(fullSolution,fulldemandMatrix,fulltravelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
    objMatrixFull(pop_i) = ffval_full.TF;
    
    objMatrix(pop_i) = ffval.TF;
    FFMatrix{pop_i} = ffval;
end
% figure(2);
tic
for genIter = 1:gen_num
    [bestFitVal, bestFitInd] = min(objMatrix); 
    bestSolution = generation{bestFitInd};
    [bestFitValFull, bestFitIndFull] = min(objMatrixFull);
    
    meanFitValFull = sum(objMatrixFull)/population;
    meanFitVal = sum(objMatrix)/population;
    disp("generation: " + num2str(genIter));
    disp("ATT: " + num2str(FFMatrix{bestFitInd}.ATT));
    disp("PT0: " + num2str(FFMatrix{bestFitInd}.PT0));
    disp("PT1: " + num2str(FFMatrix{bestFitInd}.PT1));
    disp("PT2: " + num2str(FFMatrix{bestFitInd}.PT2));
    hold on;

    figure(3)

    title("best:"+num2str(bestFitVal)+" mean:"+num2str(meanFitVal));
    plot(genIter,bestFitVal,'bo-',genIter,meanFitVal,'rd--');
    axis([1 gen_num -20 20]);
    hold on
    
    saveas(gcf,'aa.fig')
    minFit = min(objMatrix);
    maxFit = max(objMatrix);
    
    for pop_i = 1:population
        fitMatrix(pop_i) = 1/(1+3*(objMatrix(pop_i)-minFit)/(maxFit-minFit+0.0001));
    end

    fitDistribution = zeros(population+1,1);
    for pop_i = 1:population
        fitDistribution(pop_i+1) = fitDistribution(pop_i)+fitMatrix(pop_i);
    end
    fitSum = fitDistribution(end);
    fitDistribution = fitDistribution/fitSum;
    newGeneration = cell(30,1);
    containedInd = [];
    numpop = 0;
    for pop_i = 1:population
        if(rand<1-crossover_prob)  %eleitism percentage 20%, 
            if(isempty(containedInd))
                parent = bestFitInd;
            else
                parent = select_individual(fitDistribution);
                while(~isempty(find(containedInd==parent, 1)))
                    parent = select_individual(fitDistribution);
                end
            end
            containedInd = [containedInd, parent];
            newGeneration{pop_i} = generation{parent};
        else          %crossover percentage 80%
            parent = select_individual(fitDistribution);
            newParent=crossover(generation{parent},compLenMatrix,maxroutelen,minroutelen,10);
            newGeneration{pop_i} = newParent;
        end
        if rand<mutation_prob   %Mutation with percentage of 30%.
            [mutatedIndividual, ismutated] = mutation(originCompTravelTimeMatrix,newGeneration{pop_i},compLenMatrix,maxroutelen,minroutelen,10);
            if ismutated
                newGeneration{pop_i} = mutatedIndividual;
            end
        end
    end

    objMatrix = zeros(population,1);
    FFMatrix = {};
    for pop_i = 1:population
        solution = newGeneration{pop_i};
        ffval = ff_fun(solution,compDemandMatrix,originCompTravelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
        objMatrix(pop_i) = ffval.TF;
        FFMatrix{pop_i} = ffval;
    end
%     %% HCLS
    MemePop = {};
    MemePop_ind = zeros(Meme_population,1);
    for pop_i =1:Meme_population
        exist_f = 0;
        while ~exist_f
            meme_ind = randi(population);
            indd = find(MemePop_ind==meme_ind);
            exist_f = isempty(indd);
        end  
        MemePop_ind(pop_i) = meme_ind;
        MemePop{pop_i} = newGeneration{MemePop_ind(pop_i)};
    end

    [MemePopSol,fitTF]=HCLS(MemePop,objMatrix(MemePop_ind),originCompTravelTimeMatrix,originDemandMatrix,HCLS_MaxIter);
    for pop_i = 1:Meme_population
        newGeneration{MemePop_ind(pop_i)} = MemePopSol{pop_i};
        objMatrix(MemePop_ind(pop_i)) = fitTF(pop_i);
    end
    generation = newGeneration;
end
toc
%% get best full solution
fullSolution = {};
for k=1:length(bestSolution)
    full_route = [];
    comp_route = bestSolution{k};
    for kk =1:length(comp_route)
        full_route = [full_route, compressednode{comp_route(kk)}];
    end
    fullSolution{k} = full_route;
end
ffval_full  = ff_fun(fullSolution,fulldemandMatrix,fulltravelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);

for k=1:length(compressednode)
    writematrix(compressednode{k},"result.xlsx",'Sheet','compressed nodes','WriteMode','append');  
end
for k=1:length(fullSolution)
    writematrix(fullSolution{k},"result.xlsx",'Sheet','solution','WriteMode','append');
end

%% other part
Name = {'ATT';'d0';'d1';'d2';'dud'};
Value = [ffval_full.ATT;ffval_full.PT0; ...
    ffval_full.PT1;ffval_full.PT2;ffval_full.PUD];
T = table(Name,Value);
writetable(T,"result.xlsx",'Sheet', "Result Finess value");

disp("solution:");
disp(fullSolution);
disp("ATT:"+num2str(ffval_full.ATT));
disp("d0:"+num2str(ffval_full.PT0));
disp("d1:"+num2str(ffval_full.PT1));
disp("d2:"+num2str(ffval_full.PT2));

disp("end!!!");



