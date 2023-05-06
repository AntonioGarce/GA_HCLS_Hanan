numnodes = 149;

degree_p = [9.4 68.5 14.8 7.3];
degree = [1 2 3 4];
time_p = [47.3 41.7 8 2.4 0.6];
time = [1 2 3 4 6];

% for i = 1:numnodes
%     for j= 1:numnodes
%         if (i~=j) && (cracow(i,j)==0)
%             cracow(i,j)=inf;
%         end
%     end
% end

% weight=lt(cracow,1000);    %make weight matrix to detrmine the links between nodes

%% Graph Construction 
correct_graph = 0;
while(~correct_graph)
    cracow = random_graph(numnodes,degree_p,degree,time_p,time);
    G=graph(cracow,'omitselfloops');
    for k = 1:numnodes
       if(isempty(G.shortestpath(1,k,"Method","unweighted")))
          break;  
       end
       if k == numnodes
          correct_graph = 1;
       end
    end
end
hh= plot((G),'layout', 'force');
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
saveas(gcf,'cracow.fig');

for k=1:height(G.Edges)
    
end
demandData = readmatrix('KRAKOW - FULL DATA.xlsx','Sheet','TRAFFIC-DEMAND');

srcnodes = demandData(:,1);
dstnodes = demandData(:,2);
demands = demandData(:,3);

numDemands = length(demands);
numZeroDemands = sum(demands==0);
pr0 = numZeroDemands/numDemands;
numNonZeroDemands = numDemands - numZeroDemands;
lambda_e = numNonZeroDemands/(sum(demands)-numNonZeroDemands);
mu = 1/lambda_e;
newDemandData = zeros(numDemands,3);
demandMatrix = zeros(numnodes);
for k = 1:numDemands
    srcnode = randi(numnodes);
    dstnode = randi(numnodes);
    while(demandMatrix(srcnode,dstnode)~=0 ||srcnode == dstnode)
        srcnode = randi(numnodes);
        dstnode = randi(numnodes);
    end
    newDemandData(k,1) = srcnode;
    newDemandData(k,2) = dstnode;
    
    if (rand<=0.356) 
        newDemandData(k,3) = 0;
    else
        rndemand = (ceil(exprnd(mu)));
        demandMatrix(srcnode,dstnode) = rndemand;
        newDemandData(k,3) = rndemand;
    end
end
ind = find(cracow~=0);

newTravelTime = zeros(length(find(cracow~=0))/2,3);
indk = 0;
for k=1:149
    for kk=k+1:149
        if cracow(k,kk)~=0
            indk = indk+1;
            newTravelTime(indk,:) = [k, kk, cracow(k,kk)];
        end
    end
end

timeHeader = ["FromTramStop","ToTramStop","Time"];
writematrix(timeHeader, "Cracow.xlsx", 'Sheet',"Connections",'Range',"A1:C3");
writematrix(newTravelTime, "Cracow.xlsx", 'Sheet',"Connections",'Range',"A2:C"+num2str(length(newTravelTime)+1));

demandHeader = ["FromTramStop","ToTramStop","Demand"];
writematrix(demandHeader, "Cracow.xlsx", 'Sheet',"Traffic-Demand",'Range',"A1:C3");
writematrix(newDemandData, "Cracow.xlsx", 'Sheet',"Traffic-Demand",'Range',"A2:C"+num2str(numDemands+1));

save 'cracow.mat'