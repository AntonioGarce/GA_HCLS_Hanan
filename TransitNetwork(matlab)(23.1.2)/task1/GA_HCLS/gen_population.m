  travelTimeMatrix = importdata("Mumford2TravelTimes.txt");
demandMatrix = importdata("Mumford2Demand.txt");
%%find the isolated nodes.
[numnodes,~]=size(travelTimeMatrix);

ind = find(sum(isinf(travelTimeMatrix),2)>=numnodes-2);

for i=2:110
    [shortestPaths, totalCosts] = kShortestPath(travelTimeMatrix, 1, i, 100);
    if(length(shortestPaths)>=1)
        fprintf('i=%i, num=%i\n',i,length(shortestPaths));        
    end
end