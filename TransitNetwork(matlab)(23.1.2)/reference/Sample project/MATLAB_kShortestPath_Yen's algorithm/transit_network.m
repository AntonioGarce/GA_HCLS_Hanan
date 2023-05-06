clc,close all,clear 

gen_num = 2;
population=30;
totalroutes = 22;
minroutelen = 10;
maxroutelen = 35;
%% read excel file
[num1,txt1,~] = xlsread('KRAKOW - FULL DATA.xlsx','TRAFFIC-DEMAND');
[num2,txt2,~] = xlsread('KRAKOW - FULL DATA.xlsx','GRAPH');

srcnodes = num2(:,1);
dstnodes = num2(:,2);
demands = num1(:,3);
%% get nodes and make matrix
nodes = unique([srcnodes;dstnodes]);
nodenames = cell(1,length(nodes));
numnodes = length(nodes);
demandMatrix = zeros(numnodes,numnodes);

%% set elements of demand matrix.
srcnodes = num1(:,1);
dstnodes = num1(:,2);
for i=1:length(srcnodes)
    indsrc = find(nodes==srcnodes(i));
    inddst = find(nodes==dstnodes(i));
    if(~isempty(indsrc) && ~isempty(inddst))
        demandMatrix(indsrc,inddst) = demandMatrix(indsrc,inddst)+demands(i);
        demandMatrix(inddst,indsrc) = demandMatrix(indsrc,inddst);
    end
end
%% set elements of timetravel matrix
travelTimeMatrix = zeros(numnodes,numnodes);

srcnodes = num2(:,1);
dstnodes = num2(:,2);
time = num2(:,3);

for i=1:length(srcnodes)
    indsrc = find(nodes==srcnodes(i));
    inddst = find(nodes==dstnodes(i));
    travelTimeMatrix(indsrc,inddst) = time(i);
    travelTimeMatrix(inddst,indsrc) = time(i);
    nodenames{indsrc} = txt2(i+1,3);
    nodenames{inddst} = txt2(i+1,4);
end

for i = 1:numnodes
    for j = i+1:numnodes
        if(travelTimeMatrix(i,j) == 0)
            travelTimeMatrix(i,j) = Inf;
            travelTimeMatrix(j,i) = Inf;
        end
    end
end
originTravelTimeMatrix = travelTimeMatrix;

f_s = 0;

%% check routes full
numsol = 0;
% sol = zeros;
sol={};
sol_ind={};
while(numsol<population)
    travelTimeMatrix = originTravelTimeMatrix;
    numroute = 0;
    solution = {};
    iter = 0;

    while ~isFullRoutes(travelTimeMatrix )
        iter = iter +1;
        if(iter == 10000 || numroute>=totalroutes)  %if the solution is in infinitive loop, start next solution
            iter=10000;
            break;
        end
        %%find isolated nodes
        isonodes = find_isonodes(travelTimeMatrix);
        %random select of isonode
        if(~isempty(isonodes))
            stnode = isonodes(randi([1,length(isonodes)]));
        else
            for i = 1:numnodes
                if(sum(isinf(travelTimeMatrix),2) ~= (numnodes-1))
                    stnode = i;
                    break;
                end
            end
        end
        
        %% find the routes started from isolated nodes
       path =  getFullRoute(travelTimeMatrix, stnode);  %get full route 
        %  if(length(path))
       lenpath = length(path)-1;
       lowflag = 1;
        if lenpath<minroutelen              %if the length of route is not full
            % extend the route using the nodes already used
            while(length(path)<=minroutelen && lowflag==1)      
                selnode = path(end);
                linkv = originTravelTimeMatrix(selnode,:);
                linkv(selnode) = Inf;
                for i=1:length(path)
                    linkv(path(i))=Inf;
                end
                ind = find(~isinf(linkv));
                if(isempty(ind))
                    lowflag = 0;
                    break;
                end
                prevnode = selnode;
                selnode = ind(randi([1,length(ind)]));
                path = [path selnode]
            end
        end
        lenpath = length(path)-1;
        if(lowflag == 0)        %if the length of extended route is smaller than minroutelen 
            continue;
        elseif (lenpath>=minroutelen && lenpath<maxroutelen) %if the route satisfies the boundary of length
            numroute = numroute+1;
            solution{numroute} = path;
            for i = 1:length(path)
               travelTimeMatrix(path(i),path) = Inf;
               travelTimeMatrix(path(i),path(i)) = 0;
            end
        else                                                %if the length is longer than maxroutelen
            numsubroutes = floor(lenpath/maxroutelen);  
            remain = mod(lenpath,maxroutelen);
            
            for i = 1:numsubroutes-1
                route = path((i-1)*maxroutelen+1:i*maxroutelen+1);
                numroute = numroute+1;
                solution{numroute} = route;
            end
            
            if(remain==0)                  %if the length of route is integer times of maxroutelen
                numroute = numroute+1;
                solution{numroute} = path(lenpath-maxroutelen+1:end);
            else                           
                lastsubroutelen = floor((remain+maxroutelen)/2);
                numroute = numroute+1;
                solution{numroute} = path(lenpath-maxroutelen-remain+1:lenpath-lastsubroutelen+1);
                numroute = numroute+1;
                solution{numroute} = path(lenpath-lastsubroutelen+1:end);
            end
            %update the travelTimeMatrix(remove the edge that be contained
            %in routes.
            for i = 1:length(path)
               travelTimeMatrix(path(i),path) = Inf;
               travelTimeMatrix(path(i),path(i)) = 0;
            end
        end
        disp("path"+num2str(numroute));
        disp(path);
    end
    if(iter<10000)        %if the solution is gotten correctly
        if length(solution)<totalroutes
            numextrasol = totalroutes - length(solution);
            while numroute<totalroutes
                path = getRouteWithLength(originTravelTimeMatrix,randi([1,numnodes]),randi([minroutelen,maxroutelen]));
                if ~isempty(path)
                    numroute = numroute+1;
                    solution{numroute}=path;
                    disp("path"+num2str(numroute));
                    disp(path);
                end
            end
        end
        numsol = numsol + 1;
        solution_route = {};
        
        for i=1:length(solution)
            solution_route{i} = nodes(solution{i});
        end
        sol{numsol} = solution_route;
        sol_ind{numsol} = solution;
    end
end

generation = sol_ind;
for genIter = 1:gen_num
    fitMatrix = zeros(population,1);
    FFMatrix = {};
    for pop_i = 1:population
        solution = generation{pop_i};
        ffval = ff_fun(solution,demandMatrix,originTravelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
        fitMatrix(pop_i) = ffval.TF;
        FFMatrix{pop_i} = ffval;
    end
    [bestFitVal bestFitInd] = min(fitMatrix); 
    meanFitVal = sum(fitMatrix)/population;
    disp("generation: " + num2str(genIter));
    disp("ATT: " + num2str(FFMatrix{bestFitInd}.ATT));
    disp("PT0: " + num2str(FFMatrix{bestFitInd}.PT0));
    disp("PT1: " + num2str(FFMatrix{bestFitInd}.PT1));
    disp("PT2: " + num2str(FFMatrix{bestFitInd}.PT2));
    hold on;
    title("best:"+num2str(bestFitVal)+" mean:"+num2str(meanFitVal));
    plot(genIter,bestFitVal,'bo-',genIter,meanFitVal,'rd--');
    pause(0.2)
    axis([1 gen_num -20 0]);
    minFit = min(fitMatrix);
    maxFit = max(fitMatrix);
    
    for pop_i = 1:population
%         fitMatrix(pop_i) = 4*(fitMatrix(pop_i)-minFit)^2/(maxFit-minFit)^2+1;
        fitMatrix(pop_i) = 1/(1-minFit+fitMatrix(pop_i))^2;
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
        if(rand<0.2)  %eleitism percentage 20%, 
            parent = select_individual(fitDistribution);
            while(~isempty(find(containedInd==parent, 1)))
                parent = select_individual(fitDistribution);
            end
            newGeneration{pop_i} = generation{parent};
        else          %crossover percentage 80%
            parent = select_individual(fitDistribution);
            newParent=crossover(generation{parent},maxroutelen,minroutelen,10);
            newGeneration{pop_i} = newParent;
        end
        if rand<0.3
            [mutatedIndividual, ismutated] = mutation(originTravelTimeMatrix,newGeneration{pop_i},maxroutelen,minroutelen,10);
            if ismutated
                newGeneration{pop_i} = mutatedIndividual;
            end
        end
    end
    generation= newGeneration;
end

bestsolution_route = {};

for i=1:length(solution)
    bestsolution_route{i} = [nodenames(generation{bestFitInd}{i}),...
        nodes(generation{bestFitInd}{i});];
end
f_ind=2;
for cell_ind = 1:3%length(bestsolution_route)
    cell_len=length(bestsolution_route{cell_ind});
    l_ind = f_ind + cell_len;
    strLineName = cellstr(['Line',num2str(cell_ind)]);
    range=sprintf('A%i:B%i',f_ind-1,f_ind-1);
    xlswrite("bestRoute.xlsx",strLineName, 'Range',range);
    range=sprintf('A%i:A%i',f_ind,l_ind);
    xlswrite("bestRoute.xlsx",bestsolution_route{cell_ind}, 'Range',range);
    xlswrite("bestRoute.xlsx",' ', 'Range',sprintf('A%i:A%i',l_ind,l_ind));
    f_ind=f_ind+cell_len+2;
end
Name = {'ATT';'d0';'d1';'d2';'dud'};
Value = [FFMatrix{bestFitInd}.ATT;FFMatrix{bestFitInd}.PT0; ...
    FFMatrix{bestFitInd}.PT1;FFMatrix{bestFitInd}.PT2;FFMatrix{bestFitInd}.PUD];
T = table(Name,Value);
writetable(T,"bestRoute.xlsx",'Sheet', "Result Finess value");
disp("end!!!");



