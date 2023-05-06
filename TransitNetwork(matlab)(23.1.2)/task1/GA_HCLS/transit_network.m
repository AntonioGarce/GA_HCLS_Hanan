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
    nodenames{indsrc} = txt2(i+1,4);
    nodenames{inddst} = txt2(i+1,5);
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

currSolution={[125,	124,	123,	122,	121,	120,	119,	118,	117,...
    103,	116,	115,	885,	112,	452,	814,	813,	812,...
    811,	810,	809,	808,	807,	806,	821,	822,	823,...
    300,	301], [301,	300,	823,	822,	827,	833,	801,	851,...
    852,	853,	854], [702,	657,	656,	655,	649,	648,	647,...
    646,	645,	644,	643,	642,	640,	816,	815,	807,...
    806,	801,	760,	839,	381,	382,	383,	416],[125,	124,...
    123,	122,	194,	193,	107,	106,	105,	104,	103,...
    102,	101,	100,	805,	804,	803,	802,	851,	801,...
    833,	827,	828,	829,	354,	355,	356,	169,	357,...
    358,	359,	360],[125,	124,	123,	122,	121,	120,	119,...
    118,	117,	103,	102,	101,	100,	805,	804,	803,	802,...
    901,	839,	381,	382,	383,	416],[584,	583,	627,	634,...
    635,	646,	645,	644,	643,	642,	641,	556,	820,...
    819,	818,	817,	821,	822,	824,	825,	345,	453,	346,	347],...
    [562,	743,	561,	553,	559,	558,	557,	556,	820,	819,	...
    818,	817,	821,	822,	827,	828,	829,	354,	355,	356,	169,...
    357,	358,	359,	360],[139,	140,	141,	143,	161,	162,	135,	...
    118,	117,	103,	102,	101,	100,	805,	804,	803,	802,	...
    809,	619,	294,	670,	671,	674,	501,	645,	646,	647,	648,...
    649,	655,	656,	657,	702],[218,	217,	216,	215,	214,	213,	889,...
    212,	107,	106,	105,	104,	103,	102,	101,	100,	805,	804,	...
    803,	802,	851,	801,	806,	817,	818,	819,	820,	556,	557,...
    558,	559,	553],[711,	678,	674,	501,	644,	643,	642,	641,	556,...
    557,	558,	559,	553,	552,	551,	528,	526,	527,	706,	291,	...
    292,	293],[702,	657,	656,	655,	649,	648,	647,	646,	645,	644,...
    643,	642,	641,	556,	820,	819,	818,	817,	821,	822,	827,	...
    828,	829,	354,	355,	356,	169,	357],[139,	140,	141,	143,	161,...
    162,	135,	118,	117,	103,	116,	115,	885,	112,	452,	814,	...
    813,	812,	811,	810,	809,	802,	851,	801,	833,	827,	828,	...
    829,	354,	355,	356,	169,	357],[293,	292,	291,	706,	527,	526,...
    528,	551,	552,	553,	559,	558,	557,	556,	640,	816,	815,	...
    807,	808,	809,	802,	901,	839,	381],[293,	292,	291,	706,	527,...
    526,	528,	525,	524,	523,	522,	501,	872,	819,	818,	817,	...
    821,	822,	827,	833,	834,	835,	381,	382,	383,	416],[562,	743,...
    561,	553,	559,	558,	557,	556,	640,	816,	815,	807,	808,	...
    809,	802,	901,	839,	381],[711,	678,	674,	671,	670,	294,	619,...
    809,	802,	851,	801,	833,	827,	824,	825,	345,	453,	346,	...
    347],[218,	217,	216,	215,	214,	213,	889,	212,	107,	106,	105,...
    136,	118,	135,	162,	161,	144,	147,	148],[562,	743,	561,	553,...
    552,	551,	525,	524,	523,	522,	501,	872,	819,	807,	808,	...
    809,	810,	811,	812,	813,	814,	452,	112,	885,	115,	116,	...
    103,	104,	105,	106,	107,	212],[584,	583,	627,	634,	635,	646,...
    645,	644,	643,	642,	640,	816,	815,	807,	806,	801,	833,	...
    827,	828,	829,	354,	355,	356,	169,	357,	358,	359,	360],   ...
    [702,	657,	656,	655,	649,	648,	647,	646,	645,	501,	674,	...
    671,	670,	294,	619,	809,	802,	803,	804,	805],[416,	383,	382,...
    381,	839,	901,	802,	809,	619,	294,	670,	671,	674,	501,	...
    645,	646,	635,	634,	627,	583,	584],[293,	292,	291,	706,	527,...
    526,	528,	525,	524,	523,	522,	501,	872,	819,	807,	806,	...
    801,	851,	802,	803,	804,	805,	100,	101,	102,	103,	117,	...
    118,	135,	162,	161,	144,	147,	148]};

currSol_ind = {};

for rou_i = 1:length(currSolution)
    route = currSolution{rou_i};
    roulen = length(route);
    rou_ind = zeros(1,roulen);
    for node_i = 1:roulen
        rou_ii = find(nodes==route(node_i));
        rou_ind(node_i) = rou_ii(1);
    end
    currSol_ind{rou_i} = rou_ind;
end

f_s = 0;

%% check routes full
numsol = 0;
% sol = zeros;
sol={};
sol_ind={};
while(numsol<population-1)
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

sol{population} = currSolution;
sol_ind{population} = currSol_ind;

generation = sol_ind;
objMatrix = zeros(population,1);
FFMatrix = {};
for pop_i = 1:population
    solution = generation{pop_i};
    ffval = ff_fun(solution,demandMatrix,originTravelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
    objMatrix(pop_i) = ffval.TF;
    FFMatrix{pop_i} = ffval;
end

for genIter = 1:gen_num
%     fitMatrix = zeros(population,1);
%     FFMatrix = {};
%     for pop_i = 1:population
%         solution = generation{pop_i};
%         ffval = ff_fun(solution,demandMatrix,originTravelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
%         fitMatrix(pop_i) = ffval.TF;
%         FFMatrix{pop_i} = ffval;
%     end
    [bestFitVal, bestFitInd] = min(objMatrix); 
    meanFitVal = sum(objMatrix)/population;
    disp("generation: " + num2str(genIter));
    disp("ATT: " + num2str(FFMatrix{bestFitInd}.ATT));
    disp("PT0: " + num2str(FFMatrix{bestFitInd}.PT0));
    disp("PT1: " + num2str(FFMatrix{bestFitInd}.PT1));
    disp("PT2: " + num2str(FFMatrix{bestFitInd}.PT2));
    hold on;
    title("best:"+num2str(bestFitVal)+" mean:"+num2str(meanFitVal));
    plot(genIter,bestFitVal,'bo-',genIter,meanFitVal,'rd--');
    pause(0.2)
    axis([1 gen_num -30 0]);
    minFit = min(objMatrix);
    maxFit = max(objMatrix);
    
    for pop_i = 1:population
        fitMatrix(pop_i) = 1/(1+3*(objMatrix(pop_i)-minFit)/(maxFit-minFit));
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
            newParent=crossover(generation{parent},maxroutelen,minroutelen,10);
            newGeneration{pop_i} = newParent;
        end
        if rand<mutation_prob   %Mutation with percentage of 30%.
            [mutatedIndividual, ismutated] = mutation(originTravelTimeMatrix,newGeneration{pop_i},maxroutelen,minroutelen,10);
            if ismutated
                newGeneration{pop_i} = mutatedIndividual;
            end
        end
    end

    objMatrix = zeros(population,1);
    FFMatrix = {};
    for pop_i = 1:population
        solution = newGeneration{pop_i};
        ffval = ff_fun(solution,demandMatrix,originTravelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
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

    [MemePopSol,fitTF]=HCLS(MemePop,objMatrix(MemePop_ind),originTravelTimeMatrix,demandMatrix,HCLS_MaxIter);
    for pop_i = 1:Meme_population
        newGeneration{MemePop_ind(pop_i)} = MemePopSol{pop_i};
        objMatrix(MemePop_ind(pop_i)) = fitTF(pop_i);
    end
    generation= newGeneration;
end

bestsolution_route = {};
%%
xlswrite("bestRoute.xlsx",{'Line and its tram stops'},'A1:A1')
xlswrite("bestRoute.xlsx",{'TramStopCode'},'B1:B1')
xlswrite("bestRoute.xlsx",{'TravelTime'},'C1:C1')
for i=1:length(solution)
    bestsolution_route{i}{1} = nodenames(generation{bestFitInd}{i});
    bestsolution_route{i}{2} = nodes(generation{bestFitInd}{i});
    roulen = length(generation{bestFitInd}{i});
    tim = zeros(roulen-1,1);
    for j = 1:roulen-1
        tim(j) = originTravelTimeMatrix(generation{bestFitInd}{i}(j),generation{bestFitInd}{i}(j+1));
    end
    bestsolution_route{i}{3} = tim;
end

f_ind=3;

for cell_ind =1:length(bestsolution_route)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T=bestsolution_route{cell_ind};
    t=T{1};
    for ii=1:length(t)
        t1(ii)=cellstr(t{ii});
    end
    
    tn=t1';

    t2=T{2};     
    tt = T{3};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L={['Line' num2str(cell_ind)]};
    A=['A' num2str(f_ind) ':' 'A' num2str(f_ind)];
    len=length(t2);
    xlswrite("bestRoute.xlsx",L,A);
    
    xlswrite("bestRoute.xlsx",tn,['A' num2str(f_ind+1) ':A' num2str(f_ind+len)])
    xlswrite("bestRoute.xlsx",t2,['B' num2str(f_ind+1) ':B' num2str(f_ind+len)])
    xlswrite("bestRoute.xlsx",tt,['C' num2str(f_ind+2) ':C' num2str(f_ind+len)])
    f_ind=len+f_ind+2;
end

%% other part
Name = {'ATT';'d0';'d1';'d2';'dud'};
Value = [FFMatrix{bestFitInd}.ATT;FFMatrix{bestFitInd}.PT0; ...
    FFMatrix{bestFitInd}.PT1;FFMatrix{bestFitInd}.PT2;FFMatrix{bestFitInd}.PUD];
T = table(Name,Value);
writetable(T,"bestRoute.xlsx",'Sheet', "Result Finess value");
disp("end!!!");



