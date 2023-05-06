function [FFVAL] = ff_fun(solution,demandMatrix,travelTimeMatrix,w1,w2,w3,omega1,omega2,omega3,P1,P2)
    tic
    [numnodes, ~] = size(travelTimeMatrix);      %get total number of nodes
    numroutespersol = length(solution);         %number of routes per solution
    optTravelTime = Inf*ones(numnodes,numnodes); %matrix for optical travel time
    numTranOfBest = Inf*ones(numnodes,numnodes);      %matrix which shows transfer of best method
    bestRoute = {};        %best route cell
    
    for rou_i=1:numroutespersol                 %for each route
        route = solution{rou_i};                %get route
        roulen = length(route);                %get length of route
        for nodef_i = 1:roulen                 %index of first node
            nodef = route(nodef_i);             %get first node of this route
            travelTime = 0;    
            prevnode = nodef;                   %set previouse node as nodes_1
            for nodel_i = nodef_i+1:roulen        %index of last node
                nodel = route(nodel_i);         %get last node
                travelTime = travelTime + travelTimeMatrix(prevnode,nodel);   %update travel time 
                if(travelTime<optTravelTime(nodef,nodel))    
                    optTravelTime(nodef,nodel) = travelTime;                       %set travel time to the travelTime matrix
                    optTravelTime(nodel,nodef) = travelTime;                       %set travel time to the travelTime matrix
                end
                prevnode = nodel;                                             %update prevnode to reduce calcuation effort
%                 canBeT0(nodef,nodel) = 1;                                     %tick for direct trip
%                 canBeT0(nodel,nodef) = 1; 
                numTranOfBest(nodef,nodel) = 0;                               %update number of transfers of best method
                numTranOfBest(nodel,nodef) = 0;
                
                bestRoute{nodef}{nodel}.route = [rou_i];                            % update best route  
                bestRoute{nodef}{nodel}.interStop = [];
                bestRoute{nodel}{nodef} = bestRoute{nodef}{nodel};
            end
        end
    end
    % get optimal solution using bellman's opimility principle
    maxTransfer = 5;
    numTransfer = 0;
    while(sum(sum(isinf(optTravelTime)))~=numnodes && numTransfer<=maxTransfer)
        numTransfer = numTransfer + 1;
        for nodef = 1:numnodes                      %first node
            for nodel = nodef+1:numnodes            %last node
                for nodestop = 1:numnodes           %intermediate node which is intersect of between routes
                    if(nodestop==nodef || nodestop==nodel)
                        travelTime = inf;
                    else
                        travelTime = optTravelTime(nodef,nodestop)+optTravelTime(nodestop,nodel);
                    end
                    if(travelTime<optTravelTime(nodef,nodel))
                        optTravelTime(nodef,nodel) = travelTime;
                        optTravelTime(nodel,nodef) = travelTime;
                        numTranOfBest(nodef,nodel) = numTranOfBest(nodef,nodestop)+numTranOfBest(nodel,nodestop)+1;
                        numTranOfBest(nodel,nodef) = numTranOfBest(nodef,nodel);
                        bestRoute{nodef}{nodel}.route = [bestRoute{nodef}{nodestop}.route, bestRoute{nodel}{nodestop}.route];
                        bestRoute{nodef}{nodel}.interStop = [bestRoute{nodef}{nodestop}.interStop, nodestop, bestRoute{nodel}{nodestop}.interStop];
                        bestRoute{nodel}{nodef} = bestRoute{nodef}{nodel};
                    end
                end
            end
        end
    end
    
    for node = 1:numnodes
        optTravelTime(node,node) = 0;
    end
    FFVAL.TD = 0;       %total passenger demand
    FFVAL.TIVT = 0;
    FFVAL.T0 = 0;
    FFVAL.T1 = 0;
    FFVAL.T2 = 0;
    for nodef=1:numnodes
        for nodel = nodef:numnodes
            if(~isinf(optTravelTime(nodef,nodel)))
                FFVAL.TIVT = FFVAL.TIVT + demandMatrix(nodef,nodel)*optTravelTime(nodef,nodel);
            end
            FFVAL.TD = FFVAL.TD + demandMatrix(nodef,nodel);
            if numTranOfBest(nodef,nodel)==0
                FFVAL.T0 = FFVAL.T0 + demandMatrix(nodef,nodel);

            elseif numTranOfBest(nodef,nodel)==1
                FFVAL.T1 = FFVAL.T1 + demandMatrix(nodef,nodel);
            elseif numTranOfBest(nodef,nodel)==2
                FFVAL.T2 = FFVAL.T2 + demandMatrix(nodef,nodel);
            end
        end
    end
    FFVAL.ATT = (FFVAL.TIVT + FFVAL.T1*P1 + FFVAL.T2*P2)/FFVAL.TD;
    FFVAL.PT0 = 100*FFVAL.T0/FFVAL.TD;
    FFVAL.PT1 = 100*FFVAL.T1/FFVAL.TD;
    FFVAL.PT2 = 100*FFVAL.T2/FFVAL.TD;
    FFVAL.PT = -omega1*FFVAL.PT0 + omega2*FFVAL.PT1 + omega3*FFVAL.PT2;
    FFVAL.PUD = 100-FFVAL.PT0-FFVAL.PT1-FFVAL.PT2;
    FFVAL.TF = w1*FFVAL.ATT + w2*FFVAL.PT + w3*FFVAL.PUD;
    if(isnan(FFVAL.TF))
        disp('error');
    end
    toc
end