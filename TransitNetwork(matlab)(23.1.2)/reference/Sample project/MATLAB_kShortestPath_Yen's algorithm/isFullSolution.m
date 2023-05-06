function isfull = isFullSolution(solution,travelTimeMatrix,maxTransfer)
    numroutespersol = length(solution);
    optTraveTime = travelTimeMatrix;
    [numnodes, ~] = size(travelTimeMatrix);
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
                if(travelTime<optTraveTime(nodef,nodel))    
                    optTraveTime(nodef,nodel) = travelTime;                       %set travel time to the travelTime matrix
                    optTraveTime(nodel,nodef) = travelTime;                       %set travel time to the travelTime matrix
                end
                prevnode = nodel;                                             %update prevnode to reduce calcuation effort
            end
        end
    end
    % get optimal solution using bellman's opimility principle
    numTransfer = 0;
    isfull = 1;
    while(sum(sum(isinf(optTraveTime)))~=0)
        numTransfer = numTransfer + 1;
        for nodef = 1:numnodes                      %first node
            for nodel = nodef+1:numnodes            %last node
                for nodestop = 1:numnodes           %intermediate node which is intersect of between routes
                    if(nodestop==nodef || nodestop==nodel)
                        travelTime = inf;
                    else
                        travelTime = optTraveTime(nodef,nodestop)+optTraveTime(nodestop,nodel);
                    end
                    if(travelTime<optTraveTime(nodef,nodel))
                        optTraveTime(nodef,nodel) = travelTime;
                        optTraveTime(nodel,nodef) = travelTime;
                    end
                end
            end
        end
        if numTransfer>=maxTransfer
            isfull = 0;
            break;
        end
    end
end