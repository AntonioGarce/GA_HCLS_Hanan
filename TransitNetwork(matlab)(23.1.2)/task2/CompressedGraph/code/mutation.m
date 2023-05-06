function [mutatedIndividual, ismutated] = mutation(travelTimeMatrix,solution,compLenMatrix,maxRouLen,minRouLen,numTry)
    ismutated = 0;
    num = 0;
    while(num<numTry)
        num = num+1;
        numroutes = length(solution);
        route_ind = randi(numroutes);
        route = solution{route_ind};
        roulen = length(route);
        nodef_ind = randi(roulen);
        nodel_ind = nodef_ind;
        while(nodel_ind==nodef_ind)
            nodel_ind = randi(roulen);
        end
        if(nodef_ind>nodel_ind)
            [nodel_ind,nodef_ind] = deal(nodef_ind,nodel_ind);
        end
        nodef = route(nodef_ind);
        nodel = route(nodel_ind);
    %     k=2+randi(3);
        k = randi(5);
        [shortestPaths, ~] = kShortestPath(travelTimeMatrix, nodef, nodel, k);
        numpaths = length(shortestPaths);
        if numpaths > 1
            path = shortestPaths{randi(numpaths)};
            pathlen = length(path);
            mutatedIndividual = solution;
            mutatedIndividual{route_ind} = [route(1:nodef_ind-1) path route(nodel_ind+1:end)];
            newpathlen = sum(compLenMatrix(mutatedIndividual{route_ind}));
            if (~isequal(path,route(nodef_ind:nodel_ind)) && isFullSolution(mutatedIndividual,travelTimeMatrix,10) && newpathlen>=minRouLen && newpathlen<=maxRouLen) 
                ismutated = 1;
                break;
            end
        end
    end
end