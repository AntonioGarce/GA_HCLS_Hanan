function newParent = crossover(parent,compLenMatrix,maxRouLen,minRouLen,numTry)
    newParent = parent;
    num=0;
    while(num<numTry)
        num=num+1;
        parRoute_ind1 =randi(length(parent));
        parRoute_ind2 = parRoute_ind1;
        while(parRoute_ind2==parRoute_ind1)
            parRoute_ind2 = randi(length(parent));
        end
        parRoute1 = parent{parRoute_ind1};
        parRoute2 = parent{parRoute_ind2};
        [interStop, indL, indR]= intersect(parRoute1,parRoute2);
        if(~isempty(interStop))
            indind = randi(length(indL));
            ind1 = indL(indind);
            ind2 = indR(indind);
            newRoute1 = [parRoute1(1:ind1),parRoute2(ind2+1:end)];
            newRoute2 = [parRoute2(1:ind2),parRoute1(ind1+1:end)]; %end?
            newRouLen1 = sum(compLenMatrix(newRoute1));
            newRouLen2 = sum(compLenMatrix(newRoute2));
            if(max(newRouLen1,newRouLen2)<=maxRouLen && min(newRouLen1,newRouLen2)>=minRouLen)
                newParent{parRoute_ind1} = newRoute1;
                newParent{parRoute_ind2} = newRoute2;
                break;
            end
        end
    end
end