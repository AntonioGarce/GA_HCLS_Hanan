function [MemeSol_ind, fitTF] = HCLS(MemePop,fitMatrix,travelTimeMatrix,demandMatrix,maxIters)
    MemeSol_ind = MemePop;
    poplen = length(MemePop);
    for iter = 1:maxIters
        for pop_i = 1:round(poplen/2)
            sol1_ind = randi(poplen);
            sol2_ind = sol1_ind;
            while sol1_ind == sol2_ind
                sol2_ind = randi(poplen);
            end
            sol1 = MemeSol_ind(sol1_ind);
            sol2 = MemeSol_ind(sol2_ind);
            sol1 = sol1{1};
            sol2 = sol2{1};
            [~,roulen] = size(sol1);
            rou1_ind = randi(roulen);
            rou2_ind = randi(roulen);
            newsol1 = sol1;
            newsol2 = sol2;
            newsol1{rou1_ind} = sol2{rou2_ind};
            newsol2{rou2_ind} = sol1{rou1_ind};
            newMemeSol_ind = MemeSol_ind;
            if (isFullSolution(newsol1,travelTimeMatrix,10) && isFullSolution(newsol2,travelTimeMatrix,10)) 
                newffval1 = ff_fun(newsol1,demandMatrix,travelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
                newffval2 = ff_fun(newsol2,demandMatrix,travelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
                oldffval1 = fitMatrix(sol1_ind);
                oldffval2 = fitMatrix(sol2_ind);
%                 oldffval2 = ff_fun(sol2,demandMatrix,travelTimeMatrix,0.4,0.5,0.1,0.6,0.3,0.1,5,10);
                if newffval1.TF<=oldffval1 && newffval2.TF<=oldffval2
                    newMemeSol_ind{sol1_ind} = newsol1;
                    newMemeSol_ind{sol2_ind} = newsol2;
                    fitMatrix(sol1_ind) = newffval1.TF;
                    fitMatrix(sol2_ind) = newffval2.TF;
                end
            end
        end
        MemeSol_ind = newMemeSol_ind;
    end
    fitTF = fitMatrix;
end