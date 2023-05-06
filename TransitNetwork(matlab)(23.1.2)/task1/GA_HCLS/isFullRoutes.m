function output = isFullRoutes(travelTimeMatrix)
    output = 0;
    [numnodes,~] = size(travelTimeMatrix);
    if(sum(sum(isinf(travelTimeMatrix),2))==(numnodes-1)*numnodes)
        output = 1;
    end
end