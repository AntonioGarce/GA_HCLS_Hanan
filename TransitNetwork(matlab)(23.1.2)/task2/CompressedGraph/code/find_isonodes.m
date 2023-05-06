function isonodes = find_isonodes(travelTimeMatrix)
    [numnodes,~]=size(travelTimeMatrix);
    isonodes = find(sum(isinf(travelTimeMatrix),2)==numnodes-2);
end