function ord = getOrderOfNode(nd,fulltravelTimeMatrix) %nd:index of node
    travelTimeInfo = fulltravelTimeMatrix(nd,:);
    [ind_neigh] = find(~isinf(travelTimeInfo));
    ord = length(ind_neigh) - 1;
end