function compressedNode =  getCompressedNode(prevnd, nd, fulltravelTimeMatrix)
    compressedNode = [nd];
    [ numnodes,~] = size(fulltravelTimeMatrix);
    for neighbour_nd = 1:numnodes
        if isinf(fulltravelTimeMatrix(nd,neighbour_nd))
            continue;
        end
        order_neigh = numnodes - sum(isinf(fulltravelTimeMatrix(neighbour_nd,:)),2) - 1;
        if order_neigh == 2 && neighbour_nd ~= prevnd && neighbour_nd ~= nd
            compressedNode = [compressedNode, getCompressedNode(nd, neighbour_nd, fulltravelTimeMatrix)];
%         elseif order_neigh == 1
%             compressedNode = [neighbour_nd];
        end
    end
end