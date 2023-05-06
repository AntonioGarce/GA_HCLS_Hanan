function path =  getFullRoute(travelTimeMatrix, stnode)
     path=zeros(1,1);
     path(1) = stnode;
     selnode = stnode;
     while(1)
        linkv = travelTimeMatrix(selnode,:);
        linkv(selnode) = Inf;
        for i=1:length(path)
            linkv(path(i))=Inf;
        end
        ind = find(~isinf(linkv));
        if(isempty(ind))
            break;
        end
        prevnode = selnode;
        selnode = ind(randi([1,length(ind)]));
        path = [path selnode]
     end
end