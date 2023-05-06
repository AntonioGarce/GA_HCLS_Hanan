function output = isFullRoutes(solution)
    route = [];
    numRoute = length(solution);
    if numRoute~=0
        disp('d');
    end
    for k=1:numRoute
        route = [route solution{k}];
    end
    route = unique(route);
    if isequal(route,1:75)
        output = true;
    else
        output = false;
    end
end