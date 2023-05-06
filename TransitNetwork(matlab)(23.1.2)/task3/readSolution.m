function solution =  readSolution(filename, sheetname)
    [num,~] = xlsread(filename,sheetname);
    [nrow, ~] = size(num);
    route_i = 1;
    solution = {};
    temp_route = [];
    for k=1:nrow
        if isnan(num(k,1)) 
            if isnan(num(k-1,1))
                solution{route_i} = temp_route;
                temp_route = [];
                route_i = route_i + 1;
            end 
        else
            temp_route = [temp_route, num(k,1)];
        end 
    end
    solution{route_i} = temp_route;
end