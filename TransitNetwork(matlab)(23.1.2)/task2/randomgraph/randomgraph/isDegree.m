function res = isDegree(degrees)
    % there are non-positive degrees or their sum is odd
    if not(isempty(find(degrees<=0,1))) || mod(sum(degrees),2)==1
        res = false; 
        return;
    end
    % get the number of nodes
    n=length(degrees);
    % sort the degrees in descend order
    degrees=sort(degrees, 'descend');   
    % check the degrees satisfy the graph condition

    for k=1:n-1
        % sum of degrees for of group which contains k numbers of nodes with largest degrees
        sum_l = sum(degrees(1:k));
        % sum of degrees for other group
        sum_s = sum(min([k*ones(1,n-k);degrees(k+1:n)]));
        % check the graph has proper degree sequences.
        % k*(k-1) is the maximum of sum of degrees for k nodes
        if sum_l > k*(k-1) + sum_s
            res = false; 
            return; 
        end
    
    end

    res = true;
end