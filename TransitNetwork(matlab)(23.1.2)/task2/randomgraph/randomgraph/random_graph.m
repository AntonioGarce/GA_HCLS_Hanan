%       N:  number of nodes
% degrees:  degree conditions

%     adj:  adjacency matrix of generated graph

function adj = random_graph(n,degree_p,degree,time_p,time)
    % initialize the adjacent matrix
    adj=zeros(n); 

    if sum(degree_p) ~= 100 
        disp("Sum of degree_percentage should be 100%");
        return;
    end

    if sum(time_p )~= 100
        disp("Sum of time_percentage should be 100%");
        return;
    end
    
    % distribution probability of degrees
    degree_dist = zeros(1,length(degree_p)+1); 
    % number of degrees
    degree_num = zeros(1,length(degree_p)+1);
    
    for k = 1:length(degree_p)
        degree_dist(k+1) = degree_dist(k) + degree_p(k);
        degree_num(k+1) = floor(degree_dist(k+1)*n/100);
        degrees(degree_num(k)+1:degree_num(k+1)) = degree(k);
    end
    % number of edges
    num_edges = sum(degree.*(degree_num(2:end)-degree_num(1:length(degree_p))))/2;
    % distribution probability of time
    time_dist = zeros(1,length(time_p)+1);
    % time of edges
    time_num = zeros(1,length(time_p)+1);
    % 
    for k=1:length(time_p)
        time_dist(k+1) = time_dist(k) + time_p(k);
        time_num(k+1) = floor(time_dist(k+1)*num_edges/100)- sum(time_num(1:k));
    end
    
    % if degrees don't satisfy the graph condition
    if not(isDegree(degrees))
        fprintf('degrees is not a graphic sequence - select a different sequence\n'); 
        return; 
    end
    % old sum of degrees
    old_sum = 0;
    % counter for degrees have not been updating.(not updating degrees
    % means that the graph is not updated)
    cnt=0;
    degrees_org = degrees;
    while sum(degrees)>0 
      % graph is not updated for 100 times, restart algorithm
      if cnt>100
          degrees = degrees_org; 
          adj=zeros(length(degrees_org)); 
          cnt=0; 
      end
      new_sum = sum(degrees);
      % graph is not updated
      if old_sum==new_sum
          cnt = cnt+1; 
      end
      % graph is updated
      if old_sum~=new_sum
          cnt=0; 
      end
      % update old sum 
      old_sum = new_sum;
      % index of node which has maximum remain degrees
      [~,n1] = max(degrees);
      ind = find(degrees>0);
      % select node which has remain degree
      n2 = ind(randi(length(ind)));
      % if the same nodes are selected or the edge between nodes are
      % already existed
      if n1==n2 || adj(n1,n2)>0 
          continue; 
      end
      % else create the edges and updated the degrees
      findflag = 1;
      while(findflag)
          tt = randi(10000)/100;
          for k = 1:length(time_p)
            if tt<time_dist(k+1) && time_num(k+1)~=0
                traveltime = time(k);
                time_num(k+1) = time_num(k+1)-1;
                findflag = 0;
                break;
            end
          end
      end
      adj(n1,n2)=traveltime; adj(n2,n1)=traveltime;
      degrees(n1) = degrees(n1) - 1;
      degrees(n2) = degrees(n2) - 1;
    end
    
end  