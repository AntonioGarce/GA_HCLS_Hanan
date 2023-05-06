degrees = ones(1,149);
degrees(15:116) = 2;
degrees(117:138) = 3;
degrees(139:end) = 4;

cracow = random_graph(149,0.5,535,'sequence',degrees);

for i = 1:149
    for j= 1:149
        if (i~=j) && (cracow(i,j)==0)
            cracow(i,j)=inf;
        end
    end
end

save 'cracow.mat'