function selection = select_individual(fitDistribution)
    ind = find(fitDistribution>rand);
    selection = ind(1)-1;
end