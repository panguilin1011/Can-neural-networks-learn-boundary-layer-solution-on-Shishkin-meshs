function costval = cost(F,M,N,b2,W2,W3,beta,gamma)
    costval = gamma*(F{1}^2 + F{N+1}^2) + ...
        beta*(sum(b2.^2) + sum(W2.^2) + sum(W3.^2));
    for j = 2:N
        costval = costval + M{j}^2;
    end
end