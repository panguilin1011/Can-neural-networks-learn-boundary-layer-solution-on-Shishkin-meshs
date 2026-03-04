function val = F_prime(x, b2, W2, W3)
    val = 0;
    for i = 1:length(W3)
        z = W2(i)*x + b2(i);
        if z > 0
            val = val + W3(i)*W2(i);  % ReLU' = 1
        else
            val = val + 0;            % ReLU' = 0
        end
    end
end
