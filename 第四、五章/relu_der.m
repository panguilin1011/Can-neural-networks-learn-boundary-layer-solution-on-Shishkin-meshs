function y = relu_der(x)
    y = floor(heaviside(x));
end