function RE = hanshu(data, noap)

n = size(data, 2);

for i = 1:n
    first(:, i) = data(:, i) - noap;
end

RE = sum(sum(first .^ 2));
end