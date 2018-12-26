function [cor_n] = findrange(cor_n,Range)
for i = 1:length(cor_n)
    for j = 1:size(cor_n,2)
        if cor_n(i,j) <= Range(j,1) , cor_n(i,j) = Range(j,1); end
        if cor_n(i,j) >= Range(j,2) , cor_n(i,j) = Range(j,2); end
    end
end