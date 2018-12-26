function [x,lightn]=init_ffa(n,Range)
diff = Range(1,2) - Range(1,1);
x = rand(n,size(Range,1)) * diff + Range(1,1);
lightn = zeros(size(n));