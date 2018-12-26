function [cor_n]=ffa_move(cor_n,lightn,cor_0,light0,alpha,gamma,Range)
ni = size(cor_n,2); nj = size(cor_0,2);
distances = dist(cor_n');
[k,h,l]=find(distances);
distances = (distances - min(l))/max(max(distances));

for i=1:ni
for j=1:nj
	r = distances(i,j); % sqrt((cor_n(i)-cor_0(j))^2 + (yn(i)-y0(j))^2);
	if lightn(i) > light0(j) , beta0 = 1;  beta = beta0*exp(-gamma*r.^2);
		cor_n(i,:) = cor_n(j,:) .* (1-beta) + cor_0(i,:).*beta + alpha.*(rand-0.5);
%		yn(i) = yn(i) .* (1-beta) + y0(j).*beta + alpha.*(rand-0.5);
	end
end
end
[cor_n] = findrange(cor_n,Range);

