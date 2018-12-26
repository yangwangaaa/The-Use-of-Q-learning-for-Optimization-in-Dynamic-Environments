function [best] = deq(instr)
if nargin < 1 , instr=[12 50]; end
n=instr(1); maxgeneration=instr(2);
rand('state',0);
% str1 = '-sin(x)*(sin(x^2/3.14159))^20';
% str2 = '-sin(y)*(sin(2*y^2/3.14159))^20';

str1 = 'exp(-(x-4)^2-(y-4)^2)+exp(-(x+4)^2-(y-4)^2)';
str2 = '+2*exp(-x^2-(y+4)^2) + 2*exp(-x^2-y^2)';
funcstr = strcat(str1,str2);
f=vectorize(inline(funcstr));
range = [-5 5 -5 5];
% range = [0 4 0 4];
alpha = 0.2;
gamma = 1;

ngrid = 100;
dx = (range(2) - range(1))/ngrid;
dy = (range(4) - range(3))/ngrid;
[x,y] = meshgrid(range(1):dx:range(2),...
		 range(3):dy:range(4));
z=f(x,y);
figure(1); surface(x,y,z);

[xn,yn,Lighten] = init_ffa(n,range);
figure(2);
for i=1:maxgeneration
    contour(x,y,z,15); hold on;
    zn = f(xn,yn);
    [lightn,index]=sort(zn,2, 'descend');
    xn =xn(index); yn = yn(index);
    x0 = xn; y0 = yn; light0 = lightn;
    plot(xn,yn,'.','markersize',10,'markerfacecolor','g');axis(range);
    [xn,yn]=ffa_move(xn,yn,lightn,x0,y0,light0,alpha,gamma,range);
    drawnow;
    hold off;
end
best(:,1) =x0'; best(:,2) = y0'; best(:,3) = light0';

function [xn,yn,lightn]=init_ffa(n,range)
xrange = range(2) - range(1);
yrange = range(4) - range(3);
xn=rand(1,n)*xrange+range(1);
yn=rand(1,n)*yrange+range(3);
lightn = zeros(size(yn));

function [xn,yn]=ffa_move(xn,yn,lightn,x0,y0,...
	light0,alpha,gamma,range)
ni = size(yn,2); nj = size(y0,2);
for i=1:ni
for j=1:nj
	r = sqrt((xn(i)-x0(j))^2 + (yn(i)-y0(j))^2);
	if lightn(i) < light0(j) , beta0 = 1;  beta = beta0*exp(-gamma*r.^2);
		xn(i) = xn(i) .* (1-beta) + x0(j).*beta + alpha.*(rand-0.5);
		yn(i) = yn(i) .* (1-beta) + y0(j).*beta + alpha.*(rand-0.5);
	end
end
end
[xn,yn] = findrange(xn,yn,range);

function [xn,yn] = findrange(xn,yn,range)
for i = 1:length(yn)
	if xn(i) <= range(1) , xn(i) = range(1); end
	if xn(i) >= range(2) , xn(i) = range(2); end
	if yn(i) <= range(3) , yn(i) = range(3); end
	if yn(i) >= range(4) , yn(i) = range(4); end
end
	
