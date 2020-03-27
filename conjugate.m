function x = conjugate(A, b)
[r,c] = size(A);

tol = 10^-1;
n = 100;

x = zeros(r,1);
g = b-A*x;
d = g;
for k = 1 :n
    alpha= g'*g/(d'*A*d);
    xprev = x;
    x = x + alpha*d;
    % use successive point difference as stopping condition
    if abs(xprev-x) < tol
        break;
    end
    gp = g;
    g = gp -alpha*A*d;
    beta = g'*g/(gp'*gp);
    d = g+beta*d;
end
if(x - A\b > tol)
    disp('Are the dimensions more than 100?!');
end
%fprintf('%d iterations done\n', k);