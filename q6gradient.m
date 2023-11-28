% Generate random number for test, which would also be used in the next m. file

randn('state',1);
m=200;
n=100;


ALPHA = 0.01;
BETA = 0.5;
MAXITERS = 1000;
NTTOL = 1e-8;
GRADTOL = 1e-3;

% generate random problem
A = randn(m,n);
% gradient method
vals = []; steps = [];
x = zeros(n,1);
for iter = 1:MAXITERS
val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
vals = [vals, val];
d = 1./(1-A*x);
grad = A'*d - 1./(1+x) + 1./(1-x);
v = -grad;
fprime = grad'*v;
norm(grad)
if norm(grad) < GRADTOL, break; end;
t = 1;
while ((max(A*(x+t*v)) >= 1) | (max(abs(x+t*v)) >= 1)),
t = BETA*t;
end;
while ( -sum(log(1-A*(x+t*v))) - sum(log(1-(x+t*v).^2)) > ...
val + ALPHA*t*fprime )
t = BETA*t;
end;
x = x+t*v;
steps = [steps,t];
end;
figure(1)
semilogy([0:(length(vals)-2)], vals(1:length(vals)-1)-optval, '-');
xlabel('x'); ylabel('z');
figure(2)
plot([1:length(steps)], steps, ':',[1:length(steps)], steps, 'o');
xlabel('x'); ylabel('z');


% The figures show the function values and step lengths versus iteration number for an example with m = 200, n = 100.
