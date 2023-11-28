% Generate random number for test, which would also be used in the next m. file

randn('state',1);
m=200;
n=100;

% For the gradient and Hessian for the function
% f(x) = log(sum(exp(a_i' * x))),
% where a_i' is row i in the matrix A


ALPHA = 0.01;
BETA = 0.5;
MAXITERS = 1000;
NTTOL = 1e-8;
GRADTOL = 1e-3;

% Do it in CVX
% randn('seed',0); m=500; n=100; A = randn(m,n);
% cvx_begin
% variable x(n)
% minimize(log_sum_exp(A*x))
% cvx_end


% Generate random problem
A = randn(m,n);
% gradient method

% function [x,fval,tvec,k]=q930gradient(A)
% Then I used fval to represent values, and tvec for tvec

fval = []; tvec = [];
x = zeros(n,1);
for iter = 1:MAXITERS
val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
fval = [fval, val];
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
tvec = [tvec,t];
end;

% Generate Figures
% From figure 1, we can see that the relationship between function values and step lengths versus iteration number.
% The function value is gradrally decreasing as the step increases

figure(1)
semilogy([0:(length(fval)-2)], fval(1:length(fval)-1)-optval, '-');
xlabel('x'); ylabel('z');

% From figure 2, we can see that the gradient is kept roughly at several fixed values, 0.001, 0.002, 0.004, and 0.008

figure(2)
plot([1:length(tvec)], tvec, ':',[1:length(tvec)], tvec, 'o');
xlabel('x'); ylabel('z');


% Figures show an example with m = 200, n = 100.
