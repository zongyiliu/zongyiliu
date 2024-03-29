% Newton method
% function [x,fval,tvec,k]=q930gradient(A)

fval = []; tvec = [];
x = zeros(n,1);
for iter = 1:MAXITERS

val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
fval = [fval, val];
d = 1./(1-A*x);
grad = A'*d - 1./(1+x) + 1./(1-x);
hess = A'*diag(d.^2)*A + diag(1./(1+x).^2 + 1./(1-x).^2);
v = -hess\grad;
fprime = grad'*v
if abs(fprime) < NTTOL, break; end;
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
optval = fval(length(fval));

% Generate figures

% Here we can see that the function value decreases quite rapidly compared with gradient descent method
figure(3)
semilogy([0:(length(fval)-2)], fval(1:length(fval)-1)-optval, '-', ...
[0:(length(fval)-2)], fval(1:length(fval)-1)-optval, 'o');
xlabel('x'); ylabel('z');

% Each step it roughly takes value 1
figure(4)
plot([1:length(tvec)], tvec, '-', [1:length(tvec)], tvec, 'o');
axis([0, length(tvec), 0, 1.1]);
xlabel('x'); ylabel('z');


% in the figures show the function values and step lengths versus iteration number for the same example.
% Here we used alpha = 0.01 and beta = 0.5, as the same as before
