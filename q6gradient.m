ALPHA = 0.01;
BETA = 0.5;
MAXITERS = 1000;
GRADTOL = 1e-3;
x = zeros(n,1);
for iter = 1:MAXITERS
   val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
   grad =  A’*(1./(1-A*x)) - 1./(1+x) + 1./(1-x);
   if norm(grad) < GRADTOL, break; end;
   v = -grad;
       fprime = grad’*v;
       t = 1;  while ((max(A*(x+t*v)) >= 1) | (max(abs(x+t*v)) >= 1)),
          t = BETA*t;
       end;
       while ( -sum(log(1-A*(x+t*v))) - sum(log(1-(x+t*v).^2)) >  ...
               val + ALPHA*t*fprime )
          t = BETA*t;
       end;
       x = x+t*v;
    end;
