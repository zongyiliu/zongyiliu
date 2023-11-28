ALPHA = 0.01;
BETA = 0.5;
MAXITERS = 1000;
NTTOL  = 1e-8;
x = zeros(n,1);
for iter = 1:MAXITERS
   val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
   d = 1./(1-A*x);
   grad =  A'*d - 1./(1+x) + 1./(1-x);
   hess =  A'*diag(d.^2)*A + diag(1./(1+x).^2 + 1./(1-x).^2);
   v = -hess\grad;
   fprime = grad'*v;
   if abs(fprime) < NTTOL, break; end;
   t = 1;  while ((max(A*(x+t*v)) >= 1) | (max(abs(x+t*v)) >= 1)),
      t = BETA*t;
   end;
   while ( -sum(log(1-A*(x+t*v))) - sum(log(1-(x+t*v).^2)) >  ...
           val + ALPHA*t*fprime )
      t = BETA*t;
   end;
   x = x+t*v;
end;
