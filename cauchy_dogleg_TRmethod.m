function [x,fval,gnorm,niter,nfeval,ngrad] = ...
          cauchy_dogleg_TRmethod(A,x0,tol,imethod,del0,delmax,eta,maxit)

%input: 
% A coefficient for the Rosenbrock function
% x0 initial guess vector
% tol tolerance criteria for stopping the main iteration
% imethod integer parameter used to select the optimization method
% imethod = 1 implements dogleg method
% imethod ~= 1 implements Cauchy point
% del0 initial trust-region radius
% delmax maximum trust-region radius
% eta TR parameter for sufficient reduction 
% maxit the maximum number of iterations allowed

%output
% x the approximate solution to the optimization 
% fval the value of the cost function at x, that is f(x)
% gnorm the norm of the gradient at x
%REMARK: this script stores the whole sequence of iterates, that is 
% x = [x0, x1, ..., xlast], same for fval and gnorm. 
% niter the number of iterations required to reach the optimal solution
% nfeval the number of function evaluations
% ngrad the number of gradient evaluations

%delmax = 1; 
%del = 0.5; eta = 0.01; 
%maxit = 10000; 

%local parameters for TR radius adjustment
c1 = 0.25; c2 = 0.75;

%machine epsilon
uround = eps; 
%uround = 2.2204e-16;

niter = 0; nfeval = 0; ngrad = 0;

x0 = x0(:);
f0 = fcost(x0,A); nfeval = nfeval+1; 
g0 = fgrad(x0,A); ngrad = ngrad+1;
B  = fhess(x0,A);
gnorm0 = norm(g0,2);

fval = [f0]; % aray storing the cost function evolution
gnorm = [gnorm0]; %aray storing the norm of the gradient 
x = [x0(:)]; %array storing the path to optimal solution

del = del0;
while gnorm0 > tol && niter < maxit
    
    if imethod == 1
       p = dogleg(g0,B,del);
    else
       p = cauchy(B,g0,del,gnorm0);
    end
    xnew = x0+p;
    fnew = fcost(xnew,A); nfeval = nfeval+1;
    rho = (f0-fnew)/(-g0'*p - 0.5*p'*(B*p));
    
    if rho < c1
        del = c1*del;
    elseif rho > c2 && abs(norm(p,2) - del) < 10*uround
        del = min(delmax,2*del);
    end
    
    if rho > eta
        x0 = xnew;
        f0 = fnew; 
        g0 = fgrad(x0,A); gnorm0 = norm(g0,2); ngrad = ngrad+1;
        B = fhess(x0,A);
        fval = [fval; f0];
        gnorm = [gnorm; gnorm0];
        x = [x x0(:)];
    end
    
    niter = niter+1;
end


niter 
nfeval
ngrad

end


function p = dogleg(g0,B,del,nfeval)

%your job to write it
  
  
end

function p = cauchy(B,g0,del,gnorm0)

%your job to write it

if transpose(g0)*B*g0 <= 0
    tau = 1;    
else  
    tau = min((gnorm0.^3) / del*transpose(g0)*B*g0,1); 
end
p = -tau*(del/gnorm0)*g0;

end

function f = fcost(x,A) 

f = A*(x(2) - x(1)^2)^2 + (1-x(1))^2;

end

function g = fgrad(x,A)

g(1) = -4*A*(x(2)- x(1)^2)*x(1) -2*(1-x(1)); 
g(2) = 2*A*(x(2) - x(1)^2);
g=g(:);

end

function H = fhess(x,A)

H = zeros(2,2);

H(1,1) = 12*A*x(1)^2 -4*A*x(2) + 2;
H(1,2) = -4*A*x(1); 
H(2,1) = H(1,2);
H(2,2) = 2*A;

end
