function [x,fval,gnorm,niter,nfeval,ngrad] = optim_sample(A, x0, alpha0, tol, imethod, maxit)

%input: 
% A coefficient for the Rosenbrock function
% x0 initial guess vector
% alpha0 initial step size (length)
% tol tolerance criteria for stopping the main iteration
% imethod parameter used to select the optimization method
% imethod = 2 implements steepest descent
% imethod = 1 implements Newton
% for all other integer values implements scaled steepest descent
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

% parameters for line backtracking
c = 0.01; factor = 0.5; nmax = 20; 

iflag = 0;
niter = 0;
nfeval = 0; ngrad = 0;

x0 = x0(:); %column format
f0 = fcost(x0,A); nfeval = nfeval+1;
g0 = fgrad(x0,A); ngrad = ngrad + 1;
gnormk = norm(g0,2);


%trajectory storage 
x = x0;
fval = f0;
gnorm = gnormk;

iplot = 1;
if iplot
   doplot(A)
   plot(x0(1),x0(2),'ro','Linewidth',10) 
   pause(1)
end

while gnormk > tol 
    if imethod == 1 % Newton's method
        D = fhess(x0,A); 
        p0 = -D\g0;
    elseif imethod == 2 % steepest descent method
        p0 = -g0;      
    else % scaled steepest descent
        D = diag(diag(fhess(x0,A)));
        p0 = -D\g0;
    end 
    
    [alpha,iflag] = linebacktrack(alpha0,c,factor,nmax,x0,f0,p0,g0);
   % alpha = alpha0; %implements a constant step size alpha0
    if iflag ~= 0 
        break;
    else    
        
        % move to next closes postion to the solution 
        x0 = x0 + alpha*p0; 
       
        % plot the new position 
        if iplot
           plot(x0(1),x0(2),'go','Linewidth',10) 
           pause(0.5)
        end
        
        %incremment routine iterations 
        niter = niter+1;
        
        %update the current cost of function 
        f0 = fcost(x0,A); nfeval = nfeval+1;
        
        %compute new gradient 
        g0 = fgrad(x0,A); ngrad = ngrad + 1;
        
        %compute current norm of the gradient 
        gnormk = norm(g0,2);
        
        %store trajectory 
        x = [x x0]; fval = [fval; f0]; gnorm = [gnorm; gnormk]; 
        
    
        if niter == maxit 
           fprintf('Max number of iterations reached: %i',niter) 
           break; 
        end 
        
    end
end


fprintf('number of function evaluations: %i', nfeval);
fprintf('number of funciton iterations: %i',niter);


end


function [alpha,iflag] = linebacktrack(alpha0,c,factor,nmax,x0,f0,p0,g0)

%initialize vectors 
alpha = 0:nmax; rho = 0:nmax;
x = 0: nmax; f = 0: nmax;
g = 0: nmax; p = 0: nmax;


%intialize first element of vectors   
alpha(1) = alpha0; rho(1) = factor;
x(1) = x0; f(1) = f0; 
g(1) = g0; p(1) = p0;

%intialize iterator 
k = 0;

%repeat the process below until the optimal point is found
while(norm(g(k))>= 0) 
    
    p(k) = - g(k);
    
    %check if the armijo condition is satisfied 
    while f(x(k) + alpha*p(k)) > f(x(k)) + c*alpha*transpose(p(k))*g(k)
       alpha = rho*alpha;     
    end
    
    % 
    alpha(k) = alpha; 
    k = k+1;
end

%
x(k+1) = x(k) + alpha(k)*p(k);
end

%function that computes the cost of the function
function f = fcost(x,A) 

f = A*(x(2) - x(1)^2)^2 + (1-x(1))^2;
end

%function that computes 
function g = fgrad(x,A)

g = zeros(2,1);
g(1) = -4*A*(x(2)- x(1)^2)*x(1) -2*(1-x(1)); 
g(2) = 2*A*(x(2) - x(1)^2);

end

function H = fhess(x,A)

H = zeros(2,2);

H(1,1) = 12*A*x(1)^2 -4*A*x(2) + 2;
H(1,2) = -4*A*x(1); 
H(2,1) = H(1,2);
H(2,2) = 2*A;

end


function doplot(A)
[x,y] = meshgrid(-2:.1:2, -2:.1:2);
%[x,y] = meshgrid(-5:.1:5, -5.:.1:5);
f1 = (y - x.^2); 
f2 = 1-x ; 
f = A*f1.^2 + f2.^2;
figure(1)
hold on
[c,h]= contour(x,y,f,0:1:20);
clabel(c,h)
%axis square
hold on
plot(1,1,'gs','Linewidth',10)
title(strcat('Rosenbrock function: f(x,y) =', num2str(A), ... 
    '(y-x^2)^2 + (1-x)^2. Optimal point (1,1)'),'fontsize',16)

end

    


    
