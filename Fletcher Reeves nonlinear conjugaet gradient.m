function [x,fval,gnorm,niter,nfeval,ngrad] = hw4(A, x0, alpha0, tol, maxit)


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

%local parameters for line backtracking
c = 0.01; rho = 0.5; nmax = 20; 

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

% Initial search direction is steepest descent
p0 = -g0;

% Create meshgrid & plot starting point
iplot = 1;
if iplot
   doplot(A)
   plot(x0(1),x0(2),'ro','Linewidth',10) 
   pause(1)
end

while gnormk > tol 
    
    % line backtrack to determine alpha
    [alpha,iflag, nfeval] = linebacktrack(alpha0,c,rho,nmax,f0,x0,p0,g0,A,nfeval);
    if iflag ~= 0 
        break;
    else    
        % evaluate new x
        xn = x0 + alpha*p0; 
        
        % plot new x
        if iplot
           plot(xn(1),xn(2),'go','Linewidth',10) 
           pause(0.2)
        end
           
        % By breaking up beta calculation, we avoid having to store g0 
        betaDenom = transpose(g0) * g0;
        
        %re-evaluate f0, g0 with new x0 values
        fn = fcost(xn,A); nfeval = nfeval+1;
        gn = fgrad(xn,A); ngrad = ngrad + 1;
        gnormk = norm(gn,2);
        
        % store trajectory 
        x = [x xn]; fval = [fval; fn]; gnorm = [gnorm; gnormk]; 
        
        % Fletcher-Reeves
        beta = (transpose(gn) * gn) / betaDenom;
        pn = -gn + beta * p0;
        
        % cold restart
        if transpose(pn) * gn >= 0
          p0 = -gn;
          fprintf('Cold restart at: %i',niter)
        else
          p0 = pn;
        end
        
        niter = niter+1;
        if niter == maxit 
           fprintf('Max number of iterations reached: %i',niter) 
           break; 
        end 
        
        %alpha0 = alpha;
        x0 = xn;
        f0 = fn;
        g0 = gn;
        
    end
end

nfeval
ngrad
niter

end


function [alpha,iflag, nfeval] = linebacktrack(alpha0,c,rho,nmax,f0,x0,p0,g0,A, nfeval)

n = 1; i = 1; % indices
alphaK = alpha0;

while fcost((x0 + alphaK(i)*p0),A) > (f0 + c*alphaK(i)*transpose(g0)*p0)
    nfeval = nfeval + 1;
    alphaK(i+1) = alphaK(i)*rho;
    i = i+1;
    
end

n = n+1;
alpha = alphaK(i);
    if n > nmax
        iflag = 1;
    else
        iflag = 0;
    end 

end


function f = fcost(x,A) 

f = A*(x(2) - x(1)^2)^2 + (1-x(1))^2;

end

function g = fgrad(x,A)

g = zeros(2,1);
g(1) = -4*A*(x(2)- x(1)^2)*x(1) -2*(1-x(1)); 
g(2) = 2*A*(x(2) - x(1)^2);

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

    


    
