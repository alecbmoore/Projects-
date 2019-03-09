function [alpha,Ustar,lambda] = lmbvoptim(n)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%             PLEASE READ              %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% I ended up hard coding the values for u_bar since Fsolve does not 
% allow you to pass a function with more than one input variable. 
% below are the values for alpha and fcost
%
% u_bar1 = sin(pi*x)
% fcost = 0.0925
% alpha = 8.1916
%
% u_bar2 = sin(2*pi*x)
% fcost = 6.3960
% alpha = 5.5650
%
% The plots of the Lagrange multpliers oscillate, which could be   
% because the multiplier lambda scales the constraint equation, which  
% in turn helps provide an optimal solution based on the gradients of 
% the Lagrange function.        
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%intialize solution vector including alpha, u, lambda
U = zeros([2*n+1,1]);

%function call inorder to plot u_bar
[Ustar,u_bar,x] = Lgrad(U);

% this function call solves the nonlinear system of gradients 
% with respect to each variable in our lagrangian function
[Ustar,fval,exitflag] = fsolve(@Lgrad,U);

%final parsing of the vector obtained by Fsolve
alpha = Ustar(end); lambda = Ustar(n+1:2*n); u_final = Ustar(1:n);

%compute the cost of the function 
cost_f = fcost(u_final,u_bar);
cost_f
alpha


% plot of u_bar with solution 
figure(1);
plot(x,u_bar);
hold on 
plot(x,u_final);
title('plot of the opitmal solution vs data vector'); 

% plot of u_bar with lambda 
figure(2); 
plot(x,lambda);
title('plot of vector containing Lagrange multipliers'); 

end

% function to compute the cost of the 
function cost = fcost(u_final,u_bar)

cost = .5 * norm(u_final - u_bar)^2;

end


% function that solves for the grandients of the Lagrange
% function and stores them in to a nonlinear system of 
% equations 
function [F,ubar,x] = Lgrad(U)

    %sizing of vectors 
    n = 99; x = linspace(0,1,n);
    
    %initialize data vector
    ubar = zeros([n,1]);

    %build the data vector   
    for i = 1:n
        ubar(i) = sin(pi*x(i));
    end
    
%     for i = 1:n
%         ubar(i) = sin(2*pi*x(i)).^2;
%     end
    
    %build the tri-diagonal A matrix
    A = zeros(n);
    for i=1:n	
        A(i,i) = 2;	
        if i==n
        	break;
        end
    	A(i,i+1) = -1; A(i+1,i) = -1;
    end
    h = 1/(n+1); 
    A = A*(1/h^2);

    %parse vector in to desired variables that will be solved
    u = U(1:n); lambda = U(n+1:2*n); alpha = U(2*n+1);    
    
    %initialize the vector containing your solution
    F = zeros(2*n + 1,1);
    
    %nonlinear system of equations to be solved
    F(1:n) = u - ubar - A*lambda  -  3 * diag(lambda)*u.^2;
    F(n+1:2*n) = -A*u - u.^3 + alpha*ones(n,1);
    F(2*n+1) = lambda' * ones(n,1);
    
end 

 

