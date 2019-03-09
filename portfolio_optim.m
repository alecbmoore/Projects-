%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%             PLEASE READ              %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%
%Results:
%x = (3.8900,2.3624,3.3509)
%lam1 = 1.1051
%lam2 5.4689
%optimal_LR = 1.7304 (otherwise known as l)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function portfolio_optim(n)         

    %The linear space where our desired minimum point will lie 
    l = linspace(1.2,2.5,n);
    
    %intialize solution vector 
    x = zeros(3,1);
        
    %initialize covarinace matrix 
    Q = covarmtrx();
    
    %initialize rate of return vector 
    rbar = rate_of_return();
    
    %abreviated terms for solutions for lambda 
    a = .5 * transpose(ones(3,1)) * inv(Q) * ones(3,1);
    b = .5 * transpose(ones(3,1)) * inv(Q) * rbar;
    c = .5 * transpose(rbar) * inv(Q) * rbar;
    
    
    for i = 1 : n
        
        %lagrange multipliers
        lam1 = ( b * l(i) - c ) / (b^2 - a*c);
        lam2 = (b - a * l(i)) / (b^2 - a*c);
        
        %optimal solution to x' Q x 
        x = (.5 * (inv(Q)) * ones(3,1) * lam1) + (.5 * (inv(Q)) * rbar * lam2);
        
        %variance of portfolio ( i.e sigma^2 = x' Q x )
        portfolio_variance(i) = transpose(x) * Q * x;
        
    end   
    
    %plot of minimiztion problem 
    plot(portfolio_variance,l);
        
    %minimum point on x' Q x  
    minl = find(portfolio_variance == min(portfolio_variance));
    
    %optimal level of return 
    optimal_LR = l(minl);
        
    %optimal lagrange multipliers 
    lam1 = ( b * optimal_LR - c ) / (b^2 - a*c);
    lam2 = b - a * optimal_LR / (b^2 - a*c);
        
    % optimal solution 
    x = (.5 * (inv(Q)) * ones(3,1) * lam1) + (.5 * (inv(Q)) * rbar * lam2);
        
    %print results
    x
    lam1
    lam2
    optimal_LR
end 

%Function creating the covariance matrix 
function Q = covarmtrx() 
    Q = zeros(3,3);
    
    Q(1,1) = 1; Q(1,2) = 1/3; Q(1,3) = -1/3;
    
    Q(2,1) = 1/3; Q(2,2) = 2; Q(2,3) = 0;

    Q(3,1) = -1/3; Q(3,2) = 0; Q(3,3) = 3;
end 

%function creating the mean/ expected rates of return data vector  
function R = rate_of_return()

    R = zeros(3,1);
    
    R(1,1) = 1.1;
    R(2,1) = 2;
    R(3,1) = 3;
end 