function [costf,p,x,l] = optim_model(x0,y_m,m,h)

%function handle for first order optimality conditions 
fun = @(p) grad(x0,p,m,h,y_m);

%solve nonlinear system provided by first order optimality conditions
[p,fval,exitflag] = fsolve(fun,-1.2); 

% **************************************
% ****  OUTPUT  ***********************
% **************************************      

%state vector 
x = model(x0,p,h,m);
x

%vector of lagrange multipliers 
l = l_multipliers(x0,p,m,h,y_m);
l

%compute the cost of the function   
cost_f = fcost(x(m),y_m);
cost_f

%omptimal parameter P given by solving first order optimality conditions  
p

%%%%%%%%%% Results for p %%%%%%%%%%%%
%%%%% when x0 = 2 , p = -1.2186 %%%%%
%%%%% when x0 = 0.5, p = 2.3700 %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% **************************************
% **************************************
% **************************************

    %check to see which plot to use based on what user provides
    if( x0 == 0.5 ) 
       %plot of state vector based on its iterations 
       figure(1)
       plot(x);
       title('time series plot of discrete states when x0 = 0.5');

       %plot of Langrange Multipliers based on state iterations 
       figure(2)
       plot(l);
       title('Plot of Langrange multipliers when x0 = 0.5');
    else 
       %plot of state vector based on its iterations 
       figure(1)
       plot(x);
       title('time series plot of discrete states when x0 = 2');

       %plot of Langrange Multipliers based on state iterations 
       figure(2)
       plot(l);
       title('Plot of Langrange multipliers when x0 = 2');
    end 
    
end

%function that obtains each state 
function x = model(x0,p,h,m) 
  
  %initialize state vector 
  x = zeros(m+1,1);
  
  %initial state
  x(1) = x0;
  
  %obtain all other states in terms of initial state 
  for j = 2:m+1  
      
      %model equation for each discrete state
      x(j) = x(j-1) + h*(x(j-1)*(1-x(j-1)) + p);
      
  end 
  
end 

%function that obtains values for lambda
function l = l_multipliers(x0,p,m,h,y_m)
    
    %obtain state vector 
    X = model(x0,p,h,m);
    
    %initialize lambda values 
    l = zeros(m+1,1);
    
    %lambda corresponding to final state 
    l(m+1) = X(m+1) - y_m;
    
    %obtain all other lambdas using gradient equation
    for j = m : -1 : 1
        
        %first order condition on lambda 
        l(j) = l(j+1)*(1 + h - 2*h*X(j));  
    
    end
    
end

%function builds nonlinear system from first or optimality conditions 
function F = grad(x0,p,m,h,y_m)
    
    %obtain lagrange multipliers 
    l = l_multipliers(x0,p,m,h,y_m);

    %gradient with respect to the parameter P
    F = h * sum(l);
    
end 


%cost functional (cost of final state)
function f = fcost(x,y_m)
    f = .5 * (x(end) - y_m)^2 ;  
end 