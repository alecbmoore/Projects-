%load the image and store its corresponding matrix 
load dollarblur.m;
D = dollarblur;

%plot of the intial blurred image 
figure(1);
imagesc(D);
colormap(gray);

%parameters to build the matrices A and B 
L = 0.45;
B = build(L);
A = B^25;

%regularization parameter
lambda = 2e-5; 


 for i = 1:500
         
      %obtain ith column of D
      d = D(:,i);
            
      %solve for the ith column of the estimated restored matrix X 
      xest = (A' * A + lambda^2*eye(220)) \ (A' * d);

      %append xest on to X
      X(:,i) = xest;

 end

%plot of the restored image 
 figure(2);
 imagesc(X);
 colormap(gray);


%function to build the matrix B 
function B = build(L)  
  
 B = zeros(220,220);
  
 for i = 1:220
    B(i,i) = 1 - 2*L; 
    
    if i == 220
        break;
    end
    
    B(i+1,i) = L;
    B(i,i+1) = L;
 end  
end 



