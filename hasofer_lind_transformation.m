clear
clc
%IMPLEMENTING OF THE HASOFER-LIND RELIABILITY INDEX
syms o v fy

M = [45*pi/180;150;50]; % mean values for the random variable

coeff_val = [0.1; 0.25; 0.1]; %coefficient of variation

standard_dev = coeff_val.*M; %standard deviation 

D = diag(standard_dev); %diagonal matrix 

D_inv = inv(D); %inverse of the diagonal matrix

R = [1,0,0;0,1,0;0,0,1]; %covariance matrix 

L = chol(R)'; %lower triangular matrix using cholesky decomposition

L_inv = inv(L); %inverse of lower triangular matrix

DL = D*L;

x_values = [45*pi/180;150;50];

%differentiating with respect to the random variables
f = 1 - ((120*v*sin(o))/(1224*fy))^2 - ((120*v*cos(o))/(612*fy))^2;
df_o = diff(f, o);
df_v = diff(f, v);
df_fy = diff(f, fy);




%first iteration
U1 = (L_inv)*(D_inv)*(x_values - M);

g1 = subs(df_o, [o,v,fy], [x_values(1), x_values(2), x_values(3)]);
g2 = subs(df_v, [o,v,fy], [x_values(1), x_values(2), x_values(3)]);
g3 = subs(df_fy, [o,v,fy], [x_values(1), x_values(2), x_values(3)]);

grad_g = [double(g1); double(g2); double(g3)];

grad_h = (D*L)'*grad_g; %grad_h from the slide

norm1 = norm(grad_h); %normalise the grad_h

alpha = -grad_h*(1/norm1);

current_beta = transpose(alpha)*U1;

h1 = 1 - ((120*x_values(2)*sin(x_values(1)))/(1224*x_values(3)))^2 - ((120*x_values(2)*cos(x_values(1)))/(612*x_values(3)))^2;

nU = alpha * (current_beta + (h1/norm1));

nX = (DL*nU) + M;

%storing the output
nX_store = nX; % for the next iteration
x_values_stored = []; %for storing the values

U_store = nU; % for the next iteration 
U_values_stored = []; 
beta = [];
alpha1_store = [];

h2 = h1;

%previous for the loop
prev_beta = 1;

while (abs(prev_beta - current_beta) >= 0.001) || (abs(h2) >= 0.001)
    U2 = U_store;
  
    gg1 = subs(df_o, [o,v,fy], [nX_store(1), nX_store(2), nX_store(3)]);
    gg2 = subs(df_v, [o,v,fy], [nX_store(1), nX_store(2), nX_store(3)]);
    gg3 = subs(df_fy, [o,v,fy], [nX_store(1), nX_store(2), nX_store(3)]);

    prev_beta = current_beta;

    ggrad_g = [double(gg1); double(gg2); double(gg3)];

    ggrad_h = (DL)' * ggrad_g;

    norm2 = norm(ggrad_h);

    alpha1 = -ggrad_h*(1/norm2);

    alpha1_store = [alpha1_store, alpha1];

    current_beta = transpose(alpha1)*U2;

    beta = [beta, current_beta];

    h2 = 1 - ((120*nX_store(2)*sin(nX_store(1)))/(1224*nX_store(3)))^2 - ((120*nX_store(2)*cos(nX_store(1)))/(612*nX_store(3)))^2;

    nU1 = alpha1 *(current_beta + (h2/norm2));

    %storing the value of U
    U_values_stored = [U_values_stored,nU1];

    %using the value of U for the next iteration
    U_store = nU1;

    nX1 = (DL*nU1) + M;

    %using the values of x for the next iteration
    nX_store = nX1;

    %storing values
    
    x_values_stored = [x_values_stored, nX1];    
end

%IMPORTANCE VECTOR AT DESIGN POINT
La = (inv(L))'*alpha1_store(:,6);
La_norm = norm(La);

importance_vector = La/La_norm;
disp(importance_vector)
