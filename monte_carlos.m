clear
clc

n = 5e5; %number of samples/realisation
A = rand(n,1); %1st uniform random variable
Z = rand(n,1); %2nd uniform random variable
B = rand(n,1); %3rd uniform random variable

mu_A = 12; %exponential distribution (A)
mu_B = 3;
sigma_B = 0.4;
mu_Z = 0.3;
sigma_Z = 0.15;


A1 = zeros(n,1);
B1 = zeros(n,1);
Z1 = zeros(n,1);

for i = 1:n
    A1(i,1) = expinv(A(i), mu_A);
    B1(i,1) = norminv(B(i), mu_B, sigma_B);
    Z1(i,1) = logninv(Z(i), mu_Z, sigma_Z);

end

y = zeros(n,1);
for i = 1:n
    y(i,1) = A1(i)*(Z1(i)^B1(i));
end

a1 = y > 5;
a2 = y > 10;


%probability of failure when the roof displacement is greater than 5inches
pf1 = (1/n)* sum(a1);
disp(pf1)

%probability of failure when the rood displacement exceeds 10inches
pf2 = (1/n)* sum(a2);
disp(pf2)

%% correlation between A and B
A1 = randn(n,1); %1st uniform random variable
Z1 = randn(n,1); %2nd uniform random variable
B1 = randn(n,1); %3rd uniform random variable
R = [1, 0.75, 0; 0.75, 1, 0; 0,0,1];
L = chol(R, 'lower');

U = [A1';B1';Z1'];
z = L*U;

A2 = zeros(n,1);
B2 = zeros(n,1);
Z2 = zeros(n,1);
for i = 1:n
    A2(i,1) = expinv(z(1, i), mu_A);
    B2(i,1) = norminv(z(2, i), mu_B, sigma_B);
    Z2(i,1) = logninv(z(3, i), mu_Z, sigma_Z);

end

x = zeros(n,1);
for i = 1:n
    x(i,1) = A2(i)*(normcdf(Z2(i))^B2(i));
end

a3 = x > 5;
a4 = x > 10;


%probability of failure when the roof displacement is greater than 5inches
pf3 = (1/n)* sum(a3);
disp(pf3)

%probability of failure when the rood displacement exceeds 10inches
pf4 = (1/n)* sum(a4);
disp(pf4)