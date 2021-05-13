X0=0;
V0=0;
mu=1;
B=-1;
T=1;
M=2;
epsilon=0.01;
theta=1;
tic
[tau N std_err] = mlmc(X0,V0, mu, B, T, M, epsilon,theta)
toc

