function [tau, N, std_err] = mlmc(X0,V0, mu, B, T, M, epsilon,theta)
L=0;
N_initial = 200;
converged = false;

res= zeros(3,1);

while ~converged
    dif = mlmc_sim(X0,V0, mu,theta,  B, T, M, L, N_initial);
    res(1, L+1) = N_initial;
    res(2, L+1) = dif(1);
    res(3, L+1) = dif(2);
    vars = res(3, :) ./ res(1, :) - (res(2, :) ./ res(1, :)).^2;
    
    N= max(N_initial, ceil(2 * sqrt(vars ./ (M.^(0:L))) * sum(sqrt(vars .* (M.^(0:L)))) / epsilon^2));
    
    for l = 0:L
        dN =  - res(1, l+1) + N(l+1);
        if dN>0 
        dif = mlmc_sim(X0,V0, mu,theta, B, T, M,l,dN);
        res(1, l+1) = res(1, l+1) + dN;
        res(2, l+1) = res(2, l+1) + dif(1);
        res(3, l+1) = res(3, l+1) + dif(2);
        end
    end
    
    if L<=3 || T*(M^(-L)) > 0.1*theta
        L=L+1;
        continue;
    end
    converged = max( abs(res(2, L) / res(1, L)) / M, abs(res(2, L+1) / res(1, L+1)) ) < (M-1) * epsilon / sqrt(2);
    L = L + 1;
end
vars = (res(3, :) - res(2, :).^2 ./ res(1, :)) ./ (res(1, :) - 1);
std_err = sqrt(sum(vars ./ res(1, :))); 
tau = sum(res(2, :) ./ res(1, :));
end