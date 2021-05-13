function [dif,fine] =mlmc_sim(X0,V0, mu,theta,B,T,M,l,N)
    Kf = M^l; % number of points on finer level
    Kc= M^(l-1); % number of points on coarser level
    hf=T/Kf; % time step on finer level
    hc=T/Kc; % time step on coarser level
    sqrthf = sqrt(hf);
    sumf = 0; %sum of simulated first passage times
    sumdif=0;
    squaresumf = 0;%sum of squares of simulated first passage times
    squaresumdif=0;
    Xf = X0 * ones(1, N); % current value of the X (finer path)
    Xc = Xf; % current value of the X (coarser path)
    Vf = V0 * ones(1,N);% current value of the X (finer path)
    Vc = Vf;% current value of the X (coarser path)
    tauf = T * ones(1, N); % estimated first pasage time of the process (finer path)
    tauc = T * ones(1, N); %  estimated first pasage time of the process (coarser path)
    if l==0 %if initial level
    
        dW = sqrthf * randn(1, N);
        for j=1:N
            tmp = Xf(j);
            Xf(j) = Xf(j) + theta ^ (-1) * Vf(j) * hf;
            Vf(j) = Vf(j) - theta^(-2) * Vf(j) * hf - theta^(-1) * mu * PHI_prime(tmp) * hf +(2^0.5)* theta^(-1) * dW(j);
            
        end
        tauf(Xf < B) = T/2;
    
    else
    
        for i = 1:Kc
            dW_f = sqrthf * randn(M,N);
            for m=1:M 
                for j=1:N
                    tmp = Xf(j);
                    Xf(j) = Xf(j) + theta ^ (-1) * Vf(j) * hf;
                    Vf(j) = Vf(j) - theta^(-2) * Vf(j) * hf - theta^(-1) * mu * PHI_prime(tmp) * hf +(2^0.5)* theta^(-1) * dW_f(m,j);
                end
                tauf(Xf < B) = min(((i - 1) * M + m- 0.5) * hf, tauf(Xf < B));
            end
            
            dW_c = sum(dW_f);
            for j=1:N
                tmp = Xc(j);
                Xc(j) = Xc(j) + theta ^ (-1) * Vc(j) * hc;
                Vc(j) = Vc(j) - theta^(-2) * Vc(j) * hc - theta^(-1) * mu * PHI_prime(tmp) * hc +(2^0.5)* theta^(-1) * dW_c(j);
                
            end
            tauc(Xc < B) = min((i - 0.5) * hc, tauc(Xc<B));
        end
    
    end
    sumf = sumf + sum(tauf);
    sumdif = sumdif + sum(tauf - tauc);
    squaresumdif = squaresumdif + sum((tauf - tauc).^2);
    squaresumf = squaresumf + sum(tauf.^2);
    fine = [sumf, squaresumf]; 
    dif = [sumdif, squaresumdif];
    if l==0
        dif = fine;
    end 
end

  