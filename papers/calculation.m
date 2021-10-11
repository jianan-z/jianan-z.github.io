% function [Bpp, time]=calculation(M,K,lambda,mu,LoopCount)
% if nargin < 5
%     LoopCount=1
% end
clear; M = 600; K = 450; lambda = 1; mu = 1; LoopCount = 1;

t1 = cputime;
for loops=1:LoopCount
Aa = zeros(1, K+1);
Ab = ones(1, K+1);
Ac = zeros(1, K+1);
Rs = zeros(K + 1, M - K);

% calculate R^{M-K} = A_0^{M-K-1} (I - A_1^{M-K})^{-1}.
% only applies to the case K>1 (reasonable when M and K are large) to reduce computation time

    n = M - K;
    den = (M - n)*lambda + n*mu;
    Ac(1) = (M - n)*lambda/den; % first row of A1    
    for i = 2:K % 2 - (K) rows of A1
        ii = i - 1;
        den = (M - n - ii)*lambda + (n + ii)*mu;
        Ac(i) = (M - n - ii)*lambda/den;
        Aa(i) = ii*mu/den;
    end
    den = (M - n - K)*lambda + (n + K)*mu;
    Aa(K+1) = K*mu/den;
    Aa = -Aa; Ac = -Ac;  %Ac, Ab, Aa: A_1^n   
    for i = 1:K
        Ab(i+1) = Ab(i+1) - Ac(i)/Ab(i)*Aa(i+1);
    end

    rev = [zeros(1,K),1];
    rev(K) = - Aa(K+1)/Ab(K+1)*rev(K+1);
    for i = K:(-1):2
        rev(i-1) = -Aa(i)/Ab(i)*rev(i);
    end
    A0 = (M - K - n + 1)*lambda / ((M - K - n + 1)*lambda + (n - 1)*mu + K*mu); %A_0^{M-K+1}
    Rs(:,n) = A0 * rev./Ab; %R^n = A_0^{n-1} (I - A_1^n - R^{n+1}A_2^{n+1})^{-1}

% calculate R^{M-K-1} ... R^0
for n = (M - K - 1):-1:1
    if n == 35
        n
    end
    den = (M - n)*lambda + n*mu;
    Ac(1) = (M - n)*lambda/den; % first row of A1 
    Ab = ones(1,K+1);
    for i = 2:K % 2 - (K)th rows of A1
        ii = i - 1;
        den = (M - n - ii)*lambda + (n + ii)*mu;
        Ac(i) = (M - n - ii)*lambda/den;
        Aa(i) = ii*mu/den;
    end
    den = (M - n - K)*lambda + (n + K)*mu;
    Aa(K+1) = K*mu/den;
    Aa = -Aa; Ac = -Ac;  %Ac, Ab, Aa: A_1^n
    A2 = (n + 1)*mu./((M - n - 1)*lambda + (n + 1)*mu + [0:K]*(mu - lambda)); %A_2^{n+1}
    R = -Rs(:,n+1)'.*A2;  %-R^{N+1} * A_2^{n+1}
    Aa(K+1) = Aa(K+1) + R(K);
    Ab(K+1) = Ab(K+1) + R(K+1);
    for i = 1:K-1
        Ab(i+1) = Ab(i+1) - Ac(i)/Ab(i)*Aa(i+1);
        R(i+1) = R(i+1) - Ac(i)/Ab(i)*R(i);
    end
    Aa(K+1) = Aa(K+1) - Ac(K-1)/Ab(K-1)*R(K-1);
    Ab(K+1) = Ab(K+1) - Ac(K)/Ab(K)*Aa(K+1);
    
    rev = [zeros(1,K),1];
    rev(K) = - Aa(K+1)/Ab(K+1)*rev(K+1);
    for i = K:(-1):2
        rev(i-1) = -Aa(i)/Ab(i)*rev(i) - R(i-1)/Ab(K+1)*rev(K+1);
    end
    A0 = (M - K - n + 1)*lambda / ((M - K - n + 1)*lambda + (n - 1)*mu + K*mu);
    Rs(:,n) = A0 * rev./Ab; %R^n = A_0^{n-1} (I - A_1^n - R^{n+1}A_2^{n+1})^{-1}
    
%     if Ab(K+1) < 0.01
%         Ab(K+1)
%         n
%         Rs(:,n);
%     end
end

%calculate R_0
%     A1 = eye(K+1); 
%     n = 0;
%     den = M*lambda;
    Ac(1) = -1; % first row of A1    
    for i = 2:K % 2 - (K) rows of A1
        ii = i - 1;
        den = (M - ii)*lambda + ii*mu;
        Ac(i) = -(M - ii)*lambda/den;
        Aa(i) = -ii*mu/den;
    end
    den = (M - K)*lambda + K*mu;
    Aa(K+1) = -K*mu/den;
    A2 = mu./((M - 1)*lambda + mu + [0:K]*(mu - lambda)); %A_2^{n+1}
    R = -Rs(:,1)'.*A2;  %-R^{N+1} * A_2^{n+1}
    AK = R;  %AK: the last row
    AK(K) = AK(K) + Aa(K+1);  AK(K+1) = AK(K+1) + 1;

%    [U,S,V] = svd(A1'); %The last column of V spans the nullspace of A1'
%    pi1 = V(:,end);
%
    pi = zeros(1,K+1);
    pi(K+1)=10^-100;
    pi(K)  =-(pi(K+1)*AK(K+1))/Ac(K);
    pi(K-1)=-(pi(K+1)*AK(K) + pi(K)) / Ac(K-1);
    for i=K-1:-1:2
	pi(i-1)=-(pi(K+1)*AK(i) + pi(i+1)*Aa(i+1) + pi(i))/Ac(i-1);
    end
%    keyboard

    p = zeros(K+1, M-K+1);
    p(:,1) = (pi/sum(pi))';
for i = 1:(M-K)
    p(:,i+1) = p(K+1,i)*Rs(:,i); %calculate pi_i
end
p = p/sum(sum(p));

% calculate blocking probability using discrete time Markov chain

Tbd2 = 0; Tcd2 = 0; 
for j = 1:(M-K+1)
    for i = 1:(K)        
        Tcd2 = Tcd2 + p(i,j)*(M - i - j + 2)*lambda/((M-i-j+2)*lambda + (i+ j-2)*mu);
    end
end
for j = 1:(M-K+1)
    for i = K+1        
        Tbd2 = Tbd2 + p(i,j)*(M - i - j + 2)*lambda/((M-i-j+2)*lambda + (i+ j-2)*mu);
    end
end
Bpp = Tbd2/(Tbd2 + Tcd2)

% % calculate blocking probability using continuous time Markov chain
% for i = 1:(M-K+1)  %calculating pi in ctmc
%     p(i,:) = p(i,:)./((M - i + 1)*lambda + (i - 1)*mu + (mu - lambda)*[0:K]);
% end
% p = p/sum(sum(p));
% % To1 = 0; Tc1 = 0;
% % for j = 1:(K + 1)
% %     for i = 1:(M - K + 1)
% %         To1 = To1 + p(i,j) * (M - i - j + 2) * lambda / mu;
% %         Tc1 = Tc1 + p(i,j) * (j - 1);
% %     end
% % end
% % Bpp1 = (To1 - Tc1) / To1; 
% % Tb = 0; To = M*lambda/(lambda + mu);
% % for i = 1:(M-K)
% %     Tb = p(i,K+1)*(M - K - i + 1)*lambda/mu;
% % end
% % Bp = Tb/To
% To = M*lambda/(lambda + mu); Tc = 0; Tb = 0;
% for j = 1:(K + 1)
%     for i = 1:(M - K + 1)
%         %To = To + p(i,j) * (M - i - j + 2) * lambda / mu;
%         Tc = Tc + p(i,j) * (j - 1);
%         Tb = Tb + p(i,j) * (i - 1);
%     end
% end
% % To
% Bpp = (To - Tc) / To 
% Bpp1 = Tb / (Tb + Tc)
end
t2 = cputime;
time = (t2 - t1)/LoopCount
%keyboard
