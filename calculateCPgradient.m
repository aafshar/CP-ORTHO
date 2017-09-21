function [ f g ] = calculateCPgradient( X,S0,i,x,N,R )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
% Calculate objective f
%  Xn = permute(X,[i 1:i-1,i+1:N]);
%  Xn = reshape(Xn.data, size(X,i), []);
%% Upsilon and Gamma
Upsilon = cell(N,1);
for n = 1:N
    if (n==i)
        Upsilon{n} = x'*x;
    else     
        Upsilon{n} = S0{n}'*S0{n};
    end
end

Gamma = cell(N,1);
for n = 1:N
    Gamma{n} = ones(R,R);
    for m = [1:n-1,n+1:N]
        Gamma{n} = Gamma{n} .* Upsilon{m};
    end
end


% if (i==1)
%     piT=khatrirao(S0{3},S0{2},'r')';
%     %C=(S0{2}'*S0{2}).*(S0{3}'*S0{3});
% elseif(i==2)
%     piT=khatrirao(S0{3},S0{1},'r')';
%     %C=(S0{1}'*S0{1}).*(S0{3}'*S0{3});
% else
%     piT=khatrirao(S0{2},S0{1},'r')';
%     %C=(S0{2}'*S0{2}).*(S0{1}'*S0{1});
% end
% 
% f1=(0.5*norm(Xn-x*piT,'fro')^2);


f_1 = norm(X)^2;
W = Gamma{1} .* Upsilon{1};
f_3 = sum(W(:));
U = mttkrp(X,S0,i);
V = x .* U;
f_2 = sum(V(:));

%SUM
f = 0.5 * f_1 - f_2 + 0.5 * f_3;


if nargout > 1 % gradient required
    U = mttkrp(X,S0,i);
    g = -U + S0{i}*Gamma{i};
end

end

