function [ G ] = Calculate_Lambda_Gradient( A,Z,Lambda )
    %UNTITLED10 Summary of this function goes here
    %   
    % %Part of the gradient for Lambda
    N = ndims(Z);
    R = size(A{1},2);
    G = cell(1,1);
    C =cell(R,N);
    Upsilon = cell(N,1);
    for n = 1:N
        Upsilon{n} = A{n}'*A{n};
    end
    Gamma = cell(1,1);
    Gamma{1} = ones(R,R);
    for m = [2:N]
        Gamma{1} = Gamma{1} .* Upsilon{m};
    end
    Teta =  Gamma{1} .* Upsilon{1};
    for n=1:N
        for r=1:R
            C{r,n} = A{n}(:,r);
        end
    end
    G{1} = zeros(R,1);
    for r=1:R
        G{1}(r) = G{1}(r) -ttv(Z, C(r,:),1:N);
    end
    G{1} = G{1} + sum(khatrirao(Teta,Lambda'),2);
    
end

