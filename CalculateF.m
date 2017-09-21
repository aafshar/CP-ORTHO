function [ f ] = CalculateF( A,Z,Znormsqr,w,R,Lambda,norm_cons )
    %UNTITLED7 Summary of this function goes here
    %   Detailed explanation goes here
    N = ndims(Z);
    Upsilon = cell(N,1);
    for n = 1:N
        Upsilon{n} = A{n}'*A{n};
    end
    Gamma = cell(1,1);
    Gamma{1} = ones(R,R);
    for m = [2:N]
        Gamma{1} = Gamma{1} .* Upsilon{m};
    end
    LL   =  Lambda*Lambda';
    Teta =  Gamma{1} .* Upsilon{1};
    f_1 = Znormsqr;
    U = khatrirao(mttkrp(Z,A,1), Lambda');
    V = A{1} .* U;
    f_2 = sum(V(:));
    W   = Teta.*LL;
    f_3 = sum(W(:));
    f = 0.5 * f_1 - f_2 + 0.5 * f_3;
    fortho=Teta-eye(R);
    % add norm constraints to f
    for i=1:N
        for r=1:R
            f = f + 0.5* norm_cons * (norm(A{i}(:,r))-1)^2;
        end
    end
    f=f+(w/2)*sum(fortho(:).^2);
    
end

