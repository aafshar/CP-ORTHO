function [f,G] = Calculate_gradient(A,Z,Znormsqr,Lambda,w,norm_cons)
%TT_CP_FUN Calculate function and gradient for CP fit function.
%
%  [F,G] = TT_CP_FUN(X,Z) where X is a vector containing the entries of the
%  components of the model and Z is the tensor to be fit.
%
%  See also TT_CP_VEC_TO_FAC, TT_FAC_TO_VEC, TT_CP_FG, CP_OPT
%
%MATLAB Tensor Toolbox.
%Copyright 2015, Sandia Corporation.

% This is the MATLAB Tensor Toolbox by T. Kolda, B. Bader, and others.
% http://www.sandia.gov/~tgkolda/TensorToolbox.
% Copyright (2015) Sandia Corporation. Under the terms of Contract
% DE-AC04-94AL85000, there is a non-exclusive license for use of this
% work by or on behalf of the U.S. Government. Export of this data may
% require a license from the United States Government.
% The full license terms can be found in the file LICENSE.txt


%% Convert x to a cell array of matrices
%A = tt_cp_vec_to_fac(x,Z);

%% Call cp_fit and cp_gradient using cp_fg


N = ndims(Z);

if ~iscell(A) && ~isa(A,'ktensor');
    error('A must be a cell array or ktensor');
end

if isa(A,'ktensor')
    A = tocell(A);
end
R = size(A{1},2);

%% Upsilon and Gamma
Upsilon = cell(N,1);
for n = 1:N
    Upsilon{n} = A{n}'*A{n};
end

Gamma = cell(N,1);
for n = 1:N
    Gamma{n} = ones(R,R);
    for m = [1:n-1,n+1:N]
        Gamma{n} = Gamma{n} .* Upsilon{m};
    end
end

LL   =  Lambda*Lambda';
Teta =  Gamma{1} .* Upsilon{1};   

%% Calculation
%A{2}'*A{2}
%F1
if exist('Znormsqr','var')
    f_1 = Znormsqr;
else
    f_1 = norm(Z)^2;
end

%% Calculate gradient and F2 of CP part
%norm_cons*(A{1}-normc(A{1}))
G = cell(N+1,1); %gradient of two terms of objective function
G1 = cell(N+1,1);
U = khatrirao(mttkrp(Z,A,1), Lambda');
V = A{1} .* U;
f_2 = sum(V(:));
%norm_cons*(A{1}-normc(A{1}))
G1{1} = -U + A{1}*(LL.*Gamma{1})+norm_cons*(A{1}-normc(A{1}));
for n = 2:N
    U = khatrirao(mttkrp(Z,A,n),Lambda');
    %norm_cons*(A{n}-normc(A{n}))
    G1{n} = -U + A{n}*(LL.*Gamma{n})+norm_cons*(A{n}-normc(A{n}));
end

R = size(G1{1},2);
% % add norm constraints to gradient
% for i=1:N
%     for r=1:R         
%          G1{i}(:,r) = G1{i}(:,r) + norm_cons*(A{i}(:,r)-(A{i}(:,r)/norm(A{i}(:,r))));
%     end
% end



%F3
W   = Teta.*LL;
f_3 = sum(W(:));

%SUM of pieces of the loss function
f = 0.5 * f_1 - f_2 + 0.5 * f_3;

% add norm constraints to f
for i=1:N
    for r=1:R        
        f = f + 0.5* norm_cons * (norm(A{i}(:,r))-1)^2;
    end
end


% %Part of the gradient for Lambda
% C =cell(R,N);
% for n=1:N
%     for r=1:R
%         C{r,n} = A{n}(:,r);
%     end
% end
% G1{N+1} = zeros(R,1);
% for r=1:R
%      G1{N+1}(r) = G1{N+1}(r) -ttv(Z, C(r,:),1:N);
% end
% G1{N+1} = G1{N+1} + sum(khatrirao(Teta,Lambda'),2);

%% Convert a cell array to a vector
%g1 = tt_fac_to_vec(G1);

%Calculate gradient of ORTHOGONAL part
%4*((Upsilon{n}(r,r)*Gamma{n}(r,r)-1)*Gamma{n}(r,r))*A{n}(:,r);
g2=[];% vecotr of gradients
G2 = cell(N,1);
for n = 1:N
  for r=1:R  
   M=4*Upsilon{n}(r,1:R).*((Gamma{n}(r,1:R)).^2);
   M=repmat(M,size(A{n},1),1).*A{n};
   Orthog=sum(M')'-M(:,r)+4*((Upsilon{n}(r,r)*Gamma{n}(r,r)-1)*Gamma{n}(r,r))*A{n}(:,r); %subtract from the rth vector.
   G2{n}(:,r)=Orthog;
   g2=[g2; Orthog];
  end
  G{n}=G1{n}+(w*G2{n});
end
G{N+1}=G1{N+1};
fortho=Teta-eye(R);
f=f+(w/2)*sum(fortho(:).^2);
%g=g1+g2;