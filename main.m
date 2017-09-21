clear allh
clc

addpath(genpath('./tensor_toolbox'));

dim = [100,100,100]; % dimension of tensor is 10*10*10
R = 7; % rank for tensor is 5
X = Synthetic_data(dim,R);
%N=add_noise( X,0.8 );


[T, fit]= CP_ORTHO(X, R);
[ortho, A, flag, best_perm] = score(T,X)

% check orthogonality of the tensor on mode 2
(T{1}'*T{1}).*(T{2}'*T{2}).*(T{3}'*T{3})
