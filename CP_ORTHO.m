function [P,fit] = CP_ORTHO(X, R)
    % INPUT
    % X: Observed N-way tensor
    % R: rank of CP factorization
    %
    %
    % OUTPUT
    % P: includes structures of the data
    % fit
    
    PSIa = 0.1; % weight for orthogonality
    PSIb=1; %coefficient unit component
    %Initialization
    sz = size(X); %Tensor Size
    N = length(sz); %number of dimentions
    normX = norm(X); %norm of tensor
    
    S0 = cell(N,1);
    for n=1:N
        S0{n} = rand(sz(n),R);
        for j=1:R
            S0{n}(:,j) = S0{n}(:,j) / norm(S0{n}(:,j));
        end
    end
    Lambda=ones(1,R)';
    %% Fit CP using CPOPT
    normsqr = norm(X)^2;
    %out = feval(fhandle, @(x)Calculate_gradient(x,X,normsqr), tt_fac_to_vec(S0), options);
    
    %Gradient Descent
    %alpha = 0.000001;
    
    fit=0;
    Fit=[];
    %iter=0;
    
    
    while 1
        fitold = fit;
        
        P = ktensor(Lambda,S0{1},S0{2},S0{3});
        normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
        fit = 1 - (normresidual / normX) %fraction explained by model
        Fit=[Fit fit];
        fitchange = abs(fit - fitold);
        % Check for convergence
        if (fitchange < 1e-6)
            break;
        end
        [G_Lambda] = Calculate_Lambda_Gradient(S0,X,Lambda);
        Lambda_step_size=25;
        Lambda_step_size=ChooseLambdaStep(Lambda,S0,X,normsqr,normresidual,G_Lambda,Lambda_step_size);
        Lambda=max(0,Lambda-Lambda_step_size.*G_Lambda{1});
        [f,G] = Calculate_gradient(S0,X,normsqr,Lambda,PSIa,PSIb);
        G1=G(1:3);
        alpha=1;
        alpha=ChooseStep(S0,X,normsqr,f,PSIa,G1,alpha,R,Lambda,PSIb);
        for n= 1:N
            S0{n}=max(0,S0{n}-alpha.*G{n});
        end
    end
    
    % compute factors and model fit
    P = ktensor(Lambda,S0{1},S0{2},S0{3});
    
    normresidual = sqrt( normX^2 + norm(P)^2 - 2 * innerprod(X,P) );
    fit = 1 - (normresidual / normX); %fraction explained by model
    plot([1:size(Fit,2)],Fit,'*r');
    xlabel('Iterations');
    ylabel('Fit');
    
    %% Clean up final result
    % Arrange the final tensor so that the columns are normalized.
    P = arrange(P);
    % Fix the signs
    P = fixsigns(P);
