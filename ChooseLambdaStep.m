function [ Lambda_step_size ] = ChooseLambdaStep( Lambda,A,X,normsqr,f,G,Lambda_step_size )
    %UNTITLED8 Summary of this function goes here
    %   Detailed explanation goes here
    gradientLambda=tt_fac_to_vec(G);
    %etha=1e-9; %bird
    %etha=1e-4; uber
    %beta=0.4; %bird
    etha=1e-5;
    beta=0.6;
    Lambdaplus=Lambda-Lambda_step_size*gradientLambda;
    diff_Lambda=Lambdaplus-Lambda; %xplus-x
    P = ktensor(Lambdaplus,A{1},A{2},A{3});
    fplus = sqrt( normsqr + norm(P)^2 - 2 * innerprod(X,P) );
    diff_f=fplus-f;
    %alpha
    t=0;
    if(diff_f>etha*(gradientLambda'*diff_Lambda))
        while (diff_f>etha*(gradientLambda'*diff_Lambda)& t<16 & ~isequal(diff_Lambda,zeros(size(diff_Lambda,1),1)) )
            Lambda_step_size=Lambda_step_size*beta;
            Lambdaplus=Lambda-Lambda_step_size*gradientLambda;
            P = ktensor(Lambdaplus,A{1},A{2},A{3});
            fplus = sqrt( normsqr + norm(P)^2 - 2 * innerprod(X,P) );
            diff_f=fplus-f;
            diff_Lambda=Lambdaplus-Lambda;
            %diff_f
            %etha*(gradientLambda'*diff_Lambda)
            t=t+1;
        end
        
        
    end
    
    %alpha
    
end

