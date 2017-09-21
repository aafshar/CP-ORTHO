function [ alpha ] = ChooseStep( A,Z,Znormsqr,f,w,G,alpha,R,Lambda,norm_cons )
    %UNTITLED6 Summary of this function goes here
    %   Detailed explanation goes here
    x=tt_fac_to_vec(A);
    gradientf=tt_fac_to_vec(G);
    %etha=1e-8; %bird
    %etha=1e-9; uber
    %beta=0.1; %bird
    etha=1e-8;
    beta=0.5;
    xplus=max(0,x-alpha*gradientf);
    diff_x=xplus-x; %xplus-x
    fplus=CalculateF(tt_cp_vec_to_fac(xplus,Z),Z,Znormsqr,w,R,Lambda,norm_cons);
    diff_f=fplus-f;
    %alpha
    t=0;
    if(diff_f>etha*(gradientf'*diff_x))
        while (diff_f>etha*(gradientf'*diff_x)& t<16 & ~isequal(diff_x,zeros(size(diff_x,1),1)) )
            alpha=alpha*beta;
            xplus=max(0,x-alpha*gradientf);
            [ fplus ]=CalculateF(tt_cp_vec_to_fac(xplus,Z),Z,Znormsqr,w,R,Lambda,norm_cons);
            diff_f=fplus-f;
            diff_x=xplus-x;
            %diff_f
            %etha*(gradientf'*diff_x)
            t=t+1;
        end
        
        
    end
    
    %alpha
end

