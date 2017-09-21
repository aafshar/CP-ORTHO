function [ X ] = add_noise( X,percent )
    %add noise to the data.
    %For instaance, 0.1 means add noise to 10% of the samples.
    X=tensor(X);
    SIZ=size(X);
    n=SIZ(1,1);
    m=SIZ(1,2);
    k=SIZ(1,3);
    for i=1:n
        for j=1:m
            for k=1:k
                random_number=rand();
                if(random_number<percent)
                    X(i,j,k)=max(0,X(i,j,k)+random('norm', 0,1,1));
                end
            end
        end
        
    end
    
