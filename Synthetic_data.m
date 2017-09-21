function X = Synthetic_data(dim,R)
    % Produces an m x n set of orthonormal vectors,
    % (thats n vectors, each of length m)
    %
    % Inputs should be two scalars, m and n, where n is smaller than
    % or equal to m.
    A = cell(3,1);
    A{1}=randi([0,1],dim(1),R);
    A{2}=randi([0,1],dim(2),R);
    A{3}=randi([0,1],dim(3),R);
    
    for i=1:R
        for j=i+1:R
            mode=randi([1,3]);
            [row1,col1] = find(A{mode}(:,i)==1);
             A{mode}(row1,j)=0;
        end
    end
A=normc(A);
(A{1}'*A{1}).*(A{2}'*A{2}).*(A{3}'*A{3});
A = {A{1},A{2},A{3}};
Lambda=50*(rand(R,1)+1);
X = ktensor(Lambda,A);
end

