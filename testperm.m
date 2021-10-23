function P_hat=testperm(C)

sz=size(C);
n=sz(1);

A_eq=zeros(2*n,n*n);

for i=1:n
    A_eq(i,(i-1)*n+1:i*n)=1;
end

for i=1:n
    A_eq(i+n,i:n:i+(n-1)*n)=1;
end

P_hat=linprog(reshape(C,[n^2,1]),[],[],A_eq,ones(2*n,1),zeros(n*n,1),[]);
