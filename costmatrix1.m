function [C] = costmatrix1(G_1,G_2,t,K)

%This function constructs the cost matrix of two graphs G_1 and G_2, with 
%weighted adjacency matrices A and B, via the diffusion distance norm.
%Because the diffusion distance is employed, a time t is also required as
%an argument. The parameter K is a natural number to assist in the
%construction of the cost matrix decomposition.

%We assume that G_1 and G_2 are both of order n.

%Decompose adjacency matrices into C_1 and C_2 via the Laplacian
A=adjacency(G_1);
B=adjacency(G_2);
D_1=diag(sum(A,2));
D_2=diag(sum(B,2));
L_1=D_1-A;
L_2=D_2-B;

sz=size(L_1);
n=sz(1,1);

C_1=zeros(n,n);
C_2=zeros(n,n);
for i=1:K
    C_1=C_1+(1/factorial(i))*(-t*L_1)^i;
    C_2=C_2+(1/factorial(i))*(-t*L_2)^i;
end

C_1=C_1/norm(C_1,'fro');
C_2=C_2/norm(C_2,'fro');
E=[C_1 zeros(n,n); zeros(n,n) C_2];

%Generate cost matrix using diffusion distance
C=zeros(n,n);
sigma=1;
for i=1:n
    for j=1:n
        C(i,j)=diffusion_distance(E,C_1(i,:),C_2(j,:),sigma,t);
    end
end
