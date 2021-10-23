function d=graph_mmdd(G,i,j,t)

%Let t be nat. number representing a given time. G is a given graph, x and 
%y are nodes in G, d returns the diffusion distance
%between the ith and jth nodes of G. Based on paper by Maggioni and Murphy.

%Generate symmetric matrix (P_sym) that is conjugate to the transition
%matrix. Generate Gaussian kernel matrix W with parameter sigma=1.

W=full(adjacency(G));
sz=size(W);
n=sz(1);
D=diag(sum(W,2));
P_sym = (D^(-1/2))*W*(D^(-1/2));

%Construct Laplacian
L=D-W;

%Generate sorted eigenvalues and eigenvectors of P_sym
[phi,lambda]=eig(P_sym);
[Ls,ind]=sort(diag(lambda));
lambda_s=lambda(ind,ind);
phi_s=phi(:,ind);

%Normalize eigenvectors of P_sym with pi(x_i) scalars into psi matrix
pi=zeros(n,1);
for k=1:n
    pi(k)=(D(k,k))/n;
end

psi=pi.*phi_s;
lambda_bar=lambda_s.^(2*t);

%Sum the following line instead of componentwise multiplication (see Abiy
% email).
psi_bar(i,j)=sum((psi(i,:)-psi(j,:)));

d=sum(sum(lambda_bar.*psi_bar(i,j)));


%sig=zeros(1,n);
%for i=1:n
%    sig(i)=(lambda_s(i,i)^(2*t))*(psi(:,i)*x-psi(:,i)*y)^2;
%end
%sumsig=sum(sig);
%d=sqrt(sumsig);