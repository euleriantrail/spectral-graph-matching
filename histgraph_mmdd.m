function d=histgraph_mmdd(G,x,y,t)

%Let t be nat. number representing a given time. G is a given graph, x and y
%are nodes in G, d returns the diffusion distance
%between x and y. Based on paper by Maggioni and Murphy.

%Generate symmetric matrix (P_sym) that is conjugate to the transition
%matrix. Generate Gaussian kernel matrix W with parameter sigma=1.
W=adjacency(G);
sz=size(W);
n=sz(1);
D=zeros(n); 
for i=1:n
    D(i,i)=diag(sum(W,2));
end
P_sym = (D^(-1/2))*W*(D^(-1/2));

%Construct Laplacian
L=D-W;

%Generate sorted eigenvalues and eigenvectors of P_sym
[phi,lambda]=eig(P_sym);
[L,ind]=sort(diag(lambda));
lambda_s=lambda(ind,ind);
phi_s=phi(:,ind);

%Normalize eigenvectors of P_sym with pi(x_i) scalars into psi matrix
pi=zeros(n);
for i=1:n
    pi(i)=(D(i,i))/n;
end

psi=pi.*phi_s;

%Aggregate info from columns of psi matrix with histograms
for i=1:n
    h=histogram(psi(:,i));
    n=h.Values;
    psi(:,i)=n';
end

sig=zeros(1,n);
for i=1:n
    sig(i)=(lambda_s(i,i)^(2*t))*(psi(:,i)*x-psi(:,i)*y)^2;
end
sumsig=sum(sig);
d=sqrt(sumsig);