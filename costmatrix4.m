function C=costmatrix4(G1,G2,t)

%This fourth cost matrix approach uses histograms of diffusion distance
%node behavior in graphs G1,G2. The histograms each have n bins.
%The l2 distance of these histograms is then taken.

%As usual, we assume G1,G2 are both of order n.

sz=size(G1.Nodes);
n=sz(1,1);

C=zeros(n,n);
H1=zeros(n,n);
H2=zeros(n,n);
for i=1:n
    for j=1:n
        H1(i,:)=histogram(node_diffu_dist(G1,i,t),n).Values;
        H2(j,:)=histogram(node_diffu_dist(G2,j,t),n).Values;
    end
end
E=[H1 zeros(n,n); zeros(n,n) H2];
for i=1:n
    for j=1:n
        C(i,j)=norm(E(i,:)-E(j,:));
    end
end
C=C/norm(C,'fro');