% Copyright 2012 - by Magnus Ekeberg (magnus.ekeberg@gmail.com)
% All rights reserved
% 
% Permission is granted for anyone to copy, use, or modify this
% software for any uncommercial purposes, provided this copyright 
% notice is retained, and note is made of any changes that have 
% been made. This software is distributed without any warranty, 
% express or implied. In no event shall the author or contributors be 
% liable for any damage arising out of the use of this software.
% 
% The publication of research using this software, modified or not, must include an 
% appropriate citation to:
%	M. Ekeberg, C. Lövkvist, Y. Lan, M. Weigt, E. Aurell, Improved contact prediction
%	in proteins: Using pseudolikelihoods to infer Potts models, Phys. Rev. E 87, 012707 (2013)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




function [pseudoNLL,g] = pseudo_likelihood_symmetric(w,Y,Yr,q,edges,lambdah,lambdaJ)
[B,N] = size(Y);
w1=reshape(w(1:q*N),N,q); w2=reshape(w(q*N+1:end),q,q,N*(N-1)/2);

f=zeros(1,N); grad=zeros(q+q^2*(N-1),N); J_r=zeros(q^2*(N-1),N);
for r=1:N
   ls1=edges(:,2)==r;ls2=edges(:,1)==r;
   J_r(1:q^2*(r-1),r)=reshape(permute(w2(:,:,ls1),[2 1 3]),q^2*(r-1),1);
   J_r(q^2*(r-1)+1:end,r)=reshape(w2(:,:,ls2),q^2*(N-r),1);
end
parfor r=1:N
    [f(r) grad(:,r)]=marginalPL([w1(r,:)';J_r(:,r)],Y,Yr,q,edges,lambdah,lambdaJ/2,r);
end

g1=grad(1:q,:)'; grad_reshaped=reshape(grad(q+1:end,:),q,q,N-1,N); g2=zeros(q,q,N*(N-1)/2);
l=1;
for i=1:N-1
    for j=i+1:N
        g2(:,:,l)=grad_reshaped(:,:,j-1,i)+grad_reshaped(:,:,i,j)';
        l=l+1;
    end
end

pseudoNLL=sum(f); g=[g1(:);g2(:)];
