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




function [fval_r,gr] = marginalPL(wr,Y,Yr,q,edges,lambdah,halflambdaJ,r)
[B,N] = size(Y);
w1r=reshape(wr(1:q),1,q);
w2r=reshape(wr(q+1:end),q,q,N-1);

r=int32(r);
[fval_r,g1r,g2r] = marginalPLC('F',Y-1,Yr,edges-1,w1r,w2r,[lambdah;halflambdaJ],r);
gr = [g1r(:);g2r(:)];
