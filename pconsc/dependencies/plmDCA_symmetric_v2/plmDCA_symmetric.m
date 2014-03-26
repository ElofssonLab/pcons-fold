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




function plmDCA_symmetric(fastafile,outputfile,lambda_h,lambda_J,reweighting_threshold,nr_of_cores)
    options.method='cg';	%Minimization scheme. Default: 'lbfgs', 'cg' for conjugate gradient (use 'cg' if out of RAM)
    options.optTol=1e-5; 	%default: 1e-5
    options.progTol=1e-7;   	%default: 1e-9   
    addpath(genpath(pwd))
    
%Read inputfile (removing inserts), remove duplicate sequences, and calculate weights (Yr) and B_eff  
    [N,B_with_id_seq,q,Y]=return_alignment(fastafile);
    Y=unique(Y,'rows');
    [B,N]=size(Y);
    Yr = ones(B,1);
    if reweighting_threshold>0.0
        fprintf('Starting to calculate weights \n...');
        tic
        %Reweighting in MATLAB:            
        %Yr = (1./(1+sum(squareform(pdist(Y,'hamm')<=reweighting_threshold))))';       
		     
        %Reweighting in C:
        m=zeros(B,1);
        Y=int32(Y);
        calc_inverse_weights(Y-1,m,reweighting_threshold);
        Yr=1./m;
        fprintf('Finished calculating weights \n');
        toc
    end
    B_eff=sum(Yr);
    fprintf('### N = %d B_with_id_seq = %d B = %d B_eff = %.2f q = %d\n',N,B_with_id_seq,B,B_eff,q);

%Set up and run optimizer    
    field_lambda=lambda_h*B_eff;    
    coupling_lambda=lambda_J*B_eff;
    edges=[];
    for i=1:N-1
        for j=i+1:N
            edges=[edges;[i,j]];
        end
    end
    Y=int32(Y);q=int32(q);edges=int32(edges);
    funObj=@(w)pseudo_likelihood_symmetric(w,Y,Yr,q,edges,field_lambda,coupling_lambda);    
    w0=zeros(q*N+q^2*N*(N-1)/2,1);
    if nr_of_cores>1
        matlabpool('open',nr_of_cores) 
        tic
        w=minFunc(funObj,w0,options);    
        toc
        matlabpool('close') 
    else
        tic
        w=minFunc(funObj,w0,options);    
        toc
    end
        
    Jtemp=reshape(w(q*N+1:end),q,q,N*(N-1)/2);

       
%Change parameter gauge for J                      
    J=zeros(q,q,N*(N-1)/2);
    for l=1:(N*(N-1)/2)
        J(:,:,l)=Jtemp(:,:,l)-repmat(mean(Jtemp(:,:,l)),q,1)-repmat(mean(Jtemp(:,:,l),2),1,q)+mean(mean(Jtemp(:,:,l)));
    end
      

%Calculate frob. norms FN_ij
    NORMS=zeros(N,N); 
    l=1;
    for i=1:(N-1)
        for j=(i+1):N
            NORMS(i,j)=norm(J(1:end,1:end,l),'fro');
            NORMS(j,i)=NORMS(i,j);
            l=l+1;
        end
    end               
       
    
%Calculate final scores, CN_ij=FN_ij-(FN_i-)(FN_-j)/(FN_--), where '-'
%denotes average
    norm_means=mean(NORMS)*N/(N-1);
    norm_means_all=mean(mean(NORMS))*N/(N-1);
    CORRNORMS=NORMS-norm_means'*norm_means/norm_means_all;
    output=[];
    for i=1:(N-1)
        for j=(i+1):N
            output=[output;[i,j,CORRNORMS(i,j)]];
        end
    end
    dlmwrite(outputfile,output,'precision',5)
end





















function [N,B,q,Y] = return_alignment(inputfile)
% reads alignment from inputfile, removes inserts and converts into numbers
    align_full = fastaread(inputfile);
    B = length(align_full);
    ind = align_full(1).Sequence ~= '.' & align_full(1).Sequence == upper( align_full(1).Sequence );
    N = sum(ind);
    Y = zeros(B,N);

    for i=1:B
        counter = 0;
        fprintf('#!!# N = %d\n',i);
        for j=1:length(ind)
            if( ind(j) )
                counter = counter + 1;
                Y(i,counter)=letter2number( align_full(i).Sequence(j) );
            end
        end
    end
    q=max(max(Y));
end

function x=letter2number(a)
    switch(a)
        % full AA alphabet
        case '-'
             x=1;
        case 'A'    
            x=2;    
        case 'C'    
            x=3;
        case 'D'
            x=4;
        case 'E'  
            x=5;
        case 'F'
            x=6;
        case 'G'  
            x=7;
        case 'H'
            x=8;
        case 'I'  
            x=9;
        case 'K'
            x=10;
        case 'L'  
            x=11;
        case 'M'
            x=12;
        case 'N'  
            x=13;
        case 'P'
            x=14;
        case 'Q'
            x=15;
        case 'R'
            x=16;
        case 'S'  
            x=17;
        case 'T'
            x=18;
        case 'V'
            x=19;
        case 'W'
            x=20;
        case 'Y'
            x=21;
        otherwise
            x=1;
    end
end












