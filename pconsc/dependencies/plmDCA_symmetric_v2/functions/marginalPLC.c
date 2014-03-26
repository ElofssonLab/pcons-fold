#include <math.h>
#include "mex.h"

/*
This is a reworked version of the pseudolikelihood objective contained in Mark Schmidt's thesis code (http://www.di.ens.fr/~mschmidt/Software/thesis.html). The copyright conditions for the original code is included below.
---------------------------------

Copyright 2005-2012 Mark Schmidt. All rights reserved.

Redistribution and use in source and binary forms, with or without modification, are
permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice, this list of
      conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice, this list
      of conditions and the following disclaimer in the documentation and/or other materials
      provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE AUTHOR ''AS IS'' AND ANY EXPRESS OR IMPLIED
WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND
FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE AUTHOR OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    /* variables */
    char param;
    int i, s, t, n, e, n1, n2, nInstances, nNodes, nStates, nEdges,
            *edges, *y, y1, y2, *rint,r;
    double *yr, *g1r, *g2r, *logPot, *z, *fval_r, *w1r, *w2r, *nodeBel, *lambdas;
    
    /* input */
    param = *(mxChar*)mxGetData(prhs[0]);
    y = (int*)mxGetPr(prhs[1]);
    yr = mxGetPr(prhs[2]);
    edges = (int*)mxGetPr(prhs[3]);
    w1r = mxGetPr(prhs[4]);
    w2r = mxGetPr(prhs[5]);
    lambdas = mxGetPr(prhs[6]);
    rint = (int*)mxGetPr(prhs[7]);
    
    /* compute sizes */
    nInstances = mxGetDimensions(prhs[1])[0];
    nNodes = mxGetDimensions(prhs[1])[1];
    nStates = mxGetDimensions(prhs[4])[1];
    nEdges = mxGetDimensions(prhs[3])[0];
    
    /* allocate memory */
    logPot = mxCalloc(nStates, sizeof(double));
    z = mxCalloc(1, sizeof(double));
    nodeBel = mxCalloc(nStates, sizeof(double));
    
    /* output */
    plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
    fval_r = mxGetPr(plhs[0]);
    *fval_r = 0;
    plhs[1] = mxCreateDoubleMatrix(nStates, 1, mxREAL);
    g1r = mxGetPr(plhs[1]);
    plhs[2] = mxCreateDoubleMatrix(nStates*nStates*(nNodes-1), 1, mxREAL);
    g2r = mxGetPr(plhs[2]);
    
    r=rint[0]-1;

    for(i=0;i < nInstances;i++) {
        for(s=0;s < nStates;s++) {
            logPot[s] = w1r[s];
        }

        for(n2 = 0;n2 < nNodes;n2++) {
            y2 = y[i + nInstances*n2];       
            if(n2!=r) {
                for(s=0; s<nStates; s++) {
                    logPot[s] += w2r[s+nStates*(y2+nStates*(n2-(n2>r)))];               
                }
            }
            
        }
      
        z[0] = 0;
        for(s = 0; s < nStates; s++) {
            z[0] += exp(logPot[s]);
        }
        *fval_r -= yr[i]*logPot[y[i+nInstances*r]];
        *fval_r += yr[i]*log(z[0]);
 
        
        
	/*Gradient:*/
        
 
        for(s = 0; s < nStates; s++) {
            nodeBel[s] = exp(logPot[s] - log(z[0]));
        }
                       
        y1 = y[i + nInstances*r]; 
        g1r[y1] -= yr[i]*1;
            
        for(s=0; s < nStates; s++) {
            g1r[s] += yr[i]*nodeBel[s];
        }
   
        for(n2=0;n2<nNodes;n2++) {
            if(n2!=r) {
                y2 = y[i + nInstances*n2];

                g2r[y1+nStates*(y2+nStates*(n2-(n2>r)))] -= yr[i];                

                for(s=0;s<nStates;s++) {
                   g2r[s+nStates*(y2+nStates*(n2-(n2>r)))] += yr[i]*nodeBel[s];
                }	
            }
        }
        
    }

  
    
    
    
    
    /*Add contributions from R_l2*/

    for(s = 0; s < nStates; s++) {
        *fval_r += lambdas[0]*w1r[s]*w1r[s];
        g1r[s] += lambdas[0]*2*w1r[s]; 
    }

    for(n2 = 0;n2 < nNodes;n2++) {
        if(n2!=r) {
            for(s = 0; s < nStates; s++) {
                for(t = 0; t < nStates; t++) {                   
                        *fval_r += lambdas[1]*w2r[s+nStates*(t+nStates*(n2-(n2>r)))]*w2r[s+nStates*(t+nStates*(n2-(n2>r)))];
                        g2r[s+nStates*(t+nStates*(n2-(n2>r)))] += lambdas[1]*2*w2r[s+nStates*(t+nStates*(n2-(n2>r)))];                    
                }
            }
        }
    }
             
    mxFree(logPot);
    mxFree(z);
    mxFree(nodeBel);
    return;
}
