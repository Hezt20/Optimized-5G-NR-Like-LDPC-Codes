% Parity Check Matrices for Optimized 5G-NR-Like LDPC codes
% with \Phi_A = {13,26} and K=6, M<=12 in the manuscript 
% "Multi-Layer Progressive Tree-Structured Edge Growth Algorithm
% for Nested Lifting Design of 5G-NR-Like LDPC Codes"

clear; clc; close all;

% Exponent matrix with -1 entries correspond to zero entries in the base matrix
% and other entries with value in [0,Z-1] correspond to lifting values for the edges 
% in the same position of the base matrix.
A = [
    0	0	0	0	-1	-1	0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1
    23	-1	-1	11	0	0	-1	0	0	-1	-1	-1	-1	-1	-1	-1	-1	-1
    3	21	-1	18	15	-1	9	-1	0	0	-1	-1	-1	-1	-1	-1	-1	-1
    -1	20	1	-1	19	3	12	-1	-1	0	-1	-1	-1	-1	-1	-1	-1	-1
    5	7	-1	-1	-1	-1	-1	0	-1	-1	0	-1	-1	-1	-1	-1	-1	-1
    4	15	-1	-1	-1	3	-1	0	-1	-1	-1	0	-1	-1	-1	-1	-1	-1
    12	-1	-1	-1	-1	5	-1	0	-1	-1	-1	-1	0	-1	-1	-1	-1	-1
    -1	2	-1	-1	-1	21	-1	23	-1	0	-1	-1	-1	0	-1	-1	-1	-1
    19	0	-1	-1	-1	-1	-1	-1	0	-1	-1	-1	-1	-1	0	-1	-1	-1
    -1	16	-1	-1	-1	-1	20	0	-1	-1	-1	-1	-1	-1	-1	0	-1	-1
    10	0	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	-1	0	-1
    24	-1	-1	-1	-1	-1	-1	-1	-1	0	-1	-1	-1	-1	-1	-1	-1	0
];

% Parameters for the PCM, where 4<=M<=12 and Z belongs to {13,26}.
K=6; M=12; N=K+M; V=2; Rate = K/(K+M-V); Z=13;

% Construct PCM from the exponent matrix.
A = A(1:M,1:K+M);
A = mod((A>=0).*A,Z)+(A<0).*A;
nzmax = sum(sum(A~=-1))*Z; 
irow(1:nzmax) = 0; jcol(1:nzmax) = 0; 

k=1:Z;
for i1=1:M
    for j1=1:N
        if A(i1,j1)~=-1
            irow(k) = i1*Z-Z+1 + (0:Z-1);
            jcol(k) = j1*Z-Z+1 + mod(A(i1,j1)+(0:Z-1),Z);
            k = k+Z;
        end
    end
end

PCM = logical( sparse(irow',jcol',ones(1,nzmax),M*Z,N*Z,nzmax) );
