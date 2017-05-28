datadir = '~/Downloads/SrVO3_4x4x4/';
seedname = 'SrVO3';
wann_bands = 1:3;

addpath('../../');
U_11 = main(datadir,seedname,wann_bands,1,1);
U_22 = main(datadir,seedname,wann_bands,2,2);
U_33 = main(datadir,seedname,wann_bands,3,3);