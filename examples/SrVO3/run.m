datadir = 'data';
seedname = 'SrVO3';
wann_bands = 1:3;

addpath('../../');

U(iband,jband)=main(datadir,seedname,wann_bands,1,1);
U(iband,jband)=main(datadir,seedname,wann_bands,1,2);
U(iband,jband)=main(datadir,seedname,wann_bands,1,3);

U(iband,jband)=main(datadir,seedname,wann_bands,2,2);
U(iband,jband)=main(datadir,seedname,wann_bands,2,3);

U(iband,jband)=main(datadir,seedname,wann_bands,3,3);