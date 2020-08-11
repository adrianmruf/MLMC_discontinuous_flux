function N=nsamples(dx0,L)

% Calculates the number of samples at each level for mlmc, see Lemma 5.6 in
% the paper.
% N is a vector with the sample numbers, starting with the sample number for the finest level
% 

s = 1/2; % convergence rate, proved is 1/2, experimental rates are generally higher
w = 2; % w=2 is the exponent for the work estimate of the scheme
% calculate M0
denum=dx0^(2*s)*2^(-2*s*L);

j=1:L;
jsum=2.^(j*(w-s)/2);
Jsum=dx0^(s/2)*sum(jsum);
M0=(1+Jsum)/denum;

Ml=M0*(dx0^(s/2).*2.^(-j*(s+w)/2));

M0=ceil(M0);
Ml=ceil(Ml);
N=[fliplr(Ml) M0];