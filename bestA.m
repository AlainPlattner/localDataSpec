function A=bestA(S,Sp)
% A=bestA(S,Sp)
% 
% Find the best factor to multiply spec S such that it fits spec Sp
%
% Last modified by plattner-at-alumni.ethz.ch, 01/17/2020
  
% Make sure all specs and l are column
S=S(:); Sp=Sp(:);

A=exp(mean(log(Sp)-log(S)));






