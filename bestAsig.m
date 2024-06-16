function A=bestAsig(S,Sp,sig)
  % A=bestAsig(S,Sp,sig)
  %
  % Find the best factor to multiply spec S such that it fits spec Sp,
  % but weighted by the uncertainty of each degree
  % A*S=Sp;
  %
  % This calculates the optimal amplification factor to minimize
  % chi square misfit between Sp and S with standard deviation sig.
  % I have a proof for this but there is not enough space in this help
  % text to write it down (just write down chi-square formula,
  % take derivative w.r.t. amplitude, and set to zero).
  %
  % Last modified by plattner-at-alumni.ethz.ch, 06/16/2024

  S = S(:); Sp = Sp(:); sig = sig(:);

  F = S./sig; G = Sp./sig;

  A = (F'*G)/(F'*F);
