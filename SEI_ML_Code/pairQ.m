
% r = rootsFQ(f0, Q)
%
% f0 and Q are the resonance frequency and quaity factor
% r is the corresponding pair of roots (i.e., zeros or poles)
%
% if Q == 0, a real root is assumed (i.e., r = 2 * pi * f0)

function r = pairQ(f0, Q)

  r0 = 2 * pi * f0;
  nc = find(Q ~= 0);

  res = sqrt(1 - 4 * Q(nc).^2);
  mag = r0(nc) ./ (2 * Q(nc));
  r0(nc) = mag .* (1 + res);
  r = [r0, mag .* (1 - res)]/(2*pi);
end
 