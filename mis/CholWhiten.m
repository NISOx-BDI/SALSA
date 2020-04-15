function [Vhalf,spdflag] = CholWhiten(COV)
% [W,spdflag] = CholWhiten(COV)
% Whiten with Cholesky decomposition & inverse
%
  [R,spdflag]   = chol(COV);
  Vhalf         = inv(R');
  %W             = pinv(R); %FUCK!
end