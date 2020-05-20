function Usan = establishUsanThresh(epivol)
 
  num = size(epivol,1); % T?
% film_gls.cc line 159:  
%   // Residuals container:
%   Matrix residuals(sizeTS, numTS);
%
%   // Setup autocorrelation estimator:
%   AutoCorrEstimator acEst(residuals);
  
  modemode = mode(mode(epivol)); 
  %int num = epivol.Nrows();
  %Histogram hist(epivol, max(num/200,1));
  %hist.generate();
  %float mode = hist.mode();
  %cout<< "mode = " << mode << endl;
  sumsum = 0;
  % Work out standard deviation from mode for values greater than mode:
  count = 0;
  for i = 1:num 
    if(epivol(i) > modemode)
      sumsum = sumsum + ((epivol(i) - modemode)*(epivol(i) - modemode));
      count = count + 1;
    end
  end

  sig = sqrt(sumsum/num);
  %cout<< "sig = " << static_cast<int>(sig) << endl;
  Usan = sig/3;
end