function [Trends,numSpline] = GenerateTemporalTrends(T,TR,numTremd)
% Generates trends using splines. 
% bits and pieces from fmrilm.m of Keith Worsely. 
% 
% SA, Ox, 2020

   if nargin<3; numTremd = 3; end

    Y = 1:T;
    numSpline    = round(numTremd*TR*T/360);   
    trend        = ((2*Y-(max(Y)+min(Y)))./(max(Y)-min(Y)))';
    if numSpline<=3
      Trends    = (trend*ones(1,numSpline+1)).^(ones(T,1)*(0:numSpline));
    else
      Trends    = (trend*ones(1,4)).^(ones(T,1)*(0:3));
      knot      = (1:(numSpline-3))/(numSpline-2)*(max(Y)-min(Y))+min(Y);
      for k=1:length(knot)
         cut    = Y'-knot(k);
         Trends = [Trends (cut>0).*(cut./max(cut)).^3];
      end
    end
    % temporal_trend=temporal_trend*inv(chol(temporal_trend'*temporal_trend));
end