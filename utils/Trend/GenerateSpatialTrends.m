
% This came out of the fmrilm.m 
% I have never heard of the use of spatial trend + I can't work out what
% the point is. So, I'll pass it for time being. 

n_spatial = 2;
T = 900; 

Y = Y - mean(Y); 
spatial_av = mean(Y)';

if n_spatial>=1 
   trend=spatial_av-mean(spatial_av);
   spatial_trend=(trend*ones(1,n_spatial)).^(ones(T,1)*(1:n_spatial));
else
   spatial_trend=[];
end 

size(spatial_trend)