
% y = [1120 1160  963 1210 1160 1160  813 1230 1370 1140  995  935 ...
%      1110  994 1020  960 1180  799  958 1140 1100 1210 1150 1250 ...
%      1260 1220 1030 1100  774  840  874  694  940  833  701  916 ...
%      692 1020 1050  969  831 726  456  824  702 1120 1100  832  764 ...
%      821  768  845 864  862 698  845  744  796 1040  759  781  865 ...
%      845  944  984  897  822 1010  771  676  649  846  812  742  801 ...
%      1040  860  874 848  890  744  749  838 1050  918  986  797  923 ...
%      975  815 1020  906  901 1170  912  746  919  718  714  740]';

y = dY(:,40); 

options = struct('order',1,'opt',1,'winds',[0 1]);
dlm = dlmfit(y,122,[0 1.65],[],[],X,options);

y = dlm.resid2; % data are residual from the previous fit



y = residY(:,40); 
phi = [0.3]; % initial guess for AR parameter
opts = struct('order',-1,'seas',0,'arphi',phi);
dlm = dlmfit(y,0.00001,[1],[],[],X,opts); % initial DLM fit
% define object function for minimizing the -2*log(likelihood)
ofun = @(phi)getfield(dlmsmo(dlm.y,dlm.F,dlm.V,dlm.x0,assifun(dlm.G,phi,1:size(dlm.G,1),1),dlm.W,dlm.C0,dlm.XX,0),'lik');
phiopt = fminsearch(ofun,phi)