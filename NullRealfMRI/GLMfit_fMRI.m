function GLMfit_fMRI(Y,X)

addpath /Users/sorooshafyouni/Home/BCF/ARMA/AR_YW

% Get residuals
pinvX = pinv(X);
R = eye(T)-X*pinvX;

for iv = 1:nv    
    resid(:,iv) = R*Y; 
end

rho = AR_YW_voxel(r,T,p);


% [Ainvt posdef]=chol(toeplitz([1 Coradj_pix]));
% nl=size(Ainvt,1);
% A=inv(Ainvt');
% B=ones(n-nl,1)*A(nl,:);
% Vmhalf=spdiags(B,1:nl,n-nl,n);

end
