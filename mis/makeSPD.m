function K = makeSPD(K)

    [U,D] = eig(K);   % calculate the eigenvalues and eigenvectors of the GRM
    if min(diag(D)) < 0   % check whether the GRM is non-negative definite
        disp('WARNING: the GRM is not non-negative definite! Set negative eigenvalues to zero')
        D(D<0) = 0;   % set negative eigenvalues to zero
        K = U*D/U;   % reconstruct the GRM
    end

end