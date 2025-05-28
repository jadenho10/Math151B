% QR–Iteration for Eigenvalues of A

A = [ 4   1  -2   2   0;
      1   2   0   1   3;
     -2   0   3  -2   1;
      2   1  -2  -1   4;
      0   3   1   4   5 ];

maxIter = 2;     % maximum number of iterations

Ak = A;             % initialize iterate
for k = 1:maxIter
    % QR decomposition of current iterate
    [Q, R] = qr(Ak);
    
    % form next iterate
    Ak = R*Q;
end

% Print the matrix Ak
disp(eig(A))
disp(Ak);

