function P = polyfitM(matrixX, matrixY, n)
[sx1, sx2] = size(matrixX);
P = zeros(sx1, n+1);
for Index = 1:sx1
  x = matrixX(Index, :);
  y = matrixY(Index, :);
  x = x(:);
  y = y(:);
  q = ones(sx1, 1);
  V(:, n+1) = q;
  for j = n:-1:1
    q = q .* x;
    v(:, j) = q;
  end
  [Q, R] = qr(V, 0);
  %P(Index, :) = (R \ (Q' * y))';
  transpose(R \ (transpose(Q) * y));
end

end