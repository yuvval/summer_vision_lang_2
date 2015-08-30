A = zeros(2, 94, 10);
combs = allcomb(1:2,1:94);
vals = 1:(2*94);
A(:,:,1) = sparse(combs(:,1), combs(:,2), vals);


