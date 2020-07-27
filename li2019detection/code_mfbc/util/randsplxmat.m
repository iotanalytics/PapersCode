function U = randsplxmat(m,n)
U = -log(rand(m,n));
S = sum(U,1); % probably only 1 row, so specify 1 explicitly 
for j = 1:n
    U(:,j) = U(:,j) /  S(j);   
end
end