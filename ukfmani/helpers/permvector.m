function P = permvector(neworder,n)
if nargin == 2
    assert(isscalar(n),'n is the number of elements');
    left = setdiff(1:n,neworder);
    neworder = [neworder(:);left(:)];
end

assert(length(unique(neworder)) == length(neworder),'all elements');
assert(all(neworder >= 1)&all(neworder <= length(neworder)),'all elments 1..n');

P = zeros(length(neworder));
for I=1:length(neworder)
    P(I,neworder(I))= 1;
end
