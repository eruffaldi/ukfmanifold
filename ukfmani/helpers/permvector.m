function P = permvector(neworder)

assert(length(unique(neworder)) == length(neworder),'all elements');
assert(all(neworder >= 1)&all(neworder <= length(neworder)),'all elments 1..n');

P = zeros(length(neworder));
for I=1:length(neworder)
    P(I,neworder(I))= 1;
end
