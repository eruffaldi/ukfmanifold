
function Y = meannonan(X,dim)

if nargin == 1
    dim = find(size(X) > 1,1,'first');
    if isempty(dim)
        dim = 1;
    end
end

if dim == 1
    Y = zeros(size(X,2),1);
    for I=1:size(Y,1)
        Y(I) = mean(X(isnan(X(:,I)) == 0,I));
    end
else
    Y = zeros(size(X,1),1);
    for I=1:size(Y,1)
        Y(I) = mean(X(I,isnan(X(I,:)) == 0));
    end
end
