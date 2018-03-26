% fusion 
function [X,C] = manifusionn(model,steps,mus,Cs)

if isempty(steps)
    steps = 10;
end

assert(length(mus) > 0);
assert(iscell(mus));
assert(iscell(Cs));
assert(length(mus) == length(Cs));

if length(Cs) == 1
    X = mus{1};
    C = Cs{1};
    return;
elseif length(Cs) == 2
    [X,C] = manifusion(model,mus{1},Cs{1});
    return;
else
    C = zeros(size(Cs{1}));
    for I=1:length(C)
        C = C + inv(Cs{I});
    end
    C = inv(C);
    Csi = cell(size(Cs));
    for I=1:length(Csi)
        Csi{I} = C/Cs{I};
    end
        
    
    steps = 20;
    N=length(mus);

    v = zeros(N,model.alg); % preallocated
    mz = mus{1}; % first is picked as mean

    % iterative 
    for k=1:steps
        for i=1:N
            % same as: se3_logdelta but with igk once
            v(i,:) = Csi{I}*model.delta(mus{i},mz); % Csi{i}?
        end
%        mz = model.step(mz,mean(v,1));
        mz = model.step(mz,mean(v,1)'); % averaging the weighted directions
    end
    X = mz;
    
   % update v for computing covariance
    for i=1:N
%        v(i,:) = Csi{I}*model.delta(mus(i,:),X);
        v(i,:) = Csi{I}*model.delta(mus{i},X);
    end

    C = v'*v; % covariance ZZ
end
