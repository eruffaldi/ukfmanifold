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
elseif length(Cs) == 2
    [X,C] = manifusion(model,mus{1},Cs{1});
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
            v(i,:) = Csi{I}*model.delta(mus{i},mz);
        end
        mz = model.step(mz,sum(v,1)); % averaging the weighted directions
    end
end

% inv(inv(C0)+inv(C1))
C = C0 - C0/(C0+C1)*C0;

v = m.delta(x1,x0); % from x0 to x1

X = m.step(x0,C/C1*v); % idem but weightefunction [X,C] = manifusion(m,x0,x1,C0,C1)

% inv(inv(C0)+inv(C1))
C = C0 - C0/(C0+C1)*C0;

v = m.delta(x1,x0); % from x0 to x1

X = m.step(x0,C/C1*v); % idem but weightedd