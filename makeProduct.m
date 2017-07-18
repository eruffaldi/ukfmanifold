% takes a list of manifolds and builds a single manifold
%
% Emanuele Ruffaldi 2017 @ SSSA
function m = makeProduct(varargin)



if nargin == 1
    m = varargin{1};
elseif nargin == 2
    m1 = varargin{1};
    m2 = varargin{2};

    g = m1.group;
    a = m1.alg;

    m = [];
    m.inv = @(x) [m1.inv(x(1:g)), m2.inv(x(g+1:end))];
    m.prod = @(x1,x2) [m1.prod(x1(1:g),x2(1:g)), m2.prod(x1(g+1:end),x2(g+1:end))];
    if isfield(m1,'log') && isfield(m2,'log')
        m.log = @(x) [...
                m1.log(x(1:g)), ...
                m2.log(x(g+1:end))];
        m.exp = @(x) [...
                m1.exp(x(1:a)), ...
                m2.exp(x(a+1:end))];
    end
    m.delta = @(x,y) [m1.delta(x(1:g),y(1:g)), m2.delta(x(g+1:end),y(g+1:end))];
    m.step =@(X,y) [m1.step(X(1:g),y(1:a)), m2.step(X(g+1:end),y(a+1:end))]; 
    m.group = m1.group+m2.group;
    m.alg = m1.alg+m2.alg;
    m.count = m1.count+m2.count;
    m.models = {m1,m2};
    m.pack = @(x) [m1.pack(x{1}), m2.pack(x{2})]; % this works with cells
    m.unpack = @(x) {m1.unpack(x(1:g)), m2.unpack(x(g+1:end))};


else
    if 1==0
        warning('makeCom multiple to be completed');
        m = makeProduct(varargin{1},varargin{2});

        for J=3:nargin
            m = makeProduct(m,varargin{J});
        end

    else

        m = [];

        % used for speeding up the partitioning operations
        m.groupinc = [0,cumsum(cellfun(@(x) x.group,varargin))];
        m.alginc = [0,cumsum(cellfun(@(x) x.alg,varargin))];

        m.group = sum(cellfun(@(x) x.group,varargin));
        m.alg = sum(cellfun(@(x) x.alg,varargin));
        m.count = sum(cellfun(@(x) x.count,varargin));
        m.models = varargin;
        m.inv = @(x) cominv(m,x);
        m.prod = @(x1,x2) comprod(m,x1,x2);

        if all(cellfun(@(x) isfield(x,'log'),varargin))
            m.log = @(x) comlog(m,x);
            m.exp = @(x) comexp(m,x);
        end

        
        m.delta = @(x,y) comdelta(m,x,y);
        m.step =  @(X,y) comstep(m,X,y);
        m.pack = @(x) compack(m,x);
        m.unpack = @(x) comunpack(m,x);
    end

end

function z = comprod(m,x1,x2)
    z = zeros(m.group,1);
    for I=1:length(m.models)
        ab = m.groupinc(I)+1:m.groupinc(I+1);
        w2 = x2(ab);
        w1 = x1(ab);
        z(ab) = m.models{I}.prod(w1(:)',w2(:)');
    end

function z = cominv(m,X)

    z = zeros(m.group,1);
    for I=1:length(m.models)
        ab = m.groupinc(I)+1:m.groupinc(I+1); 
        z(ab) = m.models{I}.inv(X(ab));
    end

function z = comlog(m,X)

    z = zeros(m.alg,1);
    for I=1:length(m.models)
        ab = m.groupinc(I)+1:m.groupinc(I+1); 
        aab = m.alginc(I)+1:m.alginc(I+1); 
        z(aab) = m.models{I}.log(X(ab));
    end

function z = comexp(m,x)

    z = zeros(m.group,1);
    for I=1:length(m.models)
        ab = m.groupinc(I)+1:m.groupinc(I+1); 
        aab = m.alginc(I)+1:m.alginc(I+1); 
        z(ab) = m.models{I}.exp(x(aab));
    end

function z = comstep(m,X,y)

    z = zeros(m.group,1);
    for I=1:length(m.models)
        ab = m.groupinc(I)+1:m.groupinc(I+1); 
        aab = m.alginc(I)+1:m.alginc(I+1); 
        z(ab) = m.models{I}.step(X(ab)',y(aab));
    end
    
function z = comdelta(m,X,Y)

    z = zeros(m.group,1);
    for I=1:length(m.models)
        ab = m.groupinc(I)+1:m.groupinc(I+1); 
        z(ab) = m.models{I}.delta(X(ab)',Y(ab));
    end

function y = compack(m,X)

    z = zeros(m.group,1);
    for I=1:length(m.models)
        ab = m.groupinc(I)+1:m.groupinc(I+1); 
        z(ab) = m.models{I}.pack(X{I});
    end


function z = comunpack(m,x)

    z = cell(m.models,1);
    for I=1:length(m.models)
        ab = m.groupinc(I)+1:m.groupinc(I+1); 
        z{I} = m.models{I}.unpack(x(ab));
    end



