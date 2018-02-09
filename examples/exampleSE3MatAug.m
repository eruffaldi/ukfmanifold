% We build first the manifolds explicitly stating the state manifold x (Mx)
% and the augmented Mxa
Mx = {makeSE3Mat(),makeRn(3),makeRn(3)};
Mxa = [Mx,makeRn(6)];
Mz = makeSE3Mat();

% Then we proceed to the setup of the manifold in usable way, also using
% the code generation for speedup of the processing
mx = manisetup(makeGenProduct('se3mat_e3_e3',Mx{:}));
mxa = manisetup(makeGenProduct('se3mat_e3_e3_aug',Mxa{:})); % as makeProduct(....,+1)
mz = manisetup(Mz);
mzr = manisetup(makeRot()); % only rotation, loss of position
mxt = makeSE3Mat(); % helper without the need of setup

% We build these structures for supporting the Unscented transformation. We
% need them for the manifold x in the process, and the augmented xa in the
% observation
mx.wsigma = ut_mweights2(mx.group,mx.alg,0.5);
mx.wsigma.sqrt = @svdsqrt; 
mxa.wsigma = ut_mweights2(mxa.group,mxa.alg,0.5);
mxa.wsigma.sqrt = @svdsqrt; 


% initial state and noise definition
x0 = mx.step(mx.exp([0,0,0,  0,1,0,   0,0,0,   0,0,0]'),[pi/2,0.2,0,  0,0,0,   0,0,0,   0,0,0]');
P0 = 0.5*eye(mx.alg);
Q = 0.01*eye(mx.alg); % process noi!se
R = 1e-3*eye(mz.alg); % measure noise
Rr = R(1:3,1:3); % only rotations


%% Generate the Observation by simulation
z0 = mz.exp([pi/2,0,0, 0,0.1,1]');
zobsval = zeros(16,200);
v0 = zeros(6,1);
v0(4) = 0.1;
v0(1) = 0.2;
zobsval(:,1) = z0;
lzobsval = zeros(mz.alg,size(zobsval,2));
lzobsval(:,1) = mz.log(z0);
for I=2:size(zobsval,2)
    v0(2) = sin(I/100);
    zobsval(:,I) = mz.step(zobsval(:,I-1),v0);
    lzobsval(:,I) = mz.log(zobsval(:,I));
end
zobs = @(t) zobsval(:,t);


% observation is identity
% process is the integral
dt = 0.1;

%% Kalman Setup
% functions work expanding each primitive manifold with their type (e.g. matrix 4x4) 
% if the output is made of multiple e manifolds use deal for returning
% everything
f_fx = @(Tk,wk,vk) deal(mxt.step(Tk,[wk;vk]),wk,vk); % Xk = (Tk,wk,vk)
h_fx = @(Tk,wk,vk,ek) mxt.step(Tk,ek);
hr_fx = @(Tk,wk,vk,ek) mxt.step(Tk(1:3,1:3),ek(1:3));
tic
% loop
deltas = zeros(200,mz.alg);
states = zeros(size(deltas,1),mx.group);
lstates = zeros(size(deltas,1),mx.alg);
usereduxspace = 0;
Rnonadditive = eye(mxa.alg-mx.alg); % generic formulation
Rnonadditive = R;
R = zeros(size(R)); % for testing augmentation

for L=1:size(deltas,1)
    states(L,:) = x0;
    lstates(L,:) = mx.log(x0);  % [ rodriguez(R) traslazione velocitàangolar velocitàlinear ]
    
    [xp,Pp] = manistatestep(mx,x0,P0,f_fx,Q);
    assert(size(xp,2)==1);

    
    % Build the augmented state by expanding the state with noise with
    % zero mean and provided covariance
    Ppa = blkdiag(Pp,Rnonadditive);
    xpa = [xp; zeros(mxa.alg-mx.alg,1)]; % mean is zero
    [zm,Czz,Cxaz] = manievalh(mxa,mz,xpa,Ppa,h_fx);
        

    Pvv = Czz + R;
    z = zobs(L);

    K = Cxaz/Pvv;
    Ppa = (Ppa - K * Pvv * K');
    delta = mz.delta(z,zm);
    
    xpanew = mx.step(xpa,(K*delta'));
    % extract x0 from xpanew because we discard the augmentation
    % extract P0 from Ppa because we discard augmentation
    x0 = xpanew(1:mx.group);
    P0 = Ppa(1:mx.alg,1:mx.alg);
       
    deltas(L,:) = delta;
end
%%
toc
figure(1)
subplot(3,1,1);
plot(deltas(:,1:3))
title('Deltas observation-prediction');
xlabel('sample');
subplot(3,1,2);
plot(deltas(:,4:6))
title('Deltas observation-prediction');
xlabel('sample');
subplot(3,1,3);
plot(sum(deltas.^2,2));
title('Norm of error');
xlabel('sample');
figure(2)
plot(states(10:end,:))
title('All states as matrix');

%% Latency estimation
dd = zeros(1,6);
for J=1:6
    figure(2+J);
    dd(J) = finddelay(lstates(J,:),lzobsval(:,J));

    plot([lstates(J,:)';lzobsval(:,J)]);
end
disp('delay')
dd
disp('error')
sqrt(meannonan(deltas.^2,1))'
