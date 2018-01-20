addpath quaternions

%%
isquat = 0;

if isquat 
mz = manisetup(makeQuat());
else
mz = manisetup(makeRot());
end

mu0q = qomega2q([2,0,0]);
S0 = 8*[0.5 0 0.1; 0 0.2 0; 0.1 0 0.3];
wq = sampleQuats(1000,mu0q,S0); % generate 100 samples
if isquat == 0
    %w = arrayfun(@(x) reshape(q2dcm(w(x,:)),1,[]),1:size(wp,1),'UniformOutput','false');
    w = zeros(size(wq,1),9);
    for I=1:size(wq,1)
        w(I,:) = reshape(q2dcm(wq(I,:)),1,[]);
    end
    mu0 = q2dcm(mu0q);
else
    w = wq;
    mu0 = mu0q;
end

[muo,So,nv] = manimeancov(mz,w,20); % iterate 15 steps
if isquat == 0
    muo = reshape(muo,3,3);
    muoq = dcm2q(muo);
else
    muoq = muo;
end
deltamuo0 = quatdiff(muoq,mu0q)

S0
So

mu0
muo

subplot(3,1,1);
plot(nv(:,1));
xlabel('Ireration');
ylabel('Norm of step');
set(gca,'XTick',1:size(nv,1));
subplot(3,1,2);
plot(nv(:,2));
xlabel('Ireration');
ylabel('Residual Max');
set(gca,'XTick',1:size(nv,1));

subplot(3,1,3);
plot(nv(:,3));
xlabel('Ireration');
ylabel('Residual Mean');
set(gca,'XTick',1:size(nv,1));

