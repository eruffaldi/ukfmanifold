clear all

if exist('manpippo.m','file')
    delete('manpippo.m'); % remove the fil
end
q = {makeQuat(),makeRn(3),makeRot()};
mout = makeGenProduct('pippo',q{:});
%,makeSE3Mat(),makeSE3Mat());

allinputs = rand(mout.group,100,2);
tic
s = 0;
for I=1:100
z1 = allinputs(:,I,1);
z1c = mout.unpack(z1);
z1c{1}= [1,0,0,0];
z1 = mout.pack(z1c);
z2 = allinputs(:,I,2);
z2c = mout.unpack(z2);
z2c{1} = z2c{1}/norm(z2c{1});
z2 = mout.pack(z2c);
z12 = mout.prod(z1,z2);
z12d = mout.delta(z12,z2);
z12s = mout.step(z1,z12d);
s = s+sum(z12s);
end
w =toc;
disp(sprintf('gen %f time %f',s,w))


mout = makeProduct(q{:});
tic
s = 0;
for I=1:100
z1 = allinputs(:,I,1);
z1c = mout.unpack(z1);
z1c{1}= [1,0,0,0];
z1 = mout.pack(z1c);
z2 = allinputs(:,I,2);
z2c = mout.unpack(z2);
z2c{1} = z2c{1}/norm(z2c{1});
z2 = mout.pack(z2c);
z12 = mout.prod(z1,z2);
z12d = mout.delta(z12,z2);
z12s = mout.step(z1,z12d);
s = s+sum(z12s);
end
w2 =toc;
disp(sprintf('nongen %f time %f',s,w2))
disp(sprintf('non generated is %f%% slower',(w2-w)/w*100));

%Z = mout.log(z);
%Zz = mout.exp(Z);