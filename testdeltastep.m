mysetup('quaternions')
%%
m = makeQuat();

q0 = qomega2q([0.5,0.2,0]);
w1 = [0.2,0,0];
qs = m.step(q0,w1);
w2 = m.delta(qs,q0)
w1-w2


%%
m = makeSE3Mat();

x0 = m.exp([0.5,0.2,0,  0.1,0,1]); % 4x4
iix0 = m.inv(m.inv(x0));

w1 = [0.2,0,0, 0.0,0.1,2];
x1 = m.step(x0,w1);
w2 = m.delta(x1,x0)
w1-w2

%%
mr = makeRot();

R0 = mr.exp([0.5,0.2,0]);
w0 = mr.log(R0);
R0u = det(mr.unpack(R0));

w1 = [0.2,0,0];
Rs = mr.step(R0,w1);
w2 = mr.delta(Rs,R0)
w1-w2