
moL = makeSE3MatLocal();
moG = makeSE3MatGlobal();

% note that moL.exp == moG.exp
m1 = moL.exp([1,0,2,0,3,4]);
m2 = moL.exp([0,0,1,5,0,0]);
mx = cellget(moL.unpack(m1),1)*cellget(moL.unpack(m2),1); % W..L1 L1..L2

aG = moG.delta(mx,m2); % so that: exp(aG)*m2 == mx ==> exp(aG) = mx*inv(m2)
m1G = moG.step(m2,aG); 


aL = moL.delta(mx,m2); % so that: m2*exp(aL) == mx ==> exp(aL) = inv(m2)*mx
m1L = moL.step(m2,aL);
m2p = cellget(moL.unpack(m2),1);
m1p = cellget(moL.unpack(m1),1);
m1Lp = cellget(moL.unpack(m1L),1);
m1Gp = cellget(moL.unpack(m1G),1);
aLG = (se3adj(m2p)*aL')';
m1p
m1Lp
m1Gp

aL
aG
aLG


%%
se3adj(cellget(moL.unpack(moL.exp([0,0,0, 0,0,0])),1))*[0,0,0,0,2,0]'
