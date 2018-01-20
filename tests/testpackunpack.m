%%
mo = manisetup(makeProduct(makeSE3Mat()));

zi = zeros(2,mo.group);
zi(:) = 1:numel(zi);
zu = maniunpack(mo,zi);
zp = manipack(mo,zu);
size(zu)
size(zp)


%%
mo = manisetup(makeproduct(makeRot(),makeRot()));

zi = zeros(2,mo.group);
zi(:) = 1:numel(zi);
zu = maniunpack(mo,zi);
zp = manipack(mo,zu);
size(zu)
size(zp)

zai = zeros(2,mo.alg);
zai(:) = 1:numel(zai);
zau = maniunpackalg(mo,zai);
zap = manipackalg(mo,zau);
size(zau)
size(zap)