
# UKF over General Manifolds

We express first the UKF over a Manifold (Lie Group, SO3, SE3 or their combination or other) extending the existing works on the topic
(e.g. the original Quaternion UKF). 

In particular we acknowledge that

* state is in a Lie Group (or combination of there of): dimension n
* noise is expressed in the Tangent of the Lie Group: dimension m <= n
* we express a multivariate Gaussian as N(mu, Sigma) as usual

For general:
* delta(Group,Group) -> tangent
* step(Group,tangent) -> Group


* prod(Group,Group) -> Group

For Lie Group they are:
* prod(X,Y) = X * Y
* inv(X) se SO3 = R',, se SE3' [R' | -R't]
* delta(X, Y) = log(X * inv(Y))
* step(X, y) = exp(y) * X

First we define the sigma point construction and synthesis:

    Chi_i = step(mu,[ 0_alg +v1 ... +vm -v1...-vm ])
    mu = mean(Chi_i)
    Chi_not_mu_i = delta(Chi_i,mu)
    C  = 1/(2m+) Sum (Chi_not_mu_i' Chi_not_mu_i)

The mean(Chi_i) is an iterative algorithm:

    mu_0 = x0
    v_{i,k} = delta(x_i,mu_k)
    mu_{k+1} = step(mu_k,1/N Sum v_{i,k})

Then we deal with the Kalman aspect in particular the correction:

    x_i =  step(x_est_{i-1},K delta(z,z_{obs}))

# Representation

We have these representation:
- packed as vector
- tangent as vector
- unpacked as cell containing each element (e.g. 3x3 matrix for SO3 instead of 9 elements)

Example:
- SO3 quat: 4x1, 4, 3
- SO3 matrix: 3x3, 9, 3
- SE3: 4x4, 16, 7
- (SE3,SE3): {4x4,4x4}, 32, 14


# Examples of Manifold operations

Case of SO3
- group: Rotation 3x3
- algenra: vector R3
- product: R1 R2
- inverse: R'
- logarithm/exponential: rodriguez and its inverse

Case of quaternion
- group: quaternion 4
- algenra: vector R3
- product: q1 q2
- inverse: q*
- delta: axis angle between quaternions as vector (aka derivative)
- step: integral
- logarithm/exponential: rodriguez and its inverse

# TODO and Ideas
- test
- make example of two quat
- dual quaternion manifold vs SE3 matrix

# References

Riemman Manifold UKF 2013: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.646.893&rep=rep1&type=pdf
Hauberg, Søren, François Lauze, and Kim Steenstrup Pedersen. "Unscented Kalman filtering on Riemannian manifolds." Journal of mathematical imaging and vision 46.1 (2013): 103-120.

Quaternion UKF:
* Paper: http://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf (http://ieeexplore.ieee.org/document/1257247/)
* Our Paper: https://dl.dropboxusercontent.com/u/146279/papers/2017_C_ETFADiStefano.pdf

Note for Scaled UKF
* Bonus for large dimensions: http://www.cs.unc.edu/~welch/kalman/media/pdf/ACC02-IEEE1357.PDF

Note for Averaging
* Averaging: http://www.acsu.buffalo.edu/~johnc/uf_att.pdf
* Italian: https://re.public.polimi.it/retrieve/handle/11311/961634/40508/1_Manuscript.pdf

References for Beyond this:
- Bourmaud, Guillaume, et al. "Continuous-discrete extended Kalman filter on matrix Lie groups using concentrated Gaussian distributions." Journal of Mathematical Imaging and Vision 51.1 (2015): 209-228.
- Windle, Jesse, and Carlos M. Carvalho. "A tractable state-space model for symmetric positive-definite matrices." Bayesian Analysis 9.4 (2014): 759-792.
- Freifeld, Oren, Soren Hauberg, and Michael J. Black. "Model transport: Towards scalable transfer learning on manifolds." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2014.
- Hauberg, Søren. "Principal curves on Riemannian manifolds." IEEE transactions on pattern analysis and machine intelligence 38.9 (2016): 1915-1921.
- Srivatsan, Rangaprasad Arun, et al. "Estimating SE (3) elements using a dual quaternion based linear Kalman filter." Robotics: Science and Systems. 2016.

Geometric integration of quat:
- https://www.researchgate.net/profile/John_Crassidis/publication/260466470_Geometric_Integration_of_Quaternions/links/557acc0a08ae8d0481931b51/Geometric-Integration-of-Quaternions.pdf