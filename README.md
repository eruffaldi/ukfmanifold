
# Unscented Kalman Filtering (UKF) over Manifolds

This Matlab toolbox aims at supporting non-linear transformations of Multivariate Normal (MVN) variables and bayesian filtering of such variable. The toolbox relies on the concept of Unscented Transformation that is at the basis of the Unscented Kalman Filtering. For a visual comparison between Extended Kalman Filtring and Unscented Kalman filtering look at https://github.com/eruffaldi/compare-mvn-transform

Look at related publications supporting the theory of this: Quaternion UKF [1], UKF over Riemman Manifolds [2] and manifolds [3].

# Concepts

* state is in a manfold (e.g Lie Group), or combination of there of with dimension n
* noise is expressed in the tangent space (e.g. Lie Algebra): dimension m <= n
* we express a multivariate Gaussian as N(mu, Sigma) with mean over manifold and covariance tangent space

General manifolds requires this:

* delta(Group,Group) -> tangent
* step(Group,tangent) -> Group

Lie Groups require this:

* prod(X,Y) = X * Y
* inv(X)

And Lie Groups define delta and steps as:

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

## Representations

We have these representation of the manifold data

- packed as vector [G,1]: each manifold is packed as a vector, not in the minimal representation (e.g. a matrix 3x3 is 9x1, matrix 4x4 is 16x1)
- tangent as vector [A,1]
- expanded as matrix or vector (e.g. 3x3 for matrix of rotation). Products of manifolds are expressed as cell array with the contained type, eventually in nested form

When considering multiple samples we use matlab convention: [G,N] and [A,N] and cell [C,N]

%TODO: for nested manifold decide if nested cell arrays or keep cell arrays expanded

Example:
- SO3 quaternion: expanded group 4x1, packed group 4x1, algebra 3
- SO3 matrix: expanded group 3x3, packed group 9, algebra 3
- SE3 matrix: expanded group 4x4, packed group 16, algebra 6
- (SE3,SE3): expanded group {4x4,4x4}, packed group 32, algebra 12


Each Manifold is represented as a struct with:
- group as the size of the packed group
- alg as the size of the algrebra space
- count as the number of contained native manifolds, corresponding to the cell size
- function handles for each operation
- models as sub-models if present

# Usage

Initialize the toolbox running init.m 

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

# Code Generation

The makeGenProduct(name,m1...mN) produces a Matlab code that contains the efficient composition of the provided manifolds. In this way the singl operations are executed efficiently (e.g. testmakegen for a test) 

On my tests a composition of 4 manifolds give 30% speed up

TODO: support nested manifolds code generation

# Generate Documentation

- Sphinx
- m2html

m2html('mfiles','ukfmani', 'htmldir','doc', 'recursive','on', 'global','on','todo','on');%,'source','off');

# References

[1] Quaternion UKF:  http://kodlab.seas.upenn.edu/uploads/Arun/UKFpaper.pdf (http://ieeexplore.ieee.org/document/1257247/)
Kraft, E. (2003, July). A quaternion-based unscented Kalman filter for orientation tracking. In Proceedings of the Sixth International Conference of Information Fusion (Vol. 1, pp. 47-54).

[2] Riemman Manifold UKF 2013: http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.646.893&rep=rep1&type=pdf
Hauberg, Søren, François Lauze, and Kim Steenstrup Pedersen. "Unscented Kalman filtering on Riemannian manifolds." Journal of mathematical imaging and vision 46.1 (2013): 103-120.

[3] Hertzberg, Christoph, et al. "Integrating generic sensor fusion algorithms with sound state representations through encapsulation of manifolds." Information Fusion 14.1 (2013): 57-77.

[4] Our variant for robotics: http://www.eruffaldi.com/papers/2017_C_ETFADiStefano.pdf
Di Stefano, Erika, Emanuele Ruffaldi, and Carlo Alberto Avizzano. "A Multi-Camera Framework for Visual Servoing of a Collaborative Robot in Industrial Environments."

Note for Scaled UKF
* Bonus for large dimensions: http://www.cs.unc.edu/~welch/kalman/media/pdf/ACC02-IEEE1357.PDF

Note for Averaging
* Averaging: http://www.acsu.buffalo.edu/~johnc/uf_att.pdf
* Italian: https://re.public.polimi.it/retrieve/handle/11311/961634/40508/1_Manuscript.pdf

References for research starting from this:
- Bourmaud, Guillaume, et al. "Continuous-discrete extended Kalman filter on matrix Lie groups using concentrated Gaussian distributions." Journal of Mathematical Imaging and Vision 51.1 (2015): 209-228.
- Windle, Jesse, and Carlos M. Carvalho. "A tractable state-space model for symmetric positive-definite matrices." Bayesian Analysis 9.4 (2014): 759-792.
- Freifeld, Oren, Soren Hauberg, and Michael J. Black. "Model transport: Towards scalable transfer learning on manifolds." Proceedings of the IEEE Conference on Computer Vision and Pattern Recognition. 2014.
- Hauberg, Søren. "Principal curves on Riemannian manifolds." IEEE transactions on pattern analysis and machine intelligence 38.9 (2016): 1915-1921.
- Srivatsan, Rangaprasad Arun, et al. "Estimating SE (3) elements using a dual quaternion based linear Kalman filter." Robotics: Science and Systems. 2016.

Geometric integration of quat:
- https://www.researchgate.net/profile/John_Crassidis/publication/260466470_Geometric_Integration_of_Quaternions/links/557acc0a08ae8d0481931b51/Geometric-Integration-of-Quaternions.pdf