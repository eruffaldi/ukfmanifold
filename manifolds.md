
Basic manifolds:

- makeRn: euclidean R^n
- makeQuat: SO3 quaternion
- makeRot: SO3 matrix 3x3

Operations:
- makeProduct(m1...mK): combines the manifolds
- makePower(m,n): takes the manifold m and replicates n times: 

In Progress:
- makeSE3Mat: SE3 as matrix 4x4
- makeQuat1: SO3 constrained to 1 axis of rotation, quaternion

TODO
- makeRot1: as above but with 3x3 matrix
- makePositive: positive scalar (very specific case of positive definite matrices)
- makeSE3Quat: SE3 as quaternion and position

 