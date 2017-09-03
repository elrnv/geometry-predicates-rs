# geometry-predicates

[![Build Status](https://travis-ci.org/elrnv/geometry-predicates-rs.svg?branch=master)](https://travis-ci.org/elrnv/geometry-predicates-rs)

A Rust port of ["Adaptive Precision Floating-Point Arithmetic and Fast Robust
Predicates for Computational Geometry"](https://www.cs.cmu.edu/~quake/robust.html) 

This library provides a Rust solution to efficient exact geometry predicates
used widely for computational geometry.

In addition, the building blocks of these predicates, namely the adaptive precision
floating-point arithmetic primitives are also exposed to allow for extensions to
other predicates or exact geometric constructions.o

So far, the following predicates are supported:
 - orient2d
 - orient3d

## Background

These predicates have been a staple in computational geometry for many years now
and are widely used in industry.   In the context of geometry algorithms, it is
often essential to determine the orientation of a point with respect to a line (or a
plane) and whether a point lies inside a circle (or a sphere) or not.  The reason
why these tests often need to be exact is because many geometry algorithms
ask questions (to determine orientation or in-circle/sphere) about point
configurations that require consistent answers. For instance, if `a`, `b`, and
`c` are three points on a 2D plane, to ask where `b` with respect to the line
through `a` and `c` (left-of, right-of, or coincident) is the same as asking where
`a` lies with respect to the line through `c` and `b`.
Formally this condition can be written as
```
sgn(orient2d(a,c,b)) == sgn(orient2d(c,b,a))
```

Mathematically (using MATLAB-style notation), predicates like `orient2d` are
defined as
```
orient2d([ax,ay],[bx by],[cx cy]) := det([ax ay 1; bx by 1; cx cy 1])
```

It's easy to see that these predicates solve the problem of
computing the determinant of small matrices with the correct sign, regardless of how
close the matrix is to being singular.

For instance to compute the determinant of a matrix `[a b; c d]` with the
correct sign, we can invoke `orient2d([a,b], [c,d], [0,0])`.

For more details please refer to the [original
webpage](https://www.cs.cmu.edu/~quake/robust.html) for these predicates.

## Acknowledgements

This port was created by a C to Rust translator called
[Corrode](https://github.com/jameysharp/corrode). Without it, a full Rust port
of this library would have been a daunting task. With that I would specifically like to thank the
authors of Corrode for providing such a useful tool.
