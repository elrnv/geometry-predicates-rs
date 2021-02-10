# geometry-predicates

A safe Rust port of ["Adaptive Precision Floating-Point Arithmetic and Fast Robust
Predicates for Computational Geometry"](https://www.cs.cmu.edu/~quake/robust.html) 

[![Build Status](https://github.com/elrnv/geometry-predicates-rs/workflows/Test/badge.svg)](https://github.com/elrnv/geometry-predicates-rs/actions)
[![Documentation](https://docs.rs/geometry-predicates/badge.svg)](https://docs.rs/geometry-predicates)
[![Version](https://img.shields.io/crates/v/geometry-predicates.svg)](https://crates.io/crates/geometry-predicates)
[![License](https://img.shields.io/crates/l/geometry-predicates.svg)](https://github.com/elrnv/geometry-predicates-rs/blob/master/LICENSE)
[![Downloads](https://img.shields.io/crates/d/geometry-predicates.svg)](https://crates.io/crates/geometry-predicates)

This library provides a Rust solution to efficient exact geometry predicates
used widely for computational geometry.

The following predicates are exposed at the root level:
 - orient2d
 - incircle
 - orient3d
 - insphere

Along with their approximate (inexact) counterparts:
 - orient2d_fast 
 - incircle_fast
 - orient3d_fast
 - insphere_fast

In addition, the building blocks of these predicates, namely the adaptive precision
floating-point arithmetic primitives are also exposed to allow for extensions to
other predicates or exact geometric constructions.

This library supports `no-std` targets with standard IEEE 754 floats.


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

Mathematically, predicates like `orient2d` are
defined as
```
                                        ⎛⎡ax ay 1⎤⎞
orient2d([ax,ay],[bx,by],[cx,cy]) := det⎜⎢bx by 1⎥⎟
                                        ⎝⎣cx cy 1⎦⎠
```

It's easy to see that these predicates solve the problem of
computing the determinant of small matrices with the correct sign, regardless of how
close the matrix is to being singular.

For instance to compute the determinant of a matrix `[a b; c d]` with the
correct sign, we can invoke
```
orient2d([a,b], [c,d], [0,0])
```

For more details please refer to Jonathan Shewchuk's [original
webpage](https://www.cs.cmu.edu/~quake/robust.html) for these predicates.


## Caveats

These predicates do NOT handle exponent overflow [\[1\]], which means inputs with floats smaller than
`1e-142` or larger than `1e201` may not produce accurate results. This is true for the original
predicates in `predicates.c` as well as other Rust ports and bindings for these predicates.


## References

 - [\[1\] Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates][\[1\]],
   Discrete & Computational Geometry 18(3):305–363, October 1997.
 - [\[2\] Robust Adaptive Floating-Point Geometric Predicates Proceedings of the Twelfth Annual][\[2\]],
   Symposium on Computational Geometry (Philadelphia, Pennsylvania), pages 141–150, Association for
   Computing Machinery, May 1996


## Acknowledgements

This port was originally created by a C to Rust translator called
[Corrode](https://github.com/jameysharp/corrode). Without it, a full Rust port
of this library would have been a daunting task. With that I would specifically like to thank the
authors of Corrode for providing such a useful tool.
Version 0.2 of this crate used the [c2rust](https://c2rust.com/) crate. The same gratitude goes towards
the developers of C2Rust.

[\[1\]]: http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf
[\[2\]]: http://www.cs.berkeley.edu/~jrs/papers/robust-predicates.pdf

