//! A safe Rust port of the [robust adaptive floating-point geometric predicates](https://www.cs.cmu.edu/~quake/robust.html).
//!
//! This crate provides a Rust solution to efficient exact geometry predicates
//! used widely for computational geometry.
//!
//! In addition, the building blocks of these predicates, namely the adaptive precision
//! floating-point arithmetic primitives, are also exposed in [`predicates`] to allow for extensions
//! to other predicates or exact geometric constructions.
//!
//! ## Background
//!
//! These predicates have been a staple in computational geometry for many years now
//! and are widely used in industry.   In the context of geometry algorithms, it is
//! often essential to determine the orientation of a point with respect to a line (or a
//! plane) and whether a point lies inside a circle (or a sphere) or not.  The reason
//! why these tests often need to be exact is because many geometry algorithms
//! ask questions (to determine orientation or in-circle/sphere) about point
//! configurations that require consistent answers. For instance, if `a`, `b`, and
//! `c` are three points on a 2D plane, to ask where `b` with respect to the line
//! through `a` and `c` (left-of, right-of, or coincident) is the same as asking where
//! `a` lies with respect to the line through `c` and `b`.
//! In Rust, this condition can be written as
//! ```
//! # use geometry_predicates::orient2d;
//! # let a = [1.0, 2.0];
//! # let b = [3.0, 4.0];
//! # let c = [5.0, 6.0];
//! assert_eq!(orient2d(a,c,b).signum(), orient2d(c,b,a).signum());
//! ```
//!
//! Mathematically, predicates like `orient2d` are
//! defined as
//!```verbatim
//!                                        ⎛⎡ax ay 1⎤⎞
//!orient2d([ax,ay],[bx,by],[cx,cy]) := det⎜⎢bx by 1⎥⎟
//!                                        ⎝⎣cx cy 1⎦⎠
//!```
//!
//! It's easy to see that these predicates solve the problem of
//! computing the determinant of small matrices with the correct sign, regardless of how
//! close the matrix is to being singular.
//!
//! For instance to compute the determinant of a matrix `[a b; c d]` with the
//! correct sign, we can invoke
//! ```
//! # use geometry_predicates::orient2d;
//! # let a = 1.0;
//! # let b = 2.0;
//! # let c = 3.0;
//! # let d = 4.0;
//! assert_eq!(orient2d([a,b], [c,d], [0.0,0.0]), a*d - b*c);
//! ```
//!
//! For more details please refer to the [original
//! webpage](https://www.cs.cmu.edu/~quake/robust.html) for these predicates.
//!
//! ## Caveats
//!
//! These predicates do NOT handle exponent overflow [\[1\]], which means inputs with floats smaller than
//! `1e-142` or larger than `1e201` may not produce accurate results. This is true for the original
//! predicates in `predicates.c` as well as other Rust ports and bindings for these predicates.
//!
//! ## References
//!
//!  - [\[1\] Adaptive Precision Floating-Point Arithmetic and Fast Robust Geometric Predicates][\[1\]],
//!    Discrete & Computational Geometry 18(3):305–363, October 1997.
//!  - [\[2\] Robust Adaptive Floating-Point Geometric Predicates Proceedings of the Twelfth Annual][\[2\]],
//!    Symposium on Computational Geometry (Philadelphia, Pennsylvania), pages 141–150, Association for
//!    Computing Machinery, May 1996
//!
//! [\[1\]]: http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf
//! [\[2\]]: http://www.cs.berkeley.edu/~jrs/papers/robust-predicates.pdf
#![no_std]
pub mod predicates;

#[cfg(feature = "transpiled")]
mod transpiled;

// The following predicates are exposed at the top level.
// There are other alternative robust implementations (but typically slower) of these predicates.
// We use those to check the adaptive implementations.
pub use predicates::{
    // Adaptive robust predicates.
    incircle,
    // Fast inexact predicates.
    incircle_fast,
    insphere,
    insphere_fast,
    orient2d,
    orient2d_fast,
    orient3d,
    orient3d_fast,
};

#[cfg(test)]
mod tests {
    extern crate rand;
    use self::rand::{Rng, SeedableRng, StdRng};
    use super::*;

    use predicates::{
        incircle_exact, incircle_slow, insphere_exact, insphere_slow, orient2d_exact,
        orient2d_slow, orient3d_exact, orient3d_slow,
    };

    const SEED: &[usize] = &[1, 2, 3, 4];

    /*
     * Note on robustness testing.
     *
     * These predicates do NOT handle overflow or underflowo of the exponent.
     * Quoting Shechuk's predicates paper directly:
     *
     * > The four predicates implemented [in this library] will not overflow nor underflow if their
     * > inputs have exponents in the range [-142, 201] and IEEE 754 double precision arithmetic
     * > is used.
     * - Jonathan Shechuk 1997
     *
     * We will therefore be careful to test for inputs with exponents in that range and not beyond.
     */

    const EXP_BOUNDS: [i32; 2] = [-142, 201];

    /// Generate a tolerance value with exponent no less than -142.
    fn tol(rng: &mut StdRng) -> f64 {
        let max_exp = (EXP_BOUNDS[0] + 1) as f64;
        10.0_f64.powi((max_exp * rng.gen::<f64>()).round() as i32) * (rng.gen::<f64>() - 0.5)
    }

    /*
     * Many of the tests below ensure that the predicates produce expected results for all
     * permutations of the inputs. Thus, below we write an adhoc QuickPerm implementation to
     * generate swap based permutations.
     *
     * The following is a basic countdown implementation of QuickPerm (see quickperm.org).
     * This implementation can be refactored to use const generics when those stabilize.
     */
    macro_rules! quick_perm_impl {
        ($n:expr, $struct_name:ident) => {
            struct $struct_name<T> {
                perm: [T; $n],
                p: [usize; $n + 1],
                index: usize,
            }

            impl<T> $struct_name<T> {
                fn new(v: [T; $n]) -> Self {
                    let mut p = [0; $n + 1];
                    for i in 0..$n + 1 {
                        p[i] = i;
                    }
                    $struct_name {
                        perm: v,
                        p,
                        index: 0,
                    }
                }
            }

            impl<T: Clone> Iterator for $struct_name<T> {
                type Item = [T; $n];
                fn next(&mut self) -> Option<Self::Item> {
                    let $struct_name { perm, p, index } = self;
                    let mut i = *index;
                    let n = p.len() - 1;
                    if i == 0 {
                        *index += 1;
                        return Some(perm.clone());
                    } else if i >= n {
                        return None;
                    }
                    p[i] -= 1;
                    let j = if i % 2 == 0 { 0 } else { p[i] };
                    perm.swap(j, i);
                    i = 1;
                    while p[i] == 0 {
                        p[i] = i;
                        i += 1;
                    }
                    *index = i;
                    Some(perm.clone())
                }
            }
        };
    }

    quick_perm_impl!(3, QuickPerm3);
    quick_perm_impl!(4, QuickPerm4);
    quick_perm_impl!(5, QuickPerm5);

    /*
     * The tests below have specific purposes.
     *
     * - _robustness_ tests check that the adaptive predicates work in degenerate or close to degenerate
     *   configurations
     * - _regression_ tests verify the adpative predicates against their _slow and _exact variants.
     * - _transpiled_regression_ tests verify the predicates against the directly transpiled
     *   (unrefactored) version of the library.
     */

    #[test]
    fn orient2d_test() {
        let a = [0.0, 1.0];
        let b = [2.0, 3.0];
        let c = [4.0, 5.0];
        assert_eq!(orient2d(a, b, c), 0.0);
    }

    #[cfg(feature = "transpiled")]
    #[test]
    fn orient2d_transpiled_regression_test() {
        unsafe {
            transpiled::exactinit();
        }
        let mut rng: StdRng = SeedableRng::from_seed(SEED);

        let n = 99999;
        for _ in 0..n {
            let a = [tol(&mut rng), tol(&mut rng)];
            let b = [12.0, 12.0];
            let c = [24.0, 24.0];
            let o2d = orient2d;
            let o2d_transpiled = transpiled::orient2d;
            for p in QuickPerm3::new([a, b, c]) {
                assert_eq!(
                    o2d(p[0], p[1], p[2]),
                    o2d_transpiled(p[0], p[1], p[2]),
                    "{:?}",
                    a
                );
            }
        }
    }

    #[test]
    fn orient2d_robustness_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let n = 99999;

        for o2d in &[orient2d, orient2d_exact, orient2d_slow] {
            for _ in 0..n {
                let pa = [tol(&mut rng), tol(&mut rng)];
                let pb = [12.0, 12.0];
                let pc = [24.0, 24.0];

                let main = o2d(pa, pb, pc);
                if main == 0.0 {
                    for [a, b, c] in QuickPerm3::new([pa, pb, pc]) {
                        assert_eq!(o2d(a, b, c), 0.0);
                    }
                }

                let pred = main > 0.0;
                for (i, [a, b, c]) in QuickPerm3::new([pa, pb, pc]).enumerate() {
                    let t = o2d(a, b, c);
                    assert_eq!(pred, if i % 2 == 0 { t > 0.0 } else { t < 0.0 });
                }
            }
        }
    }

    // The following test verifies equivalence of all of the robust orient2d variants (including
    // exact and slow variants).
    #[test]
    fn orient2d_regression_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);

        let n = 99999;
        for _ in 0..n {
            let pa = [tol(&mut rng), tol(&mut rng)];
            let pb = [12.0, 12.0];
            let pc = [24.0, 24.0];

            let o2d = predicates::orient2d;
            let o2de = predicates::orient2d_exact;
            let o2ds = predicates::orient2d_slow;

            // Test all permutations.

            for [a, b, c] in QuickPerm3::new([pa, pb, pc]) {
                assert_eq!(
                    o2d(a, b, c).partial_cmp(&0.0),
                    o2de(a, b, c).partial_cmp(&0.0),
                    "{:?}",
                    pa
                );
                assert_eq!(
                    o2d(a, b, c).partial_cmp(&0.0),
                    o2ds(a, b, c).partial_cmp(&0.0),
                    "{:?}",
                    pa
                );
            }
        }
    }

    #[test]
    fn orient2d_fast_test() {
        let a = [0.0, 1.0];
        let b = [2.0, 3.0];
        let c = [4.0, 5.0];
        assert_eq!(orient2d_fast(a, b, c), 0.0);

        // The fast orientation test should also work when the given points are sufficiently
        // non-colinear.
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let tol = 5.0e-10; // will not work with 5.0e-14

        for _ in 0..999 {
            let a = [0.5 + tol * rng.gen::<f64>(), 0.5 + tol * rng.gen::<f64>()];
            let b = [12.0, 12.0];
            let c = [24.0, 24.0];
            assert_eq!(orient2d_fast(a, b, c) > 0.0, orient2d_fast(b, c, a) > 0.0);
            assert_eq!(orient2d_fast(b, c, a) > 0.0, orient2d_fast(c, a, b) > 0.0);
            assert_eq!(orient2d_fast(a, b, c) > 0.0, orient2d_fast(b, a, c) < 0.0);
            assert_eq!(orient2d_fast(a, b, c) > 0.0, orient2d_fast(a, c, b) < 0.0);
        }
    }

    #[test]
    fn orient3d_test() {
        let a = [0.0, 1.0, 6.0];
        let b = [2.0, 3.0, 4.0];
        let c = [4.0, 5.0, 1.0];
        let d = [6.0, 2.0, 5.3];
        assert_eq!(orient3d(a, b, c, d), 10.0);
    }

    #[cfg(feature = "transpiled")]
    #[test]
    fn orient3d_transpiled_regression_test() {
        unsafe {
            transpiled::exactinit();
        }
        let mut rng: StdRng = SeedableRng::from_seed(SEED);

        let n = 9999;
        for _ in 0..n {
            let pa = [tol(&mut rng), tol(&mut rng), tol(&mut rng)];
            let pb = [12.0, 12.0, 12.0];
            let pc = [24.0, 24.0, 24.0];
            let pd = [48.0, 48.0, 48.0];

            let o3d = predicates::orient3d;
            let o3d_transpiled = transpiled::orient3d;

            // Test all permutations.
            for [a, b, c, d] in QuickPerm4::new([pa, pb, pc, pd]) {
                assert_eq!(o3d(a, c, d, b), o3d_transpiled(a, c, d, b), "{:?}", pa);
            }
        }
    }

    // The following test verifies equivalence of all of the robust orient3d variants.
    #[test]
    fn orient3d_regression_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);

        let n = 5000;
        for _ in 0..n {
            let pa = [tol(&mut rng), tol(&mut rng), tol(&mut rng)];
            let pb = [12.0, 12.0, 12.0];
            let pc = [24.0, 24.0, 24.0];
            let pd = [48.0, 48.0, 48.0];

            let o3d = predicates::orient3d;
            let o3de = predicates::orient3d_exact;
            let o3ds = predicates::orient3d_slow;

            // Test all permutations.
            // Actually these don't need to be exactly equal, they just need to compare equally to
            // 0.0. It just so happens that they are also equal.
            for [a, b, c, d] in QuickPerm4::new([pa, pb, pc, pd]) {
                assert_eq!(o3d(a, c, d, b), o3de(a, c, d, b), "{:?}", pa);
                assert_eq!(o3d(a, c, d, b), o3ds(a, c, d, b), "{:?}", pa);
            }
        }
    }

    #[test]
    fn orient3d_robustness_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);

        let n = 999;

        for o3d in &[orient3d, orient3d_exact, orient3d_slow] {
            for _ in 0..n {
                let pa = [tol(&mut rng), tol(&mut rng), tol(&mut rng)];
                let pb = [12.0, 12.0, 12.0];
                let pc = [24.0, 24.0, 24.0];
                let pd = [48.0, 48.0, 48.0];

                // Test all permutations.

                let main = o3d(pa, pb, pc, pd);

                if main == 0.0 {
                    for [a, b, c, d] in QuickPerm4::new([pa, pb, pc, pd]) {
                        assert_eq!(o3d(a, b, c, d), 0.0);
                    }
                }

                let pred = main > 0.0;
                for (i, [a, b, c, d]) in QuickPerm4::new([pa, pb, pc, pd]).enumerate() {
                    let t = o3d(a, b, c, d);
                    assert_eq!(
                        pred,
                        if i % 2 == 0 { t > 0.0 } else { t < 0.0 },
                        "{} vs. {} at {:?}",
                        t,
                        -main,
                        pa
                    );
                }
            }
        }
    }

    #[test]
    fn orient3d_fast_test() {
        let a = [0.0, 1.0, 6.0];
        let b = [2.0, 3.0, 4.0];
        let c = [4.0, 5.0, 1.0];
        let d = [6.0, 2.0, 5.3];
        assert!(f64::abs(orient3d_fast(a, b, c, d) - 10.0) < 1.0e-4);

        // The fast orientation test should also work when the given points are sufficiently
        // non-colinear.
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let tol = 5.0;

        let a = [
            0.5 + tol * rng.gen::<f64>(),
            0.5 + tol * rng.gen::<f64>(),
            0.5 + tol * rng.gen::<f64>(),
        ];
        let b = [12.0, 12.0, 12.0];
        let c = [24.0, 24.0, 24.0];
        let d = [48.0, 48.0, 48.0];
        assert_eq!(
            orient3d_fast(a, b, c, d) > 0.0,
            orient3d_fast(b, c, a, d) > 0.0
        );
        assert_eq!(
            orient3d_fast(b, a, c, d) < 0.0,
            orient3d_fast(c, b, d, a) < 0.0
        );
        // The following orientations are expected to be inconsistent
        // (this is why we need exact predicates)
        assert_ne!(
            orient3d_fast(b, c, d, a) > 0.0,
            orient3d_fast(c, d, a, b) > 0.0
        );
        assert_ne!(
            orient3d_fast(b, d, c, a) < 0.0,
            orient3d_fast(c, d, b, a) < 0.0
        );
    }

    #[test]
    fn incircle_test() {
        let a = [0.0, 1.0];
        let b = [1.0, 0.0];
        let c = [1.0, 1.0];
        let d = [0.0, 0.0];
        assert_eq!(incircle(a, b, c, d), 0.0);
        let d = [0.1, 0.1];
        assert!(incircle(a, b, c, d) > 0.0);
        let d = [-0.1, -0.1];
        assert!(incircle(a, b, c, d) < 0.0);
    }

    #[test]
    fn incircle_robustness_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);

        let n = 999;
        for ic in &[incircle, incircle_exact, incircle_slow] {
            for _ in 0..n {
                let pa = [0.0, 1.0];
                let pb = [1.0, 0.0];
                let pc = [1.0, 1.0];
                let pd = [tol(&mut rng), tol(&mut rng)];

                let main = ic(pa, pb, pc, pd);

                if main == 0.0 {
                    for [a, b, c, d] in QuickPerm4::new([pa, pb, pc, pd]) {
                        assert_eq!(ic(a, b, c, d), 0.0);
                    }
                }

                let pred = main > 0.0;
                for (i, [a, b, c, d]) in QuickPerm4::new([pa, pb, pc, pd]).enumerate() {
                    let t = ic(a, b, c, d);
                    assert_eq!(
                        pred,
                        if i % 2 == 0 { t > 0.0 } else { t < 0.0 },
                        "{} vs. {} at {:?}",
                        t,
                        -main,
                        pd
                    );
                }
            }
        }
    }

    #[test]
    fn insphere_test() {
        let a = [0.0, 1.0, 0.0];
        let b = [1.0, 0.0, 0.0];
        let c = [0.0, 0.0, 1.0];
        let d = [1.0, 1.0, 1.0];
        let e = [0.0, 0.0, 0.0];
        assert_eq!(insphere(a, b, c, d, e), 0.0);
        let e = [0.1, 0.1, 0.1];
        assert!(insphere(a, b, c, d, e) > 0.0);
        let e = [-0.1, -0.1, -0.1];
        assert!(insphere(a, b, c, d, e) < 0.0);
    }

    #[test]
    fn insphere_robustness_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);

        let n = 99;
        for is in &[insphere, insphere_exact, insphere_slow] {
            for _ in 0..n {
                let pa = [0.0, 1.0, 0.0];
                let pb = [1.0, 0.0, 0.0];
                let pc = [0.0, 0.0, 1.0];
                let pd = [1.0, 1.0, 1.0];
                let pe = [tol(&mut rng), tol(&mut rng), tol(&mut rng)];

                let main = is(pa, pb, pc, pd, pe);

                if main == 0.0 {
                    for [a, b, c, d, e] in QuickPerm5::new([pa, pb, pc, pd, pe]) {
                        assert_eq!(is(a, b, c, d, e), 0.0);
                    }
                }

                let pred = main > 0.0;
                for (i, [a, b, c, d, e]) in QuickPerm5::new([pa, pb, pc, pd, pe]).enumerate() {
                    let t = is(a, b, c, d, e);
                    assert_eq!(
                        pred,
                        if i % 2 == 0 { t > 0.0 } else { t < 0.0 },
                        "{} vs. {} at {:?}",
                        t,
                        -main,
                        pe
                    );
                }
            }
        }
    }
}
