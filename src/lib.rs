#![no_std]
mod predicates;

#[cfg(feature = "transpiled")]
mod transpiled;

pub use predicates::{
    orient2d,
    orient2d_fast,
    orient2d_exact,
    orient2d_slow,
    orient3d,
    orient3d_fast,
    orient3d_exact,
    orient3d_slow,
    incircle,
    incircle_fast,
    incircle_exact,
    incircle_slow,
    insphere,
    insphere_fast,
    insphere_exact,
    insphere_slow,
};

#[cfg(test)]
mod tests {
    extern crate rand;
    use self::rand::{Rng, SeedableRng, StdRng};
    use super::*;

    const SEED: &[usize] = &[1, 2, 3, 4];

    /* Note on robustness testing
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

    #[test]
    fn orient2d_test() {
        let a = [0.0, 1.0];
        let b = [2.0, 3.0];
        let c = [4.0, 5.0];
        assert_eq!(orient2d(a, b, c), 0.0);
    }

    #[cfg(feature = "transpiled")]
    #[test]
    fn orient2d_regression_test() {
        unsafe { transpiled::exactinit(); }
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let tol = 5.0e-14;

        let n = 99999;
        for _ in 0..n {
            let a = [0.5 + tol * rng.gen::<f64>(), 0.5 + tol * rng.gen::<f64>()];
            let b = [12.0, 12.0];
            let c = [24.0, 24.0];
            let o2d = orient2d;
            let o2d_transpiled = transpiled::orient2d;
            assert_eq!(o2d(a, b, c), o2d_transpiled(a, b, c), "{:?}", a);
            assert_eq!(o2d(b, c, a), o2d_transpiled(b, c, a), "{:?}", a);
            assert_eq!(o2d(c, a, b), o2d_transpiled(c, a, b), "{:?}", a);
            assert_eq!(o2d(a, c, b), o2d_transpiled(a, c, b), "{:?}", a);
            assert_eq!(o2d(c, b, a), o2d_transpiled(c, b, a), "{:?}", a);
            assert_eq!(o2d(b, a, c), o2d_transpiled(b, a, c), "{:?}", a);
        }
    }

    #[test]
    fn orient2d_robustness_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let tol = 5.0e-14;

        let n = 99999;
        for _ in 0..n {
            let a = [0.5 + tol * rng.gen::<f64>(), 0.5 + tol * rng.gen::<f64>()];
            let b = [12.0, 12.0];
            let c = [24.0, 24.0];
            assert_eq!(orient2d(a, b, c) > 0.0, orient2d(b, c, a) > 0.0);
            assert_eq!(orient2d(b, c, a) > 0.0, orient2d(c, a, b) > 0.0);
            assert_eq!(orient2d(a, b, c) > 0.0, orient2d(b, a, c) < 0.0);
            assert_eq!(orient2d(a, b, c) > 0.0, orient2d(a, c, b) < 0.0);
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
        unsafe { transpiled::exactinit(); }
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let max_exp = (EXP_BOUNDS[0] + 1) as f64;

        // Generate a small non-negative number.
        let mut tol = || ((max_exp*rng.gen::<f64>()).round()).exp() * (rng.gen::<f64>() - 0.5);

        let n = 9999;
        for _ in 0..n {
            let a = [tol(), tol(), tol()];
            let b = [12.0, 12.0, 12.0];
            let c = [24.0, 24.0, 24.0];
            let d = [48.0, 48.0, 48.0];

            let o3d = predicates::orient3d;
            let o3d_transpiled = transpiled::orient3d;

            assert_eq!(o3d(a, c, d, b), o3d_transpiled(a, c, d, b), "{:?}", a);
            assert_eq!(o3d(a, d, b, c), o3d_transpiled(a, d, b, c), "{:?}", a);
            assert_eq!(o3d(b, a, d, c), o3d_transpiled(b, a, d, c), "{:?}", a);
            assert_eq!(o3d(b, c, a, d), o3d_transpiled(b, c, a, d), "{:?}", a);
            assert_eq!(o3d(b, d, c, a), o3d_transpiled(b, d, c, a), "{:?}", a);
            assert_eq!(o3d(c, a, b, d), o3d_transpiled(c, a, b, d), "{:?}", a);
            assert_eq!(o3d(c, b, d, a), o3d_transpiled(c, b, d, a), "{:?}", a);
            assert_eq!(o3d(c, d, a, b), o3d_transpiled(c, d, a, b), "{:?}", a);
            assert_eq!(o3d(d, a, c, b), o3d_transpiled(d, a, c, b), "{:?}", a);
            assert_eq!(o3d(d, b, a, c), o3d_transpiled(d, b, a, c), "{:?}", a);
            assert_eq!(o3d(d, c, b, a), o3d_transpiled(d, c, b, a), "{:?}", a);
            assert_eq!(o3d(a, b, d, c), o3d_transpiled(a, b, d, c), "{:?}", a);
            assert_eq!(o3d(a, c, b, d), o3d_transpiled(a, c, b, d), "{:?}", a);
            assert_eq!(o3d(a, d, c, b), o3d_transpiled(a, d, c, b), "{:?}", a);
            assert_eq!(o3d(b, a, c, d), o3d_transpiled(b, a, c, d), "{:?}", a);
            assert_eq!(o3d(b, c, d, a), o3d_transpiled(b, c, d, a), "{:?}", a);
            assert_eq!(o3d(b, d, a, c), o3d_transpiled(b, d, a, c), "{:?}", a);
            assert_eq!(o3d(c, a, d, b), o3d_transpiled(c, a, d, b), "{:?}", a);
            assert_eq!(o3d(c, b, a, d), o3d_transpiled(c, b, a, d), "{:?}", a);
            assert_eq!(o3d(c, d, b, a), o3d_transpiled(c, d, b, a), "{:?}", a);
            assert_eq!(o3d(d, a, b, c), o3d_transpiled(d, a, b, c), "{:?}", a);
            assert_eq!(o3d(d, b, c, a), o3d_transpiled(d, b, c, a), "{:?}", a);
            assert_eq!(o3d(d, c, a, b), o3d_transpiled(d, c, a, b), "{:?}", a);
        }
    }

    // The following test verifies equivalence of all of the robust orient3d variants.
    #[test]
    fn orient3d_regression_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);

        let max_exp = (EXP_BOUNDS[0] + 1) as f64;

        // Generate a small non-negative number.
        let mut tol = || 10.0_f64.powi((max_exp*rng.gen::<f64>()).round() as i32) * (rng.gen::<f64>() - 0.5);

        let n = 9999;
        for _ in 0..n {
            let a = [tol(), tol(), tol()];
            let b = [12.0, 12.0, 12.0];
            let c = [24.0, 24.0, 24.0];
            let d = [48.0, 48.0, 48.0];

            let o3d = predicates::orient3d;
            let o3de = predicates::orient3d_exact;
            let o3ds = predicates::orient3d_slow;

            assert_eq!(o3d(a, c, d, b), o3de(a, c, d, b), "{:?}", a);
            assert_eq!(o3d(a, d, b, c), o3de(a, d, b, c), "{:?}", a);
            assert_eq!(o3d(b, a, d, c), o3de(b, a, d, c), "{:?}", a);
            assert_eq!(o3d(b, c, a, d), o3de(b, c, a, d), "{:?}", a);
            assert_eq!(o3d(b, d, c, a), o3de(b, d, c, a), "{:?}", a);
            assert_eq!(o3d(c, a, b, d), o3de(c, a, b, d), "{:?}", a);
            assert_eq!(o3d(c, b, d, a), o3de(c, b, d, a), "{:?}", a);
            assert_eq!(o3d(c, d, a, b), o3de(c, d, a, b), "{:?}", a);
            assert_eq!(o3d(d, a, c, b), o3de(d, a, c, b), "{:?}", a);
            assert_eq!(o3d(d, b, a, c), o3de(d, b, a, c), "{:?}", a);
            assert_eq!(o3d(d, c, b, a), o3de(d, c, b, a), "{:?}", a);
            assert_eq!(o3d(a, b, d, c), o3de(a, b, d, c), "{:?}", a);
            assert_eq!(o3d(a, c, b, d), o3de(a, c, b, d), "{:?}", a);
            assert_eq!(o3d(a, d, c, b), o3de(a, d, c, b), "{:?}", a);
            assert_eq!(o3d(b, a, c, d), o3de(b, a, c, d), "{:?}", a);
            assert_eq!(o3d(b, c, d, a), o3de(b, c, d, a), "{:?}", a);
            assert_eq!(o3d(b, d, a, c), o3de(b, d, a, c), "{:?}", a);
            assert_eq!(o3d(c, a, d, b), o3de(c, a, d, b), "{:?}", a);
            assert_eq!(o3d(c, b, a, d), o3de(c, b, a, d), "{:?}", a);
            assert_eq!(o3d(c, d, b, a), o3de(c, d, b, a), "{:?}", a);
            assert_eq!(o3d(d, a, b, c), o3de(d, a, b, c), "{:?}", a);
            assert_eq!(o3d(d, b, c, a), o3de(d, b, c, a), "{:?}", a);
            assert_eq!(o3d(d, c, a, b), o3de(d, c, a, b), "{:?}", a);

            assert_eq!(o3d(a, c, d, b), o3ds(a, c, d, b), "{:?}", a);
            assert_eq!(o3d(a, d, b, c), o3ds(a, d, b, c), "{:?}", a);
            assert_eq!(o3d(b, a, d, c), o3ds(b, a, d, c), "{:?}", a);
            assert_eq!(o3d(b, c, a, d), o3ds(b, c, a, d), "{:?}", a);
            assert_eq!(o3d(b, d, c, a), o3ds(b, d, c, a), "{:?}", a);
            assert_eq!(o3d(c, a, b, d), o3ds(c, a, b, d), "{:?}", a);
            assert_eq!(o3d(c, b, d, a), o3ds(c, b, d, a), "{:?}", a);
            assert_eq!(o3d(c, d, a, b), o3ds(c, d, a, b), "{:?}", a);
            assert_eq!(o3d(d, a, c, b), o3ds(d, a, c, b), "{:?}", a);
            assert_eq!(o3d(d, b, a, c), o3ds(d, b, a, c), "{:?}", a);
            assert_eq!(o3d(d, c, b, a), o3ds(d, c, b, a), "{:?}", a);
            assert_eq!(o3d(a, b, d, c), o3ds(a, b, d, c), "{:?}", a);
            assert_eq!(o3d(a, c, b, d), o3ds(a, c, b, d), "{:?}", a);
            assert_eq!(o3d(a, d, c, b), o3ds(a, d, c, b), "{:?}", a);
            assert_eq!(o3d(b, a, c, d), o3ds(b, a, c, d), "{:?}", a);
            assert_eq!(o3d(b, c, d, a), o3ds(b, c, d, a), "{:?}", a);
            assert_eq!(o3d(b, d, a, c), o3ds(b, d, a, c), "{:?}", a);
            assert_eq!(o3d(c, a, d, b), o3ds(c, a, d, b), "{:?}", a);
            assert_eq!(o3d(c, b, a, d), o3ds(c, b, a, d), "{:?}", a);
            assert_eq!(o3d(c, d, b, a), o3ds(c, d, b, a), "{:?}", a);
            assert_eq!(o3d(d, a, b, c), o3ds(d, a, b, c), "{:?}", a);
            assert_eq!(o3d(d, b, c, a), o3ds(d, b, c, a), "{:?}", a);
            assert_eq!(o3d(d, c, a, b), o3ds(d, c, a, b), "{:?}", a);
        }
    }

    #[test]
    fn orient3d_robustness_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let max_exp = (EXP_BOUNDS[0] + 1) as f64;

        // Generate a small non-negative number.
        let mut tol = || ((max_exp*rng.gen::<f64>()).round()).exp() * (rng.gen::<f64>() - 0.5);

        let n = 999;

        for o3d in &[orient3d, orient3d_exact, orient3d_slow] {
            for _ in 0..n {
                let a = [
                    tol(),
                    tol(),
                    tol(),
                ];
                let b = [12.0, 12.0, 12.0];
                let c = [24.0, 24.0, 24.0];
                let d = [48.0, 48.0, 48.0];

                let main = o3d(a, b, c, d);

                if a[0] == 0.0 && a[1] == 0.0 {
                    assert_eq!(main, 0.0);
                    assert_eq!(o3d(a, c, d, b), 0.0);
                    assert_eq!(o3d(a, d, b, c), 0.0);
                    assert_eq!(o3d(b, a, d, c), 0.0);
                    assert_eq!(o3d(b, c, a, d), 0.0);
                    assert_eq!(o3d(b, d, c, a), 0.0);
                    assert_eq!(o3d(c, a, b, d), 0.0);
                    assert_eq!(o3d(c, b, d, a), 0.0);
                    assert_eq!(o3d(c, d, a, b), 0.0);
                    assert_eq!(o3d(d, a, c, b), 0.0);
                    assert_eq!(o3d(d, b, a, c), 0.0);
                    assert_eq!(o3d(d, c, b, a), 0.0);
                    assert_eq!(o3d(a, b, d, c), 0.0);
                    assert_eq!(o3d(a, c, b, d), 0.0);
                    assert_eq!(o3d(a, d, c, b), 0.0);
                    assert_eq!(o3d(b, a, c, d), 0.0);
                    assert_eq!(o3d(b, c, d, a), 0.0);
                    assert_eq!(o3d(b, d, a, c), 0.0);
                    assert_eq!(o3d(c, a, d, b), 0.0);
                    assert_eq!(o3d(c, b, a, d), 0.0);
                    assert_eq!(o3d(c, d, b, a), 0.0);
                    assert_eq!(o3d(d, a, b, c), 0.0);
                    assert_eq!(o3d(d, b, c, a), 0.0);
                    assert_eq!(o3d(d, c, a, b), 0.0);
                }

                let pred = main > 0.0;
                assert_eq!(pred, o3d(a, c, d, b) > 0.0, "{} vs. {} at {:?}", o3d(a, c, d, b), main, a);
                assert_eq!(pred, o3d(a, d, b, c) > 0.0, "{} vs. {} at {:?}", o3d(a, d, b, c), main, a);
                assert_eq!(pred, o3d(b, a, d, c) > 0.0, "{} vs. {} at {:?}", o3d(b, a, d, c), main, a);
                assert_eq!(pred, o3d(b, c, a, d) > 0.0, "{} vs. {} at {:?}", o3d(b, c, a, d), main, a);
                assert_eq!(pred, o3d(b, d, c, a) > 0.0, "{} vs. {} at {:?}", o3d(b, d, c, a), main, a);
                assert_eq!(pred, o3d(c, a, b, d) > 0.0, "{} vs. {} at {:?}", o3d(c, a, b, d), main, a);
                assert_eq!(pred, o3d(c, b, d, a) > 0.0, "{} vs. {} at {:?}", o3d(c, b, d, a), main, a);
                assert_eq!(pred, o3d(c, d, a, b) > 0.0, "{} vs. {} at {:?}", o3d(c, d, a, b), main, a);
                assert_eq!(pred, o3d(d, a, c, b) > 0.0, "{} vs. {} at {:?}", o3d(d, a, c, b), main, a);
                assert_eq!(pred, o3d(d, b, a, c) > 0.0, "{} vs. {} at {:?}", o3d(d, b, a, c), main, a);
                assert_eq!(pred, o3d(d, c, b, a) > 0.0, "{} vs. {} at {:?}", o3d(d, c, b, a), main, a);

                assert_eq!(pred, o3d(a, b, d, c) < 0.0, "{} vs. {} at {:?}", o3d(a, b, d, c), -main, a);
                assert_eq!(pred, o3d(a, c, b, d) < 0.0, "{} vs. {} at {:?}", o3d(a, c, b, d), -main, a);
                assert_eq!(pred, o3d(a, d, c, b) < 0.0, "{} vs. {} at {:?}", o3d(a, d, c, b), -main, a);
                assert_eq!(pred, o3d(b, a, c, d) < 0.0, "{} vs. {} at {:?}", o3d(b, a, c, d), -main, a);
                assert_eq!(pred, o3d(b, c, d, a) < 0.0, "{} vs. {} at {:?}", o3d(b, c, d, a), -main, a);
                assert_eq!(pred, o3d(b, d, a, c) < 0.0, "{} vs. {} at {:?}", o3d(b, d, a, c), -main, a);
                assert_eq!(pred, o3d(c, a, d, b) < 0.0, "{} vs. {} at {:?}", o3d(c, a, d, b), -main, a);
                assert_eq!(pred, o3d(c, b, a, d) < 0.0, "{} vs. {} at {:?}", o3d(c, b, a, d), -main, a);
                assert_eq!(pred, o3d(c, d, b, a) < 0.0, "{} vs. {} at {:?}", o3d(c, d, b, a), -main, a);
                assert_eq!(pred, o3d(d, a, b, c) < 0.0, "{} vs. {} at {:?}", o3d(d, a, b, c), -main, a);
                assert_eq!(pred, o3d(d, b, c, a) < 0.0, "{} vs. {} at {:?}", o3d(d, b, c, a), -main, a);
                assert_eq!(pred, o3d(d, c, a, b) < 0.0, "{} vs. {} at {:?}", o3d(d, c, a, b), -main, a);
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
        let tol = 5.0; // will fail (see below) with all tolerances I tried

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
        let max_exp = (EXP_BOUNDS[0] + 1) as f64;

        // Generate a small non-negative number.
        let mut tol = || ((max_exp*rng.gen::<f64>()).round()).exp() * (rng.gen::<f64>() - 0.5);

        let n = 999;
        for ic in &[incircle, incircle_exact, incircle_slow] {
            for _ in 0..n {
                let a = [0.0, 1.0];
                let b = [1.0, 0.0];
                let c = [1.0, 1.0];
                let d = [tol(), tol()];

                let main = ic(a, b, c, d);

                if d[0] == 0.0 && d[1] == 0.0 {
                    assert_eq!(main, 0.0);
                    assert_eq!(ic(a, c, d, b), 0.0);
                    assert_eq!(ic(a, d, b, c), 0.0);
                    assert_eq!(ic(b, a, d, c), 0.0);
                    assert_eq!(ic(b, c, a, d), 0.0);
                    assert_eq!(ic(b, d, c, a), 0.0);
                    assert_eq!(ic(c, a, b, d), 0.0);
                    assert_eq!(ic(c, b, d, a), 0.0);
                    assert_eq!(ic(c, d, a, b), 0.0);
                    assert_eq!(ic(d, a, c, b), 0.0);
                    assert_eq!(ic(d, b, a, c), 0.0);
                    assert_eq!(ic(d, c, b, a), 0.0);
                    assert_eq!(ic(a, b, d, c), 0.0);
                    assert_eq!(ic(a, c, b, d), 0.0);
                    assert_eq!(ic(a, d, c, b), 0.0);
                    assert_eq!(ic(b, a, c, d), 0.0);
                    assert_eq!(ic(b, c, d, a), 0.0);
                    assert_eq!(ic(b, d, a, c), 0.0);
                    assert_eq!(ic(c, a, d, b), 0.0);
                    assert_eq!(ic(c, b, a, d), 0.0);
                    assert_eq!(ic(c, d, b, a), 0.0);
                    assert_eq!(ic(d, a, b, c), 0.0);
                    assert_eq!(ic(d, b, c, a), 0.0);
                    assert_eq!(ic(d, c, a, b), 0.0);
                }

                let pred = main > 0.0;
                assert_eq!(pred, ic(a, c, d, b) > 0.0, "{} vs. {}: {:?}", ic(a, c, d, b), main, d);
                assert_eq!(pred, ic(a, d, b, c) > 0.0, "{} vs. {}: {:?}", ic(a, d, b, c), main, d);
                assert_eq!(pred, ic(b, a, d, c) > 0.0, "{} vs. {}: {:?}", ic(b, a, d, c), main, d);
                assert_eq!(pred, ic(b, c, a, d) > 0.0, "{} vs. {}: {:?}", ic(b, c, a, d), main, d);
                assert_eq!(pred, ic(b, d, c, a) > 0.0, "{} vs. {}: {:?}", ic(b, d, c, a), main, d);
                assert_eq!(pred, ic(c, a, b, d) > 0.0, "{} vs. {}: {:?}", ic(c, a, b, d), main, d);
                assert_eq!(pred, ic(c, b, d, a) > 0.0, "{} vs. {}: {:?}", ic(c, b, d, a), main, d);
                assert_eq!(pred, ic(c, d, a, b) > 0.0, "{} vs. {}: {:?}", ic(c, d, a, b), main, d);
                assert_eq!(pred, ic(d, a, c, b) > 0.0, "{} vs. {}: {:?}", ic(d, a, c, b), main, d);
                assert_eq!(pred, ic(d, b, a, c) > 0.0, "{} vs. {}: {:?}", ic(d, b, a, c), main, d);
                assert_eq!(pred, ic(d, c, b, a) > 0.0, "{} vs. {}: {:?}", ic(d, c, b, a), main, d);

                assert_eq!(pred, ic(a, b, d, c) < 0.0, "{} vs. {}: {:?}", ic(a, b, d, c), -main, d);
                assert_eq!(pred, ic(a, c, b, d) < 0.0, "{} vs. {}: {:?}", ic(a, c, b, d), -main, d);
                assert_eq!(pred, ic(a, d, c, b) < 0.0, "{} vs. {}: {:?}", ic(a, d, c, b), -main, d);
                assert_eq!(pred, ic(b, a, c, d) < 0.0, "{} vs. {}: {:?}", ic(b, a, c, d), -main, d);
                assert_eq!(pred, ic(b, c, d, a) < 0.0, "{} vs. {}: {:?}", ic(b, c, d, a), -main, d);
                assert_eq!(pred, ic(b, d, a, c) < 0.0, "{} vs. {}: {:?}", ic(b, d, a, c), -main, d);
                assert_eq!(pred, ic(c, a, d, b) < 0.0, "{} vs. {}: {:?}", ic(c, a, d, b), -main, d);
                assert_eq!(pred, ic(c, b, a, d) < 0.0, "{} vs. {}: {:?}", ic(c, b, a, d), -main, d);
                assert_eq!(pred, ic(c, d, b, a) < 0.0, "{} vs. {}: {:?}", ic(c, d, b, a), -main, d);
                assert_eq!(pred, ic(d, a, b, c) < 0.0, "{} vs. {}: {:?}", ic(d, a, b, c), -main, d);
                assert_eq!(pred, ic(d, b, c, a) < 0.0, "{} vs. {}: {:?}", ic(d, b, c, a), -main, d);
                assert_eq!(pred, ic(d, c, a, b) < 0.0, "{} vs. {}: {:?}", ic(d, c, a, b), -main, d);
            }
        }
    }
}
