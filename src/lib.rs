#![no_std]
mod predicates;

/**
 * Adaptive exact 2D orientation test.  Robust.
 *
 * Return a positive value if the points `pa`, `pb`, and `pc` occur
 * in counterclockwise order; a negative value if they occur
 * in clockwise order; and zero if they are collinear.  The
 * result is also a rough approximation of twice the signed
 * area of the triangle defined by the three points.
 *
 * The result returned is the determinant of a matrix.
 * This determinant is computed adaptively, in the sense that exact
 * arithmetic is used only to the degree it is needed to ensure that the
 * returned value has the correct sign.  Hence, `orient2d()` is usually quite
 * fast, but will run more slowly when the input points are collinear or
 * nearly so.
 */
#[inline]
pub fn orient2d(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> f64 {
    unsafe { predicates::orient2d(pa.as_ptr(), pb.as_ptr(), pc.as_ptr()) }
}

/**
 * Adaptive exact 3D orientation test. Robust.
 *
 * Return a positive value if the point `pd` lies below the
 * plane passing through `pa`, `pb`, and `pc`; "below" is defined so
 * that `pa`, `pb`, and `pc` appear in counterclockwise order when
 * viewed from above the plane.  Returns a negative value if
 * pd lies above the plane.  Returns zero if the points are
 * coplanar.  The result is also a rough approximation of six
 * times the signed volume of the tetrahedron defined by the
 * four points.
 *
 * The result returned is the determinant of a matrix.
 * This determinant is computed adaptively, in the sense that exact
 * arithmetic is used only to the degree it is needed to ensure that the
 * returned value has the correct sign.  Hence, orient3d() is usually quite
 * fast, but will run more slowly when the input points are coplanar or
 * nearly so.
 */
#[inline]
pub fn orient3d(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> f64 {
    unsafe { predicates::orient3d(pa.as_ptr(), pb.as_ptr(), pc.as_ptr(), pd.as_ptr()) }
}

/**
 * Adaptive exaact 2D incircle test. Robust.
 *
 * Return a positive value if the point `pd` lies inside the
 * circle passing through `pa`, `pb`, and `pc`; a negative value if
 * it lies outside; and zero if the four points are cocircular.
 * The points `pa`, `pb`, and `pc` must be in counterclockwise
 * order, or the sign of the result will be reversed.
 *
 * The result returned is the determinant of a matrix.
 * This determinant is computed adaptively, in the sense that exact
 * arithmetic is used only to the degree it is needed to ensure that the
 * returned value has the correct sign.  Hence, `incircle()` is usually quite
 * fast, but will run more slowly when the input points are cocircular or
 * nearly so.
 */
#[inline]
pub fn incircle(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> f64 {
    unsafe { predicates::incircle(pa.as_ptr(), pb.as_ptr(), pc.as_ptr(), pd.as_ptr()) }
}

/**
 * Adaptive exact 3D insphere test. Robust.
 *
 * Return a positive value if the point `pe` lies inside the
 * sphere passing through `pa`, `pb`, `pc`, and `pd`; a negative value
 * if it lies outside; and zero if the five points are
 * cospherical.  The points `pa`, `pb`, `pc`, and `pd` must be ordered
 * so that they have a positive orientation (as defined by
 * `orient3d()`), or the sign of the result will be reversed.
 *
 * The result returned is the determinant of a matrix.
 * this determinant is computed adaptively, in the sense that exact
 * arithmetic is used only to the degree it is needed to ensure that the
 * returned value has the correct sign.  Hence, `insphere()` is usually quite
 * fast, but will run more slowly when the input points are cospherical or
 * nearly so.
 */
#[inline]
pub fn insphere(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3], pe: [f64; 3]) -> f64 {
    unsafe {
        predicates::insphere(
            pa.as_ptr(),
            pb.as_ptr(),
            pc.as_ptr(),
            pd.as_ptr(),
            pe.as_ptr(),
        )
    }
}

/// Approximate 2D orientation test. Nonrobust version of `orient2d`.
#[inline]
pub fn orient2d_fast(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> f64 {
    unsafe { predicates::orient2dfast(pa.as_ptr(), pb.as_ptr(), pc.as_ptr()) }
}

/// Approximate 3D orientation test. Nonrobust version of `orient3d`.
#[inline]
pub fn orient3d_fast(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> f64 {
    unsafe { predicates::orient3dfast(pa.as_ptr(), pb.as_ptr(), pc.as_ptr(), pd.as_ptr()) }
}

/// Approximate 2D incircle test. Nonrobust version of `incircle`
#[inline]
pub fn incircle_fast(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> f64 {
    unsafe { predicates::incirclefast(pa.as_ptr(), pb.as_ptr(), pc.as_ptr(), pd.as_ptr()) }
}

/// Approximate 3D insphere test. Nonrobust version of `insphere`
#[inline]
pub fn insphere_fast(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3], pe: [f64; 3]) -> f64 {
    unsafe {
        predicates::inspherefast(
            pa.as_ptr(),
            pb.as_ptr(),
            pc.as_ptr(),
            pd.as_ptr(),
            pe.as_ptr(),
        )
    }
}

#[cfg(test)]
mod tests {
    extern crate rand;
    use self::rand::{Rng, SeedableRng, StdRng};
    use super::*;

    const SEED: &[usize] = &[1, 2, 3, 4];

    #[test]
    fn orient2d_test() {
        let a = [0.0, 1.0];
        let b = [2.0, 3.0];
        let c = [4.0, 5.0];
        assert_eq!(orient2d(a, b, c), 0.0);
    }

    #[test]
    fn orient2d_robustness_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let tol = 5.0e-14;

        #[cfg(miri)]
        let n = 99;
        #[cfg(not(miri))]
        let n = 9999999;
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

    #[test]
    fn orient3d_robustness_test() {
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let tol = 5.0e-14;

        #[cfg(miri)]
        let n = 99;
        #[cfg(not(miri))]
        let n = 99999;
        for _ in 0..n {
            let a = [
                0.5 + tol * rng.gen::<f64>(),
                0.5 + tol * rng.gen::<f64>(),
                0.5 + tol * rng.gen::<f64>(),
            ];
            let b = [12.0, 12.0, 12.0];
            let c = [24.0, 24.0, 24.0];
            let d = [48.0, 48.0, 48.0];
            assert_eq!(orient3d(a, b, c, d) > 0.0, orient3d(b, c, d, a) > 0.0);
            assert_eq!(orient3d(b, c, d, a) > 0.0, orient3d(c, d, a, b) > 0.0);
            assert_eq!(orient3d(b, a, c, d) < 0.0, orient3d(c, b, d, a) < 0.0);
            assert_eq!(orient3d(b, d, c, a) < 0.0, orient3d(c, d, b, a) < 0.0);
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
            orient3d_fast(b, c, d, a) > 0.0
        );
        assert_eq!(
            orient3d_fast(b, a, c, d) < 0.0,
            orient3d_fast(c, b, d, a) < 0.0
        );
        // The following orientations are expected to be inconsistent
        // (this is why we need exact predicateexact predicatess)
        assert_ne!(
            orient3d_fast(b, c, d, a) > 0.0,
            orient3d_fast(c, d, a, b) > 0.0
        );
        assert_ne!(
            orient3d_fast(b, d, c, a) < 0.0,
            orient3d_fast(c, d, b, a) < 0.0
        );
    }
}
