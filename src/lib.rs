mod predicates;

pub struct GeometryPredicates(predicates::RawGeometryPredicates);

impl GeometryPredicates {

    #[inline(always)]
    pub fn new() -> GeometryPredicates {
        GeometryPredicates(predicates::RawGeometryPredicates::exactinit())
    }

    #[inline(always)]
    pub fn orient2d(&self, pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> f64 {
        unsafe {
            self.0.orient2d(pa.as_ptr(), pb.as_ptr(), pc.as_ptr())
        }
    }

    #[inline(always)]
    pub fn orient3d(&self, pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> f64 {
        unsafe {
            self.0.orient3d(pa.as_ptr(), pb.as_ptr(), pc.as_ptr(), pd.as_ptr())
        }
    }
    
}

/// Inexact but fast oreint2d predicate. Equivalent to a floating point determinant computation.
#[inline(always)]
pub fn orient2d_fast(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> f64 {
    unsafe {
        predicates::orient2dfast(pa.as_ptr(), pb.as_ptr(), pc.as_ptr())
    }
}

/// Inexact but fast orient3d predicate. Equivalent to a floating point determinant computation.
#[inline(always)]
pub fn orient3d_fast(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> f64 {
    unsafe {
        predicates::orient3dfast(pa.as_ptr(), pb.as_ptr(), pc.as_ptr(), pd.as_ptr())
    }
}

#[cfg(test)]
mod tests {
    extern crate rand;
    use super::*;
    use self::rand::{SeedableRng, StdRng, Rng};

    static SEED: &[usize] = &[1,2,3,4];

    #[test]
    fn orient2d_test() {
        let pred = GeometryPredicates::new();
        let a = [0.0, 1.0];
        let b = [2.0, 3.0];
        let c = [4.0, 5.0];
        assert_eq!(pred.orient2d(a, b, c), 0.0);
    }

    #[test]
    fn orient2d_robustness_test() {
        let pred = GeometryPredicates::new();
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let tol = 5.0e-14;

        for _ in 0..99999 {
            let a = [0.5 + tol*rng.gen::<f64>(), 0.5 + tol*rng.gen::<f64>()];
            let b = [12.0, 12.0];
            let c = [24.0 , 24.0];
            assert_eq!(pred.orient2d(a,b,c) > 0.0, pred.orient2d(b,c,a) > 0.0);
            assert_eq!(pred.orient2d(b,c,a) > 0.0, pred.orient2d(c,a,b) > 0.0);
            assert_eq!(pred.orient2d(a,b,c) > 0.0, pred.orient2d(b,a,c) < 0.0);
            assert_eq!(pred.orient2d(a,b,c) > 0.0, pred.orient2d(a,c,b) < 0.0);
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
            let a = [0.5 + tol*rng.gen::<f64>(), 0.5 + tol*rng.gen::<f64>()];
            let b = [12.0, 12.0];
            let c = [24.0 , 24.0];
            assert_eq!(orient2d_fast(a,b,c) > 0.0, orient2d_fast(b,c,a) > 0.0);
            assert_eq!(orient2d_fast(b,c,a) > 0.0, orient2d_fast(c,a,b) > 0.0);
            assert_eq!(orient2d_fast(a,b,c) > 0.0, orient2d_fast(b,a,c) < 0.0);
            assert_eq!(orient2d_fast(a,b,c) > 0.0, orient2d_fast(a,c,b) < 0.0);
        }
    }

    #[test]
    fn orient3d_test() {
        let pred = GeometryPredicates::new();
        let a = [0.0, 1.0, 6.0];
        let b = [2.0, 3.0, 4.0];
        let c = [4.0, 5.0, 1.0];
        let d = [6.0, 2.0, 5.3];
        assert_eq!(pred.orient3d(a, b, c, d), 10.0);
    }

    #[test]
    fn orient3d_robustness_test() {
        let pred = GeometryPredicates::new();
        let mut rng: StdRng = SeedableRng::from_seed(SEED);
        let tol = 5.0e-14;

        for _ in 0..99999 {
            let a = [0.5 + tol*rng.gen::<f64>(),
                     0.5 + tol*rng.gen::<f64>(),
                     0.5 + tol*rng.gen::<f64>()];
            let b = [12.0, 12.0, 12.0];
            let c = [24.0 , 24.0, 24.0];
            let d = [48.0 , 48.0, 48.0];
            assert_eq!(pred.orient3d(a,b,c,d) > 0.0, pred.orient3d(b,c,d,a) > 0.0);
            assert_eq!(pred.orient3d(b,c,d,a) > 0.0, pred.orient3d(c,d,a,b) > 0.0);
            assert_eq!(pred.orient3d(b,a,c,d) < 0.0, pred.orient3d(c,b,d,a) < 0.0);
            assert_eq!(pred.orient3d(b,d,c,a) < 0.0, pred.orient3d(c,d,b,a) < 0.0);
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

        let a = [0.5 + tol*rng.gen::<f64>(),
                 0.5 + tol*rng.gen::<f64>(),
                 0.5 + tol*rng.gen::<f64>()];
        let b = [12.0, 12.0, 12.0];
        let c = [24.0 , 24.0, 24.0];
        let d = [48.0 , 48.0, 48.0];
        assert_eq!(orient3d_fast(a,b,c,d) > 0.0, orient3d_fast(b,c,d,a) > 0.0);
        assert_eq!(orient3d_fast(b,a,c,d) < 0.0, orient3d_fast(c,b,d,a) < 0.0);
        // The following orientations are expected to be inconsistent
        // (this is why we need exact predicateexact predicatess)
        assert_ne!(orient3d_fast(b,c,d,a) > 0.0, orient3d_fast(c,d,a,b) > 0.0);
        assert_ne!(orient3d_fast(b,d,c,a) < 0.0, orient3d_fast(c,d,b,a) < 0.0);
    }
}
