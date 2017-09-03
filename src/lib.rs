mod predicates;

pub struct GeometryPredicates(predicates::RawGeometryPredicates);

impl GeometryPredicates {

    #[inline(always)]
    pub fn exactinit() -> GeometryPredicates {
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn orient2d_test() {
        let pred = GeometryPredicates::exactinit();
        let a = [0.0, 1.0];
        let b = [2.0, 3.0];
        let c = [4.0, 5.0];
        assert_eq!(pred.orient2d(a, b, c), 0.0);
    }
    fn orient3d_test() {
        let pred = GeometryPredicates::exactinit();
        let a = [0.0, 1.0, 6.0];
        let b = [2.0, 3.0, 4.0];
        let c = [4.0, 5.0, 1.0];
        assert_eq!(pred.orient2d(a, b, c), 0.0);
    }
}
