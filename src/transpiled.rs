#![allow(dead_code)]
mod c2rust;

pub use self::c2rust::exactinit;

/*
 * Below are safe wrappers for the directly transpiled predicates located in the c2rust module.
 *
 * This module is used for testing the ported rust predicates, and is not intended to be used directly.
 */

#[inline]
pub fn orient2d(mut pa: [f64; 2], mut pb: [f64; 2], mut pc: [f64; 2]) -> f64 {
    unsafe { c2rust::orient2d(pa.as_mut_ptr(), pb.as_mut_ptr(), pc.as_mut_ptr()) }
}
#[inline]
pub fn orient2dexact(mut pa: [f64; 2], mut pb: [f64; 2], mut pc: [f64; 2]) -> f64 {
    unsafe { c2rust::orient2dexact(pa.as_mut_ptr(), pb.as_mut_ptr(), pc.as_mut_ptr()) }
}
#[inline]
pub fn orient2dslow(mut pa: [f64; 2], mut pb: [f64; 2], mut pc: [f64; 2]) -> f64 {
    unsafe { c2rust::orient2dslow(pa.as_mut_ptr(), pb.as_mut_ptr(), pc.as_mut_ptr()) }
}
#[inline]
pub fn orient3d(mut pa: [f64; 3], mut pb: [f64; 3], mut pc: [f64; 3], mut pd: [f64; 3]) -> f64 {
    unsafe {
        c2rust::orient3d(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
        )
    }
}
#[inline]
pub fn orient3dexact(
    mut pa: [f64; 3],
    mut pb: [f64; 3],
    mut pc: [f64; 3],
    mut pd: [f64; 3],
) -> f64 {
    unsafe {
        c2rust::orient3dexact(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
        )
    }
}
#[inline]
pub fn orient3dslow(mut pa: [f64; 3], mut pb: [f64; 3], mut pc: [f64; 3], mut pd: [f64; 3]) -> f64 {
    unsafe {
        c2rust::orient3dslow(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
        )
    }
}
#[inline]
pub fn incircle(mut pa: [f64; 2], mut pb: [f64; 2], mut pc: [f64; 2], mut pd: [f64; 2]) -> f64 {
    unsafe {
        c2rust::incircle(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
        )
    }
}
#[inline]
pub fn incircleexact(
    mut pa: [f64; 2],
    mut pb: [f64; 2],
    mut pc: [f64; 2],
    mut pd: [f64; 2],
) -> f64 {
    unsafe {
        c2rust::incircleexact(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
        )
    }
}
#[inline]
pub fn incircleslow(mut pa: [f64; 2], mut pb: [f64; 2], mut pc: [f64; 2], mut pd: [f64; 2]) -> f64 {
    unsafe {
        c2rust::incircleslow(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
        )
    }
}
#[inline]
pub fn insphere(
    mut pa: [f64; 3],
    mut pb: [f64; 3],
    mut pc: [f64; 3],
    mut pd: [f64; 3],
    mut pe: [f64; 3],
) -> f64 {
    unsafe {
        c2rust::insphere(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
            pe.as_mut_ptr(),
        )
    }
}
#[inline]
pub fn insphereexact(
    mut pa: [f64; 3],
    mut pb: [f64; 3],
    mut pc: [f64; 3],
    mut pd: [f64; 3],
    mut pe: [f64; 3],
) -> f64 {
    unsafe {
        c2rust::insphereexact(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
            pe.as_mut_ptr(),
        )
    }
}
#[inline]
pub fn insphereslow(
    mut pa: [f64; 3],
    mut pb: [f64; 3],
    mut pc: [f64; 3],
    mut pd: [f64; 3],
    mut pe: [f64; 3],
) -> f64 {
    unsafe {
        c2rust::insphereslow(
            pa.as_mut_ptr(),
            pb.as_mut_ptr(),
            pc.as_mut_ptr(),
            pd.as_mut_ptr(),
            pe.as_mut_ptr(),
        )
    }
}
