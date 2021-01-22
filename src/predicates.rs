#![allow(
    unused_variables,
    dead_code,
    mutable_transmutes,
    non_camel_case_types,
    non_snake_case,
    non_upper_case_globals,
    unused_assignments,
    unused_mut
)]

// f64::abs is not available in core. See https://github.com/rust-lang/rust/issues/50145
// This implementation is identical to abs on x86 but not on arm at the time of this writing.
// So for now we only use it when no_std is not enabled.
#[inline]
pub fn Absolute(a: f64) -> f64 {
    f64::from_bits(a.to_bits() & 0x7FFF_FFFF_FFFF_FFFF)
}

struct PredicateParams {
    splitter: f64, // = 2^ceiling(p / 2) + 1.  Used to split floats in half.
    /* A set of coefficients used to calculate maximum roundoff errors.          */
    resulterrbound: f64,
    ccwerrboundA: f64,
    ccwerrboundB: f64,
    ccwerrboundC: f64,
    o3derrboundA: f64,
    o3derrboundB: f64,
    o3derrboundC: f64,
    iccerrboundA: f64,
    iccerrboundB: f64,
    iccerrboundC: f64,
    isperrboundA: f64,
    isperrboundB: f64,
    isperrboundC: f64,
}

/* ***************************************************************************/
/*  The following are constants used in exact arithmetic                     */
/*                                                                           */
/*  `EPSILON' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-point error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-point numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  See exactinit() in predicates.c for the function used to generate these. */
/*                                                                           */
/* ***************************************************************************/
const EPSILON: f64 = 0.000_000_000_000_000_111_022_302_462_515_65;

const PARAMS: PredicateParams = PredicateParams {
    splitter: 134_217_729f64,
    resulterrbound: (3.0 + 8.0 * EPSILON) * EPSILON,
    ccwerrboundA: (3.0 + 16.0 * EPSILON) * EPSILON,
    ccwerrboundB: (2.0 + 12.0 * EPSILON) * EPSILON,
    ccwerrboundC: (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON,
    o3derrboundA: (7.0f64 + 56.0f64 * EPSILON) * EPSILON,
    o3derrboundB: (3.0f64 + 28.0f64 * EPSILON) * EPSILON,
    o3derrboundC: (26.0f64 + 288.0f64 * EPSILON) * EPSILON * EPSILON,
    iccerrboundA: (10.0 + 96.0 * EPSILON) * EPSILON,
    iccerrboundB: (4.0 + 48.0 * EPSILON) * EPSILON,
    iccerrboundC: (44.0 + 576.0 * EPSILON) * EPSILON * EPSILON,
    isperrboundA: (16.0f64 + 224.0f64 * EPSILON) * EPSILON,
    isperrboundB: (5.0f64 + 72.0f64 * EPSILON) * EPSILON,
    isperrboundC: (71.0f64 + 1408.0f64 * EPSILON) * EPSILON * EPSILON,
};

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */

#[inline]
pub fn Fast_Two_Sum_Tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt = x - a;
    b - bvirt
}

#[inline]
pub fn Fast_Two_Sum(a: f64, b: f64) -> [f64; 2] {
    let x = a + b;
    let y = Fast_Two_Sum_Tail(a, b, x);
    [x, y]
}

#[inline]
pub fn Fast_Two_Diff_Tail(a: f64, b: f64, x: f64) -> f64 {
    let mut bvirt: f64 = a - x;
    return bvirt - b;
}

#[inline]
pub unsafe fn Fast_Two_Diff(mut a: f64,
                                       mut b: f64,
                                       mut x: *mut f64,
                                       mut y: *mut f64) {
    *x = a - b;
    *y = Fast_Two_Diff_Tail(a, b, *x);
}

#[inline]
pub unsafe fn Two_Sum_Tail(mut a: f64,
                                      mut b: f64,
                                      mut x: f64)
 -> f64 {
    let mut bvirt: f64 = x - a;
    let mut avirt: f64 = x - bvirt;
    let mut bround: f64 = b - bvirt;
    let mut around: f64 = a - avirt;
    return around + bround;
}

#[inline]
pub unsafe fn Two_Sum(mut a: f64, mut b: f64,
                                 mut x: *mut f64,
                                 mut y: *mut f64) {
    *x = a + b;
    *y = Two_Sum_Tail(a, b, *x);
}

#[inline]
pub unsafe fn Two_Diff_Tail(mut a: f64,
                                       mut b: f64,
                                       mut x: f64)
 -> f64 {
    let mut bvirt: f64 = a - x;
    let mut avirt: f64 = x + bvirt;
    let mut bround: f64 = bvirt - b;
    let mut around: f64 = a - avirt;
    return around + bround;
}

#[inline]
pub unsafe fn Two_Diff(mut a: f64,
                                  mut b: f64,
                                  mut x: *mut f64,
                                  mut y: *mut f64) {
    *x = a - b;
    *y = Two_Diff_Tail(a, b, *x);
}

#[inline]
pub unsafe fn Split(mut a: f64,
                               mut ahi: *mut f64,
                               mut alo: *mut f64) {
    let mut c: f64 = PARAMS.splitter * a;
    let mut abig: f64 = c - a;
    *ahi = c - abig;
    *alo = a - *ahi;
}

#[inline]
pub unsafe fn Two_Product_Tail(mut a: f64,
                                          mut b: f64,
                                          mut x: f64)
 -> f64 {
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    Split(a, &mut ahi, &mut alo);
    Split(b, &mut bhi, &mut blo);
    let mut err1: f64 = x - ahi * bhi;
    let mut err2: f64 = err1 - alo * bhi;
    let mut err3: f64 = err2 - ahi * blo;
    return alo * blo - err3;
}

#[inline]
pub unsafe fn Two_Product(mut a: f64,
                                     mut b: f64,
                                     mut x: *mut f64,
                                     mut y: *mut f64) {
    *x = a * b;
    *y = Two_Product_Tail(a, b, *x);
}
/* Two_Product_Presplit() is Two_Product() where one of the inputs has       */
/*   already been split.  Avoids redundant splitting.                        */

#[inline]
pub unsafe fn Two_Product_Presplit(mut a: f64,
                                              mut b: f64,
                                              mut bhi: f64,
                                              mut blo: f64,
                                              mut x: *mut f64,
                                              mut y: *mut f64) {
    *x = a * b;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    Split(a, &mut ahi, &mut alo);
    let mut err1: f64 = *x - ahi * bhi;
    let mut err2: f64 = err1 - alo * bhi;
    let mut err3: f64 = err2 - ahi * blo;
    *y = alo * blo - err3;
}
/* Two_Product_2Presplit() is Two_Product() where both of the inputs have    */
/*   already been split.  Avoids redundant splitting.                        */

#[inline]
pub unsafe fn Two_Product_2Presplit(mut a: f64,
                                               mut ahi: f64,
                                               mut alo: f64,
                                               mut b: f64,
                                               mut bhi: f64,
                                               mut blo: f64,
                                               mut x: *mut f64,
                                               mut y: *mut f64) {
    *x = a * b;
    let mut err1: f64 = *x - ahi * bhi;
    let mut err2: f64 = err1 - alo * bhi;
    let mut err3: f64 = err2 - ahi * blo;
    *y = alo * blo - err3;
}
/* Square() can be done more quickly than Two_Product().                     */

#[inline]
pub unsafe fn Square_Tail(mut a: f64,
                                     mut x: f64)
 -> f64 {
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    Split(a, &mut ahi, &mut alo);
    let mut err1: f64 = x - ahi * ahi;
    let mut err3: f64 = err1 - (ahi + ahi) * alo;
    return alo * alo - err3;
}

#[inline]
pub unsafe fn Square(mut a: f64,
                                mut x: *mut f64,
                                mut y: *mut f64) {
    *x = a * a;
    *y = Square_Tail(a, *x);
}
/* Macros for summing expansions of various fixed lengths.  These are all    */
/*   unrolled versions of Expansion_Sum().                                   */

#[inline]
pub unsafe fn Two_One_Sum(mut a1: f64,
                                     mut a0: f64,
                                     mut b: f64,
                                     mut out: *mut f64) {
    let mut _i: f64 = 0.; // _k, _6, _5, _4, _3, _2, _1, _0, x0);
    Two_Sum(a0, b, &mut _i, &mut *out.offset(0 as i32 as isize));
    Two_Sum(a1, _i, &mut *out.offset(2 as i32 as isize),
            &mut *out.offset(1 as i32 as isize));
}

#[inline]
pub unsafe fn Two_One_Diff(mut a1: f64,
                                      mut a0: f64,
                                      mut b: f64,
                                      mut out: *mut f64) {
    let mut _i: f64 = 0.;
    Two_Diff(a0, b, &mut _i, &mut *out.offset(0 as i32 as isize));
    Two_Sum(a1, _i, &mut *out.offset(2 as i32 as isize),
            &mut *out.offset(1 as i32 as isize));
}

#[inline]
pub unsafe fn Two_Two_Sum(mut a1: f64,
                                     mut a0: f64,
                                     mut b1: f64,
                                     mut b0: f64,
                                     mut out: *mut f64) {
    let mut sum: [f64; 3] = [0.; 3];
    Two_One_Sum(a1, a0, b0, sum.as_mut_ptr());
    let mut _j: f64 = sum[2 as i32 as usize];
    let mut _0: f64 = sum[1 as i32 as usize];
    *out.offset(0 as i32 as isize) = sum[0 as i32 as usize];
    Two_One_Sum(_j, _0, b1, out.offset(1 as i32 as isize));
}

#[inline]
pub unsafe fn Two_Two_Diff(mut a1: f64,
                                      mut a0: f64,
                                      mut b1: f64,
                                      mut b0: f64,
                                      mut out: *mut f64) {
    let mut diff: [f64; 3] = [0.; 3];
    Two_One_Diff(a1, a0, b0, diff.as_mut_ptr());
    let mut _j: f64 = diff[2 as i32 as usize];
    let mut _0: f64 = diff[1 as i32 as usize];
    *out.offset(0 as i32 as isize) = diff[0 as i32 as usize];
    Two_One_Diff(_j, _0, b1, out.offset(1 as i32 as isize));
}

#[inline]
pub unsafe fn Four_One_Sum(mut a3: f64,
                                      mut a2: f64,
                                      mut a1: f64,
                                      mut a0: f64,
                                      mut b: f64,
                                      mut out: *mut f64) {
    let mut sum: [f64; 3] = [0.; 3];
    Two_One_Sum(a1, a0, b, sum.as_mut_ptr());
    let mut _j: f64 = sum[2 as i32 as usize];
    *out.offset(1 as i32 as isize) = sum[1 as i32 as usize];
    *out.offset(0 as i32 as isize) = sum[0 as i32 as usize];
    Two_One_Sum(a3, a2, _j, out.offset(2 as i32 as isize));
}

#[inline]
pub unsafe fn Four_Two_Sum(mut a3: f64,
                                      mut a2: f64,
                                      mut a1: f64,
                                      mut a0: f64,
                                      mut b1: f64,
                                      mut b0: f64,
                                      mut out: *mut f64) {
    let mut sum: [f64; 5] = [0.; 5];
    Four_One_Sum(a3, a2, a1, a0, b0, sum.as_mut_ptr());
    *out.offset(0 as i32 as isize) = sum[0 as i32 as usize];
    Four_One_Sum(sum[4 as i32 as usize],
                 sum[3 as i32 as usize],
                 sum[2 as i32 as usize],
                 sum[1 as i32 as usize], b1,
                 out.offset(1 as i32 as isize));
}

#[inline]
pub unsafe fn Four_Four_Sum(mut a3: f64,
                                       mut a2: f64,
                                       mut a1: f64,
                                       mut a0: f64,
                                       mut b4: f64,
                                       mut b3: f64,
                                       mut b1: f64,
                                       mut b0: f64,
                                       mut out: *mut f64) {
    let mut sum: [f64; 6] = [0.; 6];
    Four_Two_Sum(a3, a2, a1, a0, b1, b0, sum.as_mut_ptr());
    *out.offset(1 as i32 as isize) = sum[1 as i32 as usize];
    *out.offset(0 as i32 as isize) = sum[0 as i32 as usize];
    Four_Two_Sum(sum[5 as i32 as usize],
                 sum[4 as i32 as usize],
                 sum[3 as i32 as usize],
                 sum[2 as i32 as usize], b4, b3,
                 out.offset(2 as i32 as isize));
}

#[inline]
pub unsafe fn Eight_One_Sum(mut a7: f64,
                                       mut a6: f64,
                                       mut a5: f64,
                                       mut a4: f64,
                                       mut a3: f64,
                                       mut a2: f64,
                                       mut a1: f64,
                                       mut a0: f64,
                                       mut b: f64,
                                       mut out: *mut f64) {
    let mut sum: [f64; 5] = [0.; 5];
    Four_One_Sum(a3, a2, a1, a0, b, sum.as_mut_ptr());
    *out.offset(3 as i32 as isize) = sum[3 as i32 as usize];
    *out.offset(2 as i32 as isize) = sum[2 as i32 as usize];
    *out.offset(1 as i32 as isize) = sum[1 as i32 as usize];
    *out.offset(0 as i32 as isize) = sum[0 as i32 as usize];
    Four_One_Sum(a7, a6, a5, a4, sum[4 as i32 as usize],
                 out.offset(4 as i32 as isize));
}

#[inline]
pub unsafe fn Eight_Two_Sum(mut a7: f64,
                                       mut a6: f64,
                                       mut a5: f64,
                                       mut a4: f64,
                                       mut a3: f64,
                                       mut a2: f64,
                                       mut a1: f64,
                                       mut a0: f64,
                                       mut b1: f64,
                                       mut b0: f64,
                                       mut out: *mut f64) {
    let mut sum: [f64; 9] = [0.; 9];
    Eight_One_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b0, sum.as_mut_ptr());
    *out.offset(0 as i32 as isize) = sum[0 as i32 as usize];
    Eight_One_Sum(sum[8 as i32 as usize],
                  sum[7 as i32 as usize],
                  sum[6 as i32 as usize],
                  sum[5 as i32 as usize],
                  sum[4 as i32 as usize],
                  sum[3 as i32 as usize],
                  sum[2 as i32 as usize],
                  sum[1 as i32 as usize], b1,
                  out.offset(1 as i32 as isize));
}

#[inline]
pub unsafe fn Eight_Four_Sum(mut a7: f64,
                                        mut a6: f64,
                                        mut a5: f64,
                                        mut a4: f64,
                                        mut a3: f64,
                                        mut a2: f64,
                                        mut a1: f64,
                                        mut a0: f64,
                                        mut b4: f64,
                                        mut b3: f64,
                                        mut b1: f64,
                                        mut b0: f64,
                                        mut out: *mut f64) {
    Eight_Two_Sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0, out);
    Eight_Two_Sum(*out.offset(9 as i32 as isize),
                  *out.offset(8 as i32 as isize),
                  *out.offset(7 as i32 as isize),
                  *out.offset(6 as i32 as isize),
                  *out.offset(5 as i32 as isize),
                  *out.offset(4 as i32 as isize),
                  *out.offset(3 as i32 as isize),
                  *out.offset(2 as i32 as isize), b4, b3,
                  out.offset(2 as i32 as isize));
}
/* Macros for multiplying expansions of various fixed lengths.               */

#[inline]
pub unsafe fn Two_One_Product(a1: f64, a0: f64, b: f64, out: *mut f64) {
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    let mut _k: f64 = 0.;
    Split(b, &mut bhi, &mut blo);
    Two_Product_Presplit(a0, b, bhi, blo, &mut _i,
                         &mut *out.offset(0 as i32 as isize));
    Two_Product_Presplit(a1, b, bhi, blo, &mut _j, &mut _0);
    Two_Sum(_i, _0, &mut _k, &mut *out.offset(1 as i32 as isize));
    let [x, y] = Fast_Two_Sum(_j, _k);
    *out.offset(3 as i32 as isize) = x;
    *out.offset(2 as i32 as isize) = y;
}

#[inline]
pub unsafe fn Four_One_Product(a3: f64,
                                          a2: f64,
                                          a1: f64,
                                          a0: f64,
                                          b: f64,
                                          out: *mut f64) {
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    let mut _k: f64 = 0.;
    Split(b, &mut bhi, &mut blo);
    Two_Product_Presplit(a0, b, bhi, blo, &mut _i,
                         &mut *out.offset(0 as i32 as isize));
    Two_Product_Presplit(a1, b, bhi, blo, &mut _j, &mut _0);
    Two_Sum(_i, _0, &mut _k, &mut *out.offset(1 as i32 as isize));
    let [_i, y] = Fast_Two_Sum(_j, _k);
    *out.offset(2 as i32 as isize) = y;
    Two_Product_Presplit(a2, b, bhi, blo, &mut _j, &mut _0);
    Two_Sum(_i, _0, &mut _k, &mut *out.offset(3 as i32 as isize));
    let [_i, y] = Fast_Two_Sum(_j, _k);
    *out.offset(4 as i32 as isize) = y;
    Two_Product_Presplit(a3, b, bhi, blo, &mut _j, &mut _0);
    Two_Sum(_i, _0, &mut _k, &mut *out.offset(5 as i32 as isize));
    let [x,y] = Fast_Two_Sum(_j, _k);
    *out.offset(7 as i32 as isize) = x;
    *out.offset(6 as i32 as isize) = y;
}

#[inline]
pub unsafe fn Two_Two_Product(mut a1: f64,
                                         mut a0: f64,
                                         mut b1: f64,
                                         mut b0: f64,
                                         mut out: *mut f64) {
    let mut a0hi: f64 = 0.;
    let mut a0lo: f64 = 0.;
    let mut a1hi: f64 = 0.;
    let mut a1lo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    Split(a0, &mut a0hi, &mut a0lo);
    Split(b0, &mut bhi, &mut blo);
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _k: f64 = 0.;
    let mut _l: f64 = 0.;
    let mut _m: f64 = 0.;
    let mut _n: f64 = 0.;
    let mut _0: f64 = 0.;
    let mut _1: f64 = 0.;
    let mut _2: f64 = 0.;
    Two_Product_2Presplit(a0, a0hi, a0lo, b0, bhi, blo, &mut _i,
                          &mut *out.offset(0 as i32 as isize));
    Split(a1, &mut a1hi, &mut a1lo);
    Two_Product_2Presplit(a1, a1hi, a1lo, b0, bhi, blo, &mut _j, &mut _0);
    Two_Sum(_i, _0, &mut _k, &mut _1);
    let [mut _l, mut _2] = Fast_Two_Sum(_j, _k);
    Split(b1, &mut bhi, &mut blo);
    Two_Product_2Presplit(a0, a0hi, a0lo, b1, bhi, blo, &mut _i, &mut _0);
    Two_Sum(_1, _0, &mut _k, &mut *out.offset(1 as i32 as isize));
    Two_Sum(_2, _k, &mut _j, &mut _1);
    Two_Sum(_l, _j, &mut _m, &mut _2);
    Two_Product_2Presplit(a1, a1hi, a1lo, b1, bhi, blo, &mut _j, &mut _0);
    Two_Sum(_i, _0, &mut _n, &mut _0);
    Two_Sum(_1, _0, &mut _i, &mut *out.offset(2 as i32 as isize));
    Two_Sum(_2, _i, &mut _k, &mut _1);
    Two_Sum(_m, _k, &mut _l, &mut _2);
    Two_Sum(_j, _n, &mut _k, &mut _0);
    Two_Sum(_1, _0, &mut _j, &mut *out.offset(3 as i32 as isize));
    Two_Sum(_2, _j, &mut _i, &mut _1);
    Two_Sum(_l, _i, &mut _m, &mut _2);
    Two_Sum(_1, _k, &mut _i, &mut *out.offset(4 as i32 as isize));
    Two_Sum(_2, _i, &mut _k, &mut *out.offset(5 as i32 as isize));
    Two_Sum(_m, _k, &mut *out.offset(7 as i32 as isize),
            &mut *out.offset(6 as i32 as isize));
}

/* An expansion of length two can be squared more quickly than finding the   */
/*   product of two different expansions of length two, and the result is    */
/*   guaranteed to have no more than six (rather than eight) components.     */

#[inline]
pub unsafe fn Two_Square(mut a1: f64,
                                    mut a0: f64,
                                    mut out: *mut f64) {
    let mut _j: f64 = 0.;
    Square(a0, &mut _j, &mut *out.offset(0 as i32 as isize));
    let mut _0: f64 = a0 + a0;
    let mut _k: f64 = 0.;
    let mut _1: f64 = 0.;
    let mut _2: f64 = 0.;
    let mut _l: f64 = 0.;
    Two_Product(a1, _0, &mut _k, &mut _1);
    let mut sum: [f64; 3] = [0.; 3];
    Two_One_Sum(_k, _1, _j, sum.as_mut_ptr());
    _l = sum[2 as i32 as usize];
    _2 = sum[1 as i32 as usize];
    *out.offset(1 as i32 as isize) = sum[0 as i32 as usize];
    Square(a1, &mut _j, &mut _1);
    Two_Two_Sum(_j, _1, _l, _2, out.offset(2 as i32 as isize));
}

/* ****************************************************************************/
/*                                                                           */
/*  grow_expansion()   Add a scalar to an expansion.                         */
/*                                                                           */
/*  Sets h = e + b.  See the long version of my paper for details.           */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn grow_expansion(mut elen: i32,
                                        mut e: *const f64,
                                        mut b: f64,
                                        mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut enow: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    Q = b;
    eindex = 0 as i32;
    while eindex < elen {
        enow = *e.offset(eindex as isize);
        Two_Sum(Q, enow, &mut Qnew, &mut *h.offset(eindex as isize));
        Q = Qnew;
        eindex += 1
    }
    *h.offset(eindex as isize) = Q;
    return eindex + 1 as i32;
}
/* ****************************************************************************/
/*                                                                           */
/*  grow_expansion_zeroelim()   Add a scalar to an expansion, eliminating    */
/*                              zero components from the output expansion.   */
/*                                                                           */
/*  Sets h = e + b.  See the long version of my paper for details.           */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn grow_expansion_zeroelim(mut elen: i32,
                                                 mut e: *const f64,
                                                 mut b: f64,
                                                 mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut hh: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut enow: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    hindex = 0 as i32;
    Q = b;
    eindex = 0 as i32;
    while eindex < elen {
        enow = *e.offset(eindex as isize);
        Two_Sum(Q, enow, &mut Qnew, &mut hh);
        Q = Qnew;
        if hh != 0.0f64 {
            let fresh0 = hindex;
            hindex = hindex + 1;
            *h.offset(fresh0 as isize) = hh
        }
        eindex += 1
    }
    if Q != 0.0f64 || hindex == 0 as i32 {
        let fresh1 = hindex;
        hindex = hindex + 1;
        *h.offset(fresh1 as isize) = Q
    }
    return hindex;
}
/* ****************************************************************************/
/*                                                                           */
/*  expansion_sum()   Sum two expansions.                                    */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
/*  strongly nonoverlapping property.                                        */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn expansion_sum(mut elen: i32,
                                       mut e: *const f64,
                                       mut flen: i32,
                                       mut f: *const f64,
                                       mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut findex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut hlast: i32 = 0;
    let mut hnow: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    Q = *f.offset(0 as i32 as isize);
    hindex = 0 as i32;
    while hindex < elen {
        hnow = *e.offset(hindex as isize);
        Two_Sum(Q, hnow, &mut Qnew, &mut *h.offset(hindex as isize));
        Q = Qnew;
        hindex += 1
    }
    *h.offset(hindex as isize) = Q;
    hlast = hindex;
    findex = 1 as i32;
    while findex < flen {
        Q = *f.offset(findex as isize);
        hindex = findex;
        while hindex <= hlast {
            hnow = *h.offset(hindex as isize);
            Two_Sum(Q, hnow, &mut Qnew, &mut *h.offset(hindex as isize));
            Q = Qnew;
            hindex += 1
        }
        hlast += 1;
        *h.offset(hlast as isize) = Q;
        findex += 1
    }
    return hlast + 1 as i32;
}
/* ****************************************************************************/
/*                                                                           */
/*  expansion_sum_zeroelim1()   Sum two expansions, eliminating zero         */
/*                              components from the output expansion.        */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
/*  strongly nonoverlapping property.                                        */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn expansion_sum_zeroelim1(mut elen: i32,
                                                 mut e: *const f64,
                                                 mut flen: i32,
                                                 mut f: *const f64,
                                                 mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut index: i32 = 0;
    let mut findex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut hlast: i32 = 0;
    let mut hnow: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    Q = *f.offset(0 as i32 as isize);
    hindex = 0 as i32;
    while hindex < elen {
        hnow = *e.offset(hindex as isize);
        Two_Sum(Q, hnow, &mut Qnew, &mut *h.offset(hindex as isize));
        Q = Qnew;
        hindex += 1
    }
    *h.offset(hindex as isize) = Q;
    hlast = hindex;
    findex = 1 as i32;
    while findex < flen {
        Q = *f.offset(findex as isize);
        hindex = findex;
        while hindex <= hlast {
            hnow = *h.offset(hindex as isize);
            Two_Sum(Q, hnow, &mut Qnew, &mut *h.offset(hindex as isize));
            Q = Qnew;
            hindex += 1
        }
        hlast += 1;
        *h.offset(hlast as isize) = Q;
        findex += 1
    }
    hindex = -(1 as i32);
    index = 0 as i32;
    while index <= hlast {
        hnow = *h.offset(index as isize);
        if hnow != 0.0f64 { hindex += 1; *h.offset(hindex as isize) = hnow }
        index += 1
    }
    if hindex == -(1 as i32) {
        return 1 as i32
    } else { return hindex + 1 as i32 };
}
/* ****************************************************************************/
/*                                                                           */
/*  expansion_sum_zeroelim2()   Sum two expansions, eliminating zero         */
/*                              components from the output expansion.        */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the nonadjacent property as well.  (That is,   */
/*  if e has one of these properties, so will h.)  Does NOT maintain the     */
/*  strongly nonoverlapping property.                                        */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn expansion_sum_zeroelim2(mut elen: i32,
                                                 mut e: *const f64,
                                                 mut flen: i32,
                                                 mut f: *const f64,
                                                 mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut hh: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut findex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut hlast: i32 = 0;
    let mut enow: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    hindex = 0 as i32;
    Q = *f.offset(0 as i32 as isize);
    eindex = 0 as i32;
    while eindex < elen {
        enow = *e.offset(eindex as isize);
        Two_Sum(Q, enow, &mut Qnew, &mut hh);
        Q = Qnew;
        if hh != 0.0f64 {
            let fresh2 = hindex;
            hindex = hindex + 1;
            *h.offset(fresh2 as isize) = hh
        }
        eindex += 1
    }
    *h.offset(hindex as isize) = Q;
    hlast = hindex;
    findex = 1 as i32;
    while findex < flen {
        hindex = 0 as i32;
        Q = *f.offset(findex as isize);
        eindex = 0 as i32;
        while eindex <= hlast {
            enow = *h.offset(eindex as isize);
            Two_Sum(Q, enow, &mut Qnew, &mut hh);
            Q = Qnew;
            if hh != 0 as i32 as f64 {
                let fresh3 = hindex;
                hindex = hindex + 1;
                *h.offset(fresh3 as isize) = hh
            }
            eindex += 1
        }
        *h.offset(hindex as isize) = Q;
        hlast = hindex;
        findex += 1
    }
    return hlast + 1 as i32;
}
/* ****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum()   Sum two expansions.                               */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn fast_expansion_sum(mut elen: i32,
                                            mut e: *const f64,
                                            mut flen: i32,
                                            mut f: *const f64,
                                            mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut findex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut enow: f64 = 0.;
    let mut fnow: f64 = 0.;
    enow = *e.offset(0 as i32 as isize);
    fnow = *f.offset(0 as i32 as isize);
    findex = 0 as i32;
    eindex = findex;
    if (fnow > enow) as i32 == (fnow > -enow) as i32 {
        Q = enow;
        eindex += 1;
        enow = *e.offset(eindex as isize)
    } else { Q = fnow; findex += 1; fnow = *f.offset(findex as isize) }
    hindex = 0 as i32;
    if eindex < elen && findex < flen {
        if (fnow > enow) as i32 == (fnow > -enow) as i32 {
            let [x, y] = Fast_Two_Sum(enow, Q);
            Qnew = x;
            *h.offset(0 as i32 as isize) = y;
            eindex += 1;
            enow = *e.offset(eindex as isize)
        } else {
            let [x, y] = Fast_Two_Sum(fnow, Q);
            Qnew = x;
            *h.offset(0 as i32 as isize) = y;
            findex += 1;
            fnow = *f.offset(findex as isize)
        }
        Q = Qnew;
        hindex = 1 as i32;
        while eindex < elen && findex < flen {
            if (fnow > enow) as i32 == (fnow > -enow) as i32 {
                Two_Sum(Q, enow, &mut Qnew, &mut *h.offset(hindex as isize));
                eindex += 1;
                enow = *e.offset(eindex as isize)
            } else {
                Two_Sum(Q, fnow, &mut Qnew, &mut *h.offset(hindex as isize));
                findex += 1;
                fnow = *f.offset(findex as isize)
            }
            Q = Qnew;
            hindex += 1
        }
    }
    while eindex < elen {
        Two_Sum(Q, enow, &mut Qnew, &mut *h.offset(hindex as isize));
        eindex += 1;
        enow = *e.offset(eindex as isize);
        Q = Qnew;
        hindex += 1
    }
    while findex < flen {
        Two_Sum(Q, fnow, &mut Qnew, &mut *h.offset(hindex as isize));
        findex += 1;
        fnow = *f.offset(findex as isize);
        Q = Qnew;
        hindex += 1
    }
    *h.offset(hindex as isize) = Q;
    return hindex + 1 as i32;
}
/* ****************************************************************************/
/*                                                                           */
/*  fast_expansion_sum_zeroelim()   Sum two expansions, eliminating zero     */
/*                                  components from the output expansion.    */
/*                                                                           */
/*  Sets h = e + f.  See the long version of my paper for details.           */
/*                                                                           */
/*  If round-to-even is used (as with IEEE 754), maintains the strongly      */
/*  nonoverlapping property.  (That is, if e is strongly nonoverlapping, h   */
/*  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent      */
/*  properties.                                                              */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn fast_expansion_sum_zeroelim(mut elen: i32,
                                                     mut e:
                                                         *const f64,
                                                     mut flen: i32,
                                                     mut f:
                                                         *const f64,
                                                     mut h:
                                                         *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut hh: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut findex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut enow: f64 = 0.;
    let mut fnow: f64 = 0.;
    enow = *e.offset(0 as i32 as isize);
    fnow = *f.offset(0 as i32 as isize);
    findex = 0 as i32;
    eindex = findex;
    if (fnow > enow) as i32 == (fnow > -enow) as i32 {
        Q = enow;
        eindex += 1;
        enow = *e.offset(eindex as isize)
    } else { Q = fnow; findex += 1; fnow = *f.offset(findex as isize) }
    hindex = 0 as i32;
    if eindex < elen && findex < flen {
        if (fnow > enow) as i32 == (fnow > -enow) as i32 {
            let [x, y] = Fast_Two_Sum(enow, Q);
            Qnew = x;
            hh = y;
            eindex += 1;
            enow = *e.offset(eindex as isize)
        } else {
            let [x, y] = Fast_Two_Sum(fnow, Q);
            Qnew = x;
            hh = y;
            findex += 1;
            fnow = *f.offset(findex as isize)
        }
        Q = Qnew;
        if hh != 0.0f64 {
            let fresh4 = hindex;
            hindex = hindex + 1;
            *h.offset(fresh4 as isize) = hh
        }
        while eindex < elen && findex < flen {
            if (fnow > enow) as i32 == (fnow > -enow) as i32 {
                Two_Sum(Q, enow, &mut Qnew, &mut hh);
                eindex += 1;
                enow = *e.offset(eindex as isize)
            } else {
                Two_Sum(Q, fnow, &mut Qnew, &mut hh);
                findex += 1;
                fnow = *f.offset(findex as isize)
            }
            Q = Qnew;
            if hh != 0.0f64 {
                let fresh5 = hindex;
                hindex = hindex + 1;
                *h.offset(fresh5 as isize) = hh
            }
        }
    }
    while eindex < elen {
        Two_Sum(Q, enow, &mut Qnew, &mut hh);
        eindex += 1;
        enow = *e.offset(eindex as isize);
        Q = Qnew;
        if hh != 0.0f64 {
            let fresh6 = hindex;
            hindex = hindex + 1;
            *h.offset(fresh6 as isize) = hh
        }
    }
    while findex < flen {
        Two_Sum(Q, fnow, &mut Qnew, &mut hh);
        findex += 1;
        fnow = *f.offset(findex as isize);
        Q = Qnew;
        if hh != 0.0f64 {
            let fresh7 = hindex;
            hindex = hindex + 1;
            *h.offset(fresh7 as isize) = hh
        }
    }
    if Q != 0.0f64 || hindex == 0 as i32 {
        let fresh8 = hindex;
        hindex = hindex + 1;
        *h.offset(fresh8 as isize) = Q
    }
    return hindex;
}
/* ****************************************************************************/
/*                                                                           */
/*  linear_expansion_sum()   Sum two expansions.                             */
/*                                                                           */
/*  Sets h = e + f.  See either version of my paper for details.             */
/*                                                                           */
/*  Maintains the nonoverlapping property.  (That is, if e is                */
/*  nonoverlapping, h will be also.)                                         */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn linear_expansion_sum(mut elen: i32,
                                              mut e: *const f64,
                                              mut flen: i32,
                                              mut f: *const f64,
                                              mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut q: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut R: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut findex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut enow: f64 = 0.;
    let mut fnow: f64 = 0.;
    let mut g0: f64 = 0.;
    enow = *e.offset(0 as i32 as isize);
    fnow = *f.offset(0 as i32 as isize);
    findex = 0 as i32;
    eindex = findex;
    if (fnow > enow) as i32 == (fnow > -enow) as i32 {
        g0 = enow;
        eindex += 1;
        enow = *e.offset(eindex as isize)
    } else { g0 = fnow; findex += 1; fnow = *f.offset(findex as isize) }
    if eindex < elen &&
           (findex >= flen ||
                (fnow > enow) as i32 == (fnow > -enow) as i32)
       {
        let [x, y] = Fast_Two_Sum(enow, g0);
        Qnew = x;
        q = y;
        eindex += 1;
        enow = *e.offset(eindex as isize)
    } else {
        let [x, y] = Fast_Two_Sum(fnow, g0);
        Qnew = x;
        q = y;
        findex += 1;
        fnow = *f.offset(findex as isize)
    }
    Q = Qnew;
    hindex = 0 as i32;
    while hindex < elen + flen - 2 as i32 {
        if eindex < elen &&
               (findex >= flen ||
                    (fnow > enow) as i32 ==
                        (fnow > -enow) as i32) {
            let [x, y] = Fast_Two_Sum(enow, q);
            R = x;
            *h.offset(hindex as isize) = y;
            eindex += 1;
            enow = *e.offset(eindex as isize)
        } else {
            let [x, y] = Fast_Two_Sum(fnow, q);
            R = x;
            *h.offset(hindex as isize) = y;
            findex += 1;
            fnow = *f.offset(findex as isize)
        }
        Two_Sum(Q, R, &mut Qnew, &mut q);
        Q = Qnew;
        hindex += 1
    }
    *h.offset(hindex as isize) = q;
    *h.offset((hindex + 1 as i32) as isize) = Q;
    return hindex + 2 as i32;
}
/* ****************************************************************************/
/*                                                                           */
/*  linear_expansion_sum_zeroelim()   Sum two expansions, eliminating zero   */
/*                                    components from the output expansion.  */
/*                                                                           */
/*  Sets h = e + f.  See either version of my paper for details.             */
/*                                                                           */
/*  Maintains the nonoverlapping property.  (That is, if e is                */
/*  nonoverlapping, h will be also.)                                         */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn linear_expansion_sum_zeroelim(mut elen: i32,
                                                       mut e:
                                                           *const f64,
                                                       mut flen: i32,
                                                       mut f:
                                                           *const f64,
                                                       mut h:
                                                           *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut q: f64 = 0.;
    let mut hh: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut R: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut findex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut count: i32 = 0;
    let mut enow: f64 = 0.;
    let mut fnow: f64 = 0.;
    let mut g0: f64 = 0.;
    enow = *e.offset(0 as i32 as isize);
    fnow = *f.offset(0 as i32 as isize);
    findex = 0 as i32;
    eindex = findex;
    hindex = 0 as i32;
    if (fnow > enow) as i32 == (fnow > -enow) as i32 {
        g0 = enow;
        eindex += 1;
        enow = *e.offset(eindex as isize)
    } else { g0 = fnow; findex += 1; fnow = *f.offset(findex as isize) }
    if eindex < elen &&
           (findex >= flen ||
                (fnow > enow) as i32 == (fnow > -enow) as i32)
       {
        let [x, y] = Fast_Two_Sum(enow, g0);
        Qnew = x;
        q = y;
        eindex += 1;
        enow = *e.offset(eindex as isize)
    } else {
        let [x, y] = Fast_Two_Sum(fnow, g0);
        Qnew = x;
        q = y;
        findex += 1;
        fnow = *f.offset(findex as isize)
    }
    Q = Qnew;
    count = 2 as i32;
    while count < elen + flen {
        if eindex < elen &&
               (findex >= flen ||
                    (fnow > enow) as i32 ==
                        (fnow > -enow) as i32) {
            let [x, y] = Fast_Two_Sum(enow, q);
            R = x;
            hh = y;
            eindex += 1;
            enow = *e.offset(eindex as isize)
        } else {
            let [x, y] = Fast_Two_Sum(fnow, q);
            R = x;
            hh = y;
            findex += 1;
            fnow = *f.offset(findex as isize)
        }
        Two_Sum(Q, R, &mut Qnew, &mut q);
        Q = Qnew;
        if hh != 0 as i32 as f64 {
            let fresh9 = hindex;
            hindex = hindex + 1;
            *h.offset(fresh9 as isize) = hh
        }
        count += 1
    }
    if q != 0 as i32 as f64 {
        let fresh10 = hindex;
        hindex = hindex + 1;
        *h.offset(fresh10 as isize) = q
    }
    if Q != 0.0f64 || hindex == 0 as i32 {
        let fresh11 = hindex;
        hindex = hindex + 1;
        *h.offset(fresh11 as isize) = Q
    }
    return hindex;
}
/* ****************************************************************************/
/*                                                                           */
/*  scale_expansion()   Multiply an expansion by a scalar.                   */
/*                                                                           */
/*  Sets h = be.  See either version of my paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn scale_expansion(mut elen: i32,
                                         mut e: *const f64,
                                         mut b: f64,
                                         mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut sum: f64 = 0.;
    let mut product1: f64 = 0.;
    let mut product0: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut enow: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    Split(b, &mut bhi, &mut blo);
    Two_Product_Presplit(*e.offset(0 as i32 as isize), b, bhi, blo,
                         &mut Q, &mut *h.offset(0 as i32 as isize));
    hindex = 1 as i32;
    eindex = 1 as i32;
    while eindex < elen {
        enow = *e.offset(eindex as isize);
        Two_Product_Presplit(enow, b, bhi, blo, &mut product1, &mut product0);
        Two_Sum(Q, product0, &mut sum, &mut *h.offset(hindex as isize));
        hindex += 1;
        Two_Sum(product1, sum, &mut Q, &mut *h.offset(hindex as isize));
        hindex += 1;
        eindex += 1
    }
    *h.offset(hindex as isize) = Q;
    return elen + elen;
}
/* ****************************************************************************/
/*                                                                           */
/*  scale_expansion_zeroelim()   Multiply an expansion by a scalar,          */
/*                               eliminating zero components from the        */
/*                               output expansion.                           */
/*                                                                           */
/*  Sets h = be.  See either version of my paper for details.                */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent    */
/*  properties as well.  (That is, if e has one of these properties, so      */
/*  will h.)                                                                 */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn scale_expansion_zeroelim(mut elen: i32,
                                                  mut e:
                                                      *const f64,
                                                  mut b: f64,
                                                  mut h: *mut f64)
 -> i32 {
    let mut Q: f64 = 0.;
    let mut sum: f64 = 0.;
    let mut hh: f64 = 0.;
    let mut product1: f64 = 0.;
    let mut product0: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut enow: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    Split(b, &mut bhi, &mut blo);
    Two_Product_Presplit(*e.offset(0 as i32 as isize), b, bhi, blo,
                         &mut Q, &mut hh);
    hindex = 0 as i32;
    if hh != 0 as i32 as f64 {
        let fresh12 = hindex;
        hindex = hindex + 1;
        *h.offset(fresh12 as isize) = hh
    }
    eindex = 1 as i32;
    while eindex < elen {
        enow = *e.offset(eindex as isize);
        Two_Product_Presplit(enow, b, bhi, blo, &mut product1, &mut product0);
        Two_Sum(Q, product0, &mut sum, &mut hh);
        if hh != 0 as i32 as f64 {
            let fresh13 = hindex;
            hindex = hindex + 1;
            *h.offset(fresh13 as isize) = hh
        }
        let [x, y] = Fast_Two_Sum(product1, sum);
        Q = x;
        hh = y;
        if hh != 0 as i32 as f64 {
            let fresh14 = hindex;
            hindex = hindex + 1;
            *h.offset(fresh14 as isize) = hh
        }
        eindex += 1
    }
    if Q != 0.0f64 || hindex == 0 as i32 {
        let fresh15 = hindex;
        hindex = hindex + 1;
        *h.offset(fresh15 as isize) = Q
    }
    return hindex;
}
/* ****************************************************************************/
/*                                                                           */
/*  compress()   Compress an expansion.                                      */
/*                                                                           */
/*  See the long version of my paper for details.                            */
/*                                                                           */
/*  Maintains the nonoverlapping property.  If round-to-even is used (as     */
/*  with IEEE 754), then any nonoverlapping expansion is converted to a      */
/*  nonadjacent expansion.                                                   */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn compress(mut elen: i32,
                                  mut e: *const f64,
                                  mut h: *mut f64) -> i32 {
    let mut Q: f64 = 0.;
    let mut q: f64 = 0.;
    let mut Qnew: f64 = 0.;
    let mut eindex: i32 = 0;
    let mut hindex: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut enow: f64 = 0.;
    let mut hnow: f64 = 0.;
    let mut top: i32 = 0;
    let mut bottom: i32 = 0;
    bottom = elen - 1 as i32;
    Q = *e.offset(bottom as isize);
    eindex = elen - 2 as i32;
    while eindex >= 0 as i32 {
        enow = *e.offset(eindex as isize);
        let [x, y] = Fast_Two_Sum(Q, enow);
        Qnew = x;
        q = y;
        if q != 0 as i32 as f64 {
            let fresh16 = bottom;
            bottom = bottom - 1;
            *h.offset(fresh16 as isize) = Qnew;
            Q = q
        } else { Q = Qnew }
        eindex -= 1
    }
    top = 0 as i32;
    hindex = bottom + 1 as i32;
    while hindex < elen {
        hnow = *h.offset(hindex as isize);
        let mut sum: [f64; 2] = [0.; 2];
        let [x, y] = Fast_Two_Sum(hnow, Q);
        Qnew = x;
        q = y;
        if q != 0 as i32 as f64 {
            let fresh17 = top;
            top = top + 1;
            *h.offset(fresh17 as isize) = q
        }
        Q = Qnew;
        hindex += 1
    }
    *h.offset(top as isize) = Q;
    return top + 1 as i32;
}
/* ****************************************************************************/
/*                                                                           */
/*  estimate()   Produce a one-word estimate of an expansion's value.        */
/*                                                                           */
/*  See either version of my paper for details.                              */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn estimate(mut elen: i32,
                                  mut e: *const f64)
 -> f64 {
    let mut Q: f64 = 0.;
    let mut eindex: i32 = 0;
    Q = *e.offset(0 as i32 as isize);
    eindex = 1 as i32;
    while eindex < elen { Q += *e.offset(eindex as isize); eindex += 1 }
    return Q;
}
/* ****************************************************************************/
/*                                                                           */
/*  orient2dfast()   Approximate 2D orientation test.  Nonrobust.            */
/*  orient2dexact()   Exact 2D orientation test.  Robust.                    */
/*  orient2dslow()   Another exact 2D orientation test.  Robust.             */
/*  orient2d()   Adaptive exact 2D orientation test.  Robust.                */
/*                                                                           */
/*               Return a positive value if the points pa, pb, and pc occur  */
/*               in counterclockwise order; a negative value if they occur   */
/*               in clockwise order; and zero if they are collinear.  The    */
/*               result is also a rough approximation of twice the signed    */
/*               area of the triangle defined by the three points.           */
/*                                                                           */
/*  Only the first and last routine should be used; the middle two are for   */
/*  timings.                                                                 */
/*                                                                           */
/*  The last three use exact arithmetic to ensure a correct answer.  The     */
/*  result returned is the determinant of a matrix.  In orient2d() only,     */
/*  this determinant is computed adaptively, in the sense that exact         */
/*  arithmetic is used only to the degree it is needed to ensure that the    */
/*  returned value has the correct sign.  Hence, orient2d() is usually quite */
/*  fast, but will run more slowly when the input points are collinear or    */
/*  nearly so.                                                               */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn orient2dfast(mut pa: *const f64,
                                      mut pb: *const f64,
                                      mut pc: *const f64)
 -> f64 {
    let mut acx: f64 = 0.;
    let mut bcx: f64 = 0.;
    let mut acy: f64 = 0.;
    let mut bcy: f64 = 0.;
    acx =
        *pa.offset(0 as i32 as isize) -
            *pc.offset(0 as i32 as isize);
    bcx =
        *pb.offset(0 as i32 as isize) -
            *pc.offset(0 as i32 as isize);
    acy =
        *pa.offset(1 as i32 as isize) -
            *pc.offset(1 as i32 as isize);
    bcy =
        *pb.offset(1 as i32 as isize) -
            *pc.offset(1 as i32 as isize);
    return acx * bcy - acy * bcx;
}

pub unsafe fn orient2dexact(mut pa: *const f64,
                                       mut pb: *const f64,
                                       mut pc: *const f64)
 -> f64 {
    let mut axby1: f64 = 0.;
    let mut axcy1: f64 = 0.;
    let mut bxcy1: f64 = 0.;
    let mut bxay1: f64 = 0.;
    let mut cxay1: f64 = 0.;
    let mut cxby1: f64 = 0.;
    let mut axby0: f64 = 0.;
    let mut axcy0: f64 = 0.;
    let mut bxcy0: f64 = 0.;
    let mut bxay0: f64 = 0.;
    let mut cxay0: f64 = 0.;
    let mut cxby0: f64 = 0.;
    let mut aterms: [f64; 4] = [0.; 4];
    let mut bterms: [f64; 4] = [0.; 4];
    let mut cterms: [f64; 4] = [0.; 4];
    let mut aterms3: f64 = 0.;
    let mut bterms3: f64 = 0.;
    let mut cterms3: f64 = 0.;
    let mut v: [f64; 8] = [0.; 8];
    let mut w: [f64; 12] = [0.; 12];
    let mut vlength: i32 = 0;
    let mut wlength: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    Two_Product(*pa.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut axby1,
                &mut axby0);
    Two_Product(*pa.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut axcy1,
                &mut axcy0);
    Two_Two_Diff(axby1, axby0, axcy1, axcy0, aterms.as_mut_ptr());
    Two_Product(*pb.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut bxcy1,
                &mut bxcy0);
    Two_Product(*pb.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut bxay1,
                &mut bxay0);
    Two_Two_Diff(bxcy1, bxcy0, bxay1, bxay0, bterms.as_mut_ptr());
    Two_Product(*pc.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut cxay1,
                &mut cxay0);
    Two_Product(*pc.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut cxby1,
                &mut cxby0);
    Two_Two_Diff(cxay1, cxay0, cxby1, cxby0, cterms.as_mut_ptr());
    vlength =
        fast_expansion_sum_zeroelim(4 as i32, aterms.as_mut_ptr(),
                                    4 as i32, bterms.as_mut_ptr(),
                                    v.as_mut_ptr());
    wlength =
        fast_expansion_sum_zeroelim(vlength, v.as_mut_ptr(), 4 as i32,
                                    cterms.as_mut_ptr(), w.as_mut_ptr());
    return w[(wlength - 1 as i32) as usize];
}

pub unsafe fn orient2dslow(mut pa: *const f64,
                                      mut pb: *const f64,
                                      mut pc: *const f64)
 -> f64 {
    let mut acx: f64 = 0.;
    let mut acy: f64 = 0.;
    let mut bcx: f64 = 0.;
    let mut bcy: f64 = 0.;
    let mut acxtail: f64 = 0.;
    let mut acytail: f64 = 0.;
    let mut bcxtail: f64 = 0.;
    let mut bcytail: f64 = 0.;
    let mut negate: f64 = 0.;
    let mut negatetail: f64 = 0.;
    let mut axby: [f64; 8] = [0.; 8];
    let mut bxay: [f64; 8] = [0.; 8];
    let mut deter: [f64; 16] = [0.; 16];
    let mut deterlen: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut a0hi: f64 = 0.;
    let mut a0lo: f64 = 0.;
    let mut a1hi: f64 = 0.;
    let mut a1lo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _k: f64 = 0.;
    let mut _l: f64 = 0.;
    let mut _m: f64 = 0.;
    let mut _n: f64 = 0.;
    let mut _0: f64 = 0.;
    let mut _1: f64 = 0.;
    let mut _2: f64 = 0.;
    Two_Diff(*pa.offset(0 as i32 as isize),
             *pc.offset(0 as i32 as isize), &mut acx, &mut acxtail);
    Two_Diff(*pa.offset(1 as i32 as isize),
             *pc.offset(1 as i32 as isize), &mut acy, &mut acytail);
    Two_Diff(*pb.offset(0 as i32 as isize),
             *pc.offset(0 as i32 as isize), &mut bcx, &mut bcxtail);
    Two_Diff(*pb.offset(1 as i32 as isize),
             *pc.offset(1 as i32 as isize), &mut bcy, &mut bcytail);
    Two_Two_Product(acx, acxtail, bcy, bcytail, axby.as_mut_ptr());
    negate = -acy;
    negatetail = -acytail;
    Two_Two_Product(bcx, bcxtail, negate, negatetail, bxay.as_mut_ptr());
    deterlen =
        fast_expansion_sum_zeroelim(8 as i32, axby.as_mut_ptr(),
                                    8 as i32, bxay.as_mut_ptr(),
                                    deter.as_mut_ptr());
    return deter[(deterlen - 1 as i32) as usize];
}

pub unsafe fn orient2dadapt(mut pa: *const f64,
                                       mut pb: *const f64,
                                       mut pc: *const f64,
                                       mut detsum: f64)
 -> f64 {
    let mut acx: f64 = 0.;
    let mut acy: f64 = 0.;
    let mut bcx: f64 = 0.;
    let mut bcy: f64 = 0.;
    let mut acxtail: f64 = 0.;
    let mut acytail: f64 = 0.;
    let mut bcxtail: f64 = 0.;
    let mut bcytail: f64 = 0.;
    let mut detleft: f64 = 0.;
    let mut detright: f64 = 0.;
    let mut detlefttail: f64 = 0.;
    let mut detrighttail: f64 = 0.;
    let mut det: f64 = 0.;
    let mut errbound: f64 = 0.;
    let mut B: [f64; 4] = [0.; 4];
    let mut C1: [f64; 8] = [0.; 8];
    let mut C2: [f64; 12] = [0.; 12];
    let mut D: [f64; 16] = [0.; 16];
    let mut B3: f64 = 0.;
    let mut C1length: i32 = 0;
    let mut C2length: i32 = 0;
    let mut Dlength: i32 = 0;
    let mut u: [f64; 4] = [0.; 4];
    let mut u3: f64 = 0.;
    let mut s1: f64 = 0.;
    let mut t1: f64 = 0.;
    let mut s0: f64 = 0.;
    let mut t0: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    acx =
        *pa.offset(0 as i32 as isize) -
            *pc.offset(0 as i32 as isize);
    bcx =
        *pb.offset(0 as i32 as isize) -
            *pc.offset(0 as i32 as isize);
    acy =
        *pa.offset(1 as i32 as isize) -
            *pc.offset(1 as i32 as isize);
    bcy =
        *pb.offset(1 as i32 as isize) -
            *pc.offset(1 as i32 as isize);
    Two_Product(acx, bcy, &mut detleft, &mut detlefttail);
    Two_Product(acy, bcx, &mut detright, &mut detrighttail);
    Two_Two_Diff(detleft, detlefttail, detright, detrighttail,
                 B.as_mut_ptr());
    det = estimate(4 as i32, B.as_mut_ptr());
    errbound = PARAMS.ccwerrboundB * detsum;
    if det >= errbound || -det >= errbound { return det }
    acxtail =
        Two_Diff_Tail(*pa.offset(0 as i32 as isize),
                      *pc.offset(0 as i32 as isize), acx);
    bcxtail =
        Two_Diff_Tail(*pb.offset(0 as i32 as isize),
                      *pc.offset(0 as i32 as isize), bcx);
    acytail =
        Two_Diff_Tail(*pa.offset(1 as i32 as isize),
                      *pc.offset(1 as i32 as isize), acy);
    bcytail =
        Two_Diff_Tail(*pb.offset(1 as i32 as isize),
                      *pc.offset(1 as i32 as isize), bcy);
    if acxtail == 0.0f64 && acytail == 0.0f64 && bcxtail == 0.0f64 &&
           bcytail == 0.0f64 {
        return det
    }
    errbound = PARAMS.ccwerrboundC * detsum + PARAMS.resulterrbound * Absolute(det);
    det += acx * bcytail + bcy * acxtail - (acy * bcxtail + bcx * acytail);
    if det >= errbound || -det >= errbound { return det }
    Two_Product(acxtail, bcy, &mut s1, &mut s0);
    Two_Product(acytail, bcx, &mut t1, &mut t0);
    Two_Two_Diff(s1, s0, t1, t0, u.as_mut_ptr());
    C1length =
        fast_expansion_sum_zeroelim(4 as i32, B.as_mut_ptr(),
                                    4 as i32, u.as_mut_ptr(),
                                    C1.as_mut_ptr());
    Two_Product(acx, bcytail, &mut s1, &mut s0);
    Two_Product(acy, bcxtail, &mut t1, &mut t0);
    Two_Two_Diff(s1, s0, t1, t0, u.as_mut_ptr());
    C2length =
        fast_expansion_sum_zeroelim(C1length, C1.as_mut_ptr(),
                                    4 as i32, u.as_mut_ptr(),
                                    C2.as_mut_ptr());
    Two_Product(acxtail, bcytail, &mut s1, &mut s0);
    Two_Product(acytail, bcxtail, &mut t1, &mut t0);
    Two_Two_Diff(s1, s0, t1, t0, u.as_mut_ptr());
    Dlength =
        fast_expansion_sum_zeroelim(C2length, C2.as_mut_ptr(),
                                    4 as i32, u.as_mut_ptr(),
                                    D.as_mut_ptr());
    return D[(Dlength - 1 as i32) as usize];
}

pub unsafe fn orient2d(mut pa: *const f64,
                                  mut pb: *const f64,
                                  mut pc: *const f64)
 -> f64 {
    let mut detleft: f64 = 0.;
    let mut detright: f64 = 0.;
    let mut det: f64 = 0.;
    let mut detsum: f64 = 0.;
    let mut errbound: f64 = 0.;
    detleft =
        (*pa.offset(0 as i32 as isize) -
             *pc.offset(0 as i32 as isize)) *
            (*pb.offset(1 as i32 as isize) -
                 *pc.offset(1 as i32 as isize));
    detright =
        (*pa.offset(1 as i32 as isize) -
             *pc.offset(1 as i32 as isize)) *
            (*pb.offset(0 as i32 as isize) -
                 *pc.offset(0 as i32 as isize));
    det = detleft - detright;
    if detleft > 0.0f64 {
        if detright <= 0.0f64 {
            return det
        } else { detsum = detleft + detright }
    } else if detleft < 0.0f64 {
        if detright >= 0.0f64 {
            return det
        } else { detsum = -detleft - detright }
    } else { return det }
    errbound = PARAMS.ccwerrboundA * detsum;
    if det >= errbound || -det >= errbound { return det }
    return orient2dadapt(pa, pb, pc, detsum);
}
/* ****************************************************************************/
/*                                                                           */
/*  orient3dfast()   Approximate 3D orientation test.  Nonrobust.            */
/*  orient3dexact()   Exact 3D orientation test.  Robust.                    */
/*  orient3dslow()   Another exact 3D orientation test.  Robust.             */
/*  orient3d()   Adaptive exact 3D orientation test.  Robust.                */
/*                                                                           */
/*               Return a positive value if the point pd lies below the      */
/*               plane passing through pa, pb, and pc; "below" is defined so */
/*               that pa, pb, and pc appear in counterclockwise order when   */
/*               viewed from above the plane.  Returns a negative value if   */
/*               pd lies above the plane.  Returns zero if the points are    */
/*               coplanar.  The result is also a rough approximation of six  */
/*               times the signed volume of the tetrahedron defined by the   */
/*               four points.                                                */
/*                                                                           */
/*  Only the first and last routine should be used; the middle two are for   */
/*  timings.                                                                 */
/*                                                                           */
/*  The last three use exact arithmetic to ensure a correct answer.  The     */
/*  result returned is the determinant of a matrix.  In orient3d() only,     */
/*  this determinant is computed adaptively, in the sense that exact         */
/*  arithmetic is used only to the degree it is needed to ensure that the    */
/*  returned value has the correct sign.  Hence, orient3d() is usually quite */
/*  fast, but will run more slowly when the input points are coplanar or     */
/*  nearly so.                                                               */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn orient3dfast(mut pa: *const f64,
                                      mut pb: *const f64,
                                      mut pc: *const f64,
                                      mut pd: *const f64)
 -> f64 {
    let mut adx: f64 = 0.;
    let mut bdx: f64 = 0.;
    let mut cdx: f64 = 0.;
    let mut ady: f64 = 0.;
    let mut bdy: f64 = 0.;
    let mut cdy: f64 = 0.;
    let mut adz: f64 = 0.;
    let mut bdz: f64 = 0.;
    let mut cdz: f64 = 0.;
    adx =
        *pa.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    bdx =
        *pb.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    cdx =
        *pc.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    ady =
        *pa.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    bdy =
        *pb.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    cdy =
        *pc.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    adz =
        *pa.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    bdz =
        *pb.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    cdz =
        *pc.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    return adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) +
               cdx * (ady * bdz - adz * bdy);
}

pub unsafe fn orient3dexact(mut pa: *const f64,
                                       mut pb: *const f64,
                                       mut pc: *const f64,
                                       mut pd: *const f64)
 -> f64 {
    let mut axby1: f64 = 0.;
    let mut bxcy1: f64 = 0.;
    let mut cxdy1: f64 = 0.;
    let mut dxay1: f64 = 0.;
    let mut axcy1: f64 = 0.;
    let mut bxdy1: f64 = 0.;
    let mut bxay1: f64 = 0.;
    let mut cxby1: f64 = 0.;
    let mut dxcy1: f64 = 0.;
    let mut axdy1: f64 = 0.;
    let mut cxay1: f64 = 0.;
    let mut dxby1: f64 = 0.;
    let mut axby0: f64 = 0.;
    let mut bxcy0: f64 = 0.;
    let mut cxdy0: f64 = 0.;
    let mut dxay0: f64 = 0.;
    let mut axcy0: f64 = 0.;
    let mut bxdy0: f64 = 0.;
    let mut bxay0: f64 = 0.;
    let mut cxby0: f64 = 0.;
    let mut dxcy0: f64 = 0.;
    let mut axdy0: f64 = 0.;
    let mut cxay0: f64 = 0.;
    let mut dxby0: f64 = 0.;
    let mut ab: [f64; 4] = [0.; 4];
    let mut bc: [f64; 4] = [0.; 4];
    let mut cd: [f64; 4] = [0.; 4];
    let mut da: [f64; 4] = [0.; 4];
    let mut ac: [f64; 4] = [0.; 4];
    let mut bd: [f64; 4] = [0.; 4];
    let mut temp8: [f64; 8] = [0.; 8];
    let mut templen: i32 = 0;
    let mut abc: [f64; 12] = [0.; 12];
    let mut bcd: [f64; 12] = [0.; 12];
    let mut cda: [f64; 12] = [0.; 12];
    let mut dab: [f64; 12] = [0.; 12];
    let mut abclen: i32 = 0;
    let mut bcdlen: i32 = 0;
    let mut cdalen: i32 = 0;
    let mut dablen: i32 = 0;
    let mut adet: [f64; 24] = [0.; 24];
    let mut bdet: [f64; 24] = [0.; 24];
    let mut cdet: [f64; 24] = [0.; 24];
    let mut ddet: [f64; 24] = [0.; 24];
    let mut alen: i32 = 0;
    let mut blen: i32 = 0;
    let mut clen: i32 = 0;
    let mut dlen: i32 = 0;
    let mut abdet: [f64; 48] = [0.; 48];
    let mut cddet: [f64; 48] = [0.; 48];
    let mut ablen: i32 = 0;
    let mut cdlen: i32 = 0;
    let mut deter: [f64; 96] = [0.; 96];
    let mut deterlen: i32 = 0;
    let mut i: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    Two_Product(*pa.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut axby1,
                &mut axby0);
    Two_Product(*pb.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut bxay1,
                &mut bxay0);
    Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab.as_mut_ptr());
    Two_Product(*pb.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut bxcy1,
                &mut bxcy0);
    Two_Product(*pc.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut cxby1,
                &mut cxby0);
    Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc.as_mut_ptr());
    Two_Product(*pc.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut cxdy1,
                &mut cxdy0);
    Two_Product(*pd.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut dxcy1,
                &mut dxcy0);
    Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd.as_mut_ptr());
    Two_Product(*pd.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut dxay1,
                &mut dxay0);
    Two_Product(*pa.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut axdy1,
                &mut axdy0);
    Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da.as_mut_ptr());
    Two_Product(*pa.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut axcy1,
                &mut axcy0);
    Two_Product(*pc.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut cxay1,
                &mut cxay0);
    Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac.as_mut_ptr());
    Two_Product(*pb.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut bxdy1,
                &mut bxdy0);
    Two_Product(*pd.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut dxby1,
                &mut dxby0);
    Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd.as_mut_ptr());
    templen =
        fast_expansion_sum_zeroelim(4 as i32, cd.as_mut_ptr(),
                                    4 as i32, da.as_mut_ptr(),
                                    temp8.as_mut_ptr());
    cdalen =
        fast_expansion_sum_zeroelim(templen, temp8.as_mut_ptr(),
                                    4 as i32, ac.as_mut_ptr(),
                                    cda.as_mut_ptr());
    templen =
        fast_expansion_sum_zeroelim(4 as i32, da.as_mut_ptr(),
                                    4 as i32, ab.as_mut_ptr(),
                                    temp8.as_mut_ptr());
    dablen =
        fast_expansion_sum_zeroelim(templen, temp8.as_mut_ptr(),
                                    4 as i32, bd.as_mut_ptr(),
                                    dab.as_mut_ptr());
    i = 0 as i32;
    while i < 4 as i32 {
        bd[i as usize] = -bd[i as usize];
        ac[i as usize] = -ac[i as usize];
        i += 1
    }
    templen =
        fast_expansion_sum_zeroelim(4 as i32, ab.as_mut_ptr(),
                                    4 as i32, bc.as_mut_ptr(),
                                    temp8.as_mut_ptr());
    abclen =
        fast_expansion_sum_zeroelim(templen, temp8.as_mut_ptr(),
                                    4 as i32, ac.as_mut_ptr(),
                                    abc.as_mut_ptr());
    templen =
        fast_expansion_sum_zeroelim(4 as i32, bc.as_mut_ptr(),
                                    4 as i32, cd.as_mut_ptr(),
                                    temp8.as_mut_ptr());
    bcdlen =
        fast_expansion_sum_zeroelim(templen, temp8.as_mut_ptr(),
                                    4 as i32, bd.as_mut_ptr(),
                                    bcd.as_mut_ptr());
    alen =
        scale_expansion_zeroelim(bcdlen, bcd.as_mut_ptr(),
                                 *pa.offset(2 as i32 as isize),
                                 adet.as_mut_ptr());
    blen =
        scale_expansion_zeroelim(cdalen, cda.as_mut_ptr(),
                                 -*pb.offset(2 as i32 as isize),
                                 bdet.as_mut_ptr());
    clen =
        scale_expansion_zeroelim(dablen, dab.as_mut_ptr(),
                                 *pc.offset(2 as i32 as isize),
                                 cdet.as_mut_ptr());
    dlen =
        scale_expansion_zeroelim(abclen, abc.as_mut_ptr(),
                                 -*pd.offset(2 as i32 as isize),
                                 ddet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    cdlen =
        fast_expansion_sum_zeroelim(clen, cdet.as_mut_ptr(), dlen,
                                    ddet.as_mut_ptr(), cddet.as_mut_ptr());
    deterlen =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), cdlen,
                                    cddet.as_mut_ptr(), deter.as_mut_ptr());
    return deter[(deterlen - 1 as i32) as usize];
}

pub unsafe fn orient3dslow(mut pa: *const f64,
                                      mut pb: *const f64,
                                      mut pc: *const f64,
                                      mut pd: *const f64)
 -> f64 {
    let mut adx: f64 = 0.;
    let mut ady: f64 = 0.;
    let mut adz: f64 = 0.;
    let mut bdx: f64 = 0.;
    let mut bdy: f64 = 0.;
    let mut bdz: f64 = 0.;
    let mut cdx: f64 = 0.;
    let mut cdy: f64 = 0.;
    let mut cdz: f64 = 0.;
    let mut adxtail: f64 = 0.;
    let mut adytail: f64 = 0.;
    let mut adztail: f64 = 0.;
    let mut bdxtail: f64 = 0.;
    let mut bdytail: f64 = 0.;
    let mut bdztail: f64 = 0.;
    let mut cdxtail: f64 = 0.;
    let mut cdytail: f64 = 0.;
    let mut cdztail: f64 = 0.;
    let mut negate: f64 = 0.;
    let mut negatetail: f64 = 0.;
    let mut axby7: f64 = 0.;
    let mut bxcy7: f64 = 0.;
    let mut axcy7: f64 = 0.;
    let mut bxay7: f64 = 0.;
    let mut cxby7: f64 = 0.;
    let mut cxay7: f64 = 0.;
    let mut axby: [f64; 8] = [0.; 8];
    let mut bxcy: [f64; 8] = [0.; 8];
    let mut axcy: [f64; 8] = [0.; 8];
    let mut bxay: [f64; 8] = [0.; 8];
    let mut cxby: [f64; 8] = [0.; 8];
    let mut cxay: [f64; 8] = [0.; 8];
    let mut temp16: [f64; 16] = [0.; 16];
    let mut temp32: [f64; 32] = [0.; 32];
    let mut temp32t: [f64; 32] = [0.; 32];
    let mut temp16len: i32 = 0;
    let mut temp32len: i32 = 0;
    let mut temp32tlen: i32 = 0;
    let mut adet: [f64; 64] = [0.; 64];
    let mut bdet: [f64; 64] = [0.; 64];
    let mut cdet: [f64; 64] = [0.; 64];
    let mut alen: i32 = 0;
    let mut blen: i32 = 0;
    let mut clen: i32 = 0;
    let mut abdet: [f64; 128] = [0.; 128];
    let mut ablen: i32 = 0;
    let mut deter: [f64; 192] = [0.; 192];
    let mut deterlen: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut a0hi: f64 = 0.;
    let mut a0lo: f64 = 0.;
    let mut a1hi: f64 = 0.;
    let mut a1lo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _k: f64 = 0.;
    let mut _l: f64 = 0.;
    let mut _m: f64 = 0.;
    let mut _n: f64 = 0.;
    let mut _0: f64 = 0.;
    let mut _1: f64 = 0.;
    let mut _2: f64 = 0.;
    Two_Diff(*pa.offset(0 as i32 as isize),
             *pd.offset(0 as i32 as isize), &mut adx, &mut adxtail);
    Two_Diff(*pa.offset(1 as i32 as isize),
             *pd.offset(1 as i32 as isize), &mut ady, &mut adytail);
    Two_Diff(*pa.offset(2 as i32 as isize),
             *pd.offset(2 as i32 as isize), &mut adz, &mut adztail);
    Two_Diff(*pb.offset(0 as i32 as isize),
             *pd.offset(0 as i32 as isize), &mut bdx, &mut bdxtail);
    Two_Diff(*pb.offset(1 as i32 as isize),
             *pd.offset(1 as i32 as isize), &mut bdy, &mut bdytail);
    Two_Diff(*pb.offset(2 as i32 as isize),
             *pd.offset(2 as i32 as isize), &mut bdz, &mut bdztail);
    Two_Diff(*pc.offset(0 as i32 as isize),
             *pd.offset(0 as i32 as isize), &mut cdx, &mut cdxtail);
    Two_Diff(*pc.offset(1 as i32 as isize),
             *pd.offset(1 as i32 as isize), &mut cdy, &mut cdytail);
    Two_Diff(*pc.offset(2 as i32 as isize),
             *pd.offset(2 as i32 as isize), &mut cdz, &mut cdztail);
    Two_Two_Product(adx, adxtail, bdy, bdytail, axby.as_mut_ptr());
    negate = -ady;
    negatetail = -adytail;
    Two_Two_Product(bdx, bdxtail, negate, negatetail, bxay.as_mut_ptr());
    Two_Two_Product(bdx, bdxtail, cdy, cdytail, bxcy.as_mut_ptr());
    negate = -bdy;
    negatetail = -bdytail;
    Two_Two_Product(cdx, cdxtail, negate, negatetail, cxby.as_mut_ptr());
    Two_Two_Product(cdx, cdxtail, ady, adytail, cxay.as_mut_ptr());
    negate = -cdy;
    negatetail = -cdytail;
    Two_Two_Product(adx, adxtail, negate, negatetail, axcy.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(8 as i32, bxcy.as_mut_ptr(),
                                    8 as i32, cxby.as_mut_ptr(),
                                    temp16.as_mut_ptr());
    temp32len =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), adz,
                                 temp32.as_mut_ptr());
    temp32tlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), adztail,
                                 temp32t.as_mut_ptr());
    alen =
        fast_expansion_sum_zeroelim(temp32len, temp32.as_mut_ptr(),
                                    temp32tlen, temp32t.as_mut_ptr(),
                                    adet.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(8 as i32, cxay.as_mut_ptr(),
                                    8 as i32, axcy.as_mut_ptr(),
                                    temp16.as_mut_ptr());
    temp32len =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), bdz,
                                 temp32.as_mut_ptr());
    temp32tlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), bdztail,
                                 temp32t.as_mut_ptr());
    blen =
        fast_expansion_sum_zeroelim(temp32len, temp32.as_mut_ptr(),
                                    temp32tlen, temp32t.as_mut_ptr(),
                                    bdet.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(8 as i32, axby.as_mut_ptr(),
                                    8 as i32, bxay.as_mut_ptr(),
                                    temp16.as_mut_ptr());
    temp32len =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), cdz,
                                 temp32.as_mut_ptr());
    temp32tlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), cdztail,
                                 temp32t.as_mut_ptr());
    clen =
        fast_expansion_sum_zeroelim(temp32len, temp32.as_mut_ptr(),
                                    temp32tlen, temp32t.as_mut_ptr(),
                                    cdet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    deterlen =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), clen,
                                    cdet.as_mut_ptr(), deter.as_mut_ptr());
    return deter[(deterlen - 1 as i32) as usize];
}

pub unsafe fn orient3dadapt(mut pa: *const f64,
                                       mut pb: *const f64,
                                       mut pc: *const f64,
                                       mut pd: *const f64,
                                       mut permanent: f64)
 -> f64 {
    let mut adx: f64 = 0.;
    let mut bdx: f64 = 0.;
    let mut cdx: f64 = 0.;
    let mut ady: f64 = 0.;
    let mut bdy: f64 = 0.;
    let mut cdy: f64 = 0.;
    let mut adz: f64 = 0.;
    let mut bdz: f64 = 0.;
    let mut cdz: f64 = 0.;
    let mut det: f64 = 0.;
    let mut errbound: f64 = 0.;
    let mut bdxcdy1: f64 = 0.;
    let mut cdxbdy1: f64 = 0.;
    let mut cdxady1: f64 = 0.;
    let mut adxcdy1: f64 = 0.;
    let mut adxbdy1: f64 = 0.;
    let mut bdxady1: f64 = 0.;
    let mut bdxcdy0: f64 = 0.;
    let mut cdxbdy0: f64 = 0.;
    let mut cdxady0: f64 = 0.;
    let mut adxcdy0: f64 = 0.;
    let mut adxbdy0: f64 = 0.;
    let mut bdxady0: f64 = 0.;
    let mut bc: [f64; 4] = [0.; 4];
    let mut ca: [f64; 4] = [0.; 4];
    let mut ab: [f64; 4] = [0.; 4];
    let mut bc3: f64 = 0.;
    let mut ca3: f64 = 0.;
    let mut ab3: f64 = 0.;
    let mut adet: [f64; 8] = [0.; 8];
    let mut bdet: [f64; 8] = [0.; 8];
    let mut cdet: [f64; 8] = [0.; 8];
    let mut alen: i32 = 0;
    let mut blen: i32 = 0;
    let mut clen: i32 = 0;
    let mut abdet: [f64; 16] = [0.; 16];
    let mut ablen: i32 = 0;
    let mut finnow: *mut f64 = 0 as *mut f64;
    let mut finother: *mut f64 = 0 as *mut f64;
    let mut finswap: *mut f64 = 0 as *mut f64;
    let mut fin1: [f64; 192] = [0.; 192];
    let mut fin2: [f64; 192] = [0.; 192];
    let mut finlength: i32 = 0;
    let mut adxtail: f64 = 0.;
    let mut bdxtail: f64 = 0.;
    let mut cdxtail: f64 = 0.;
    let mut adytail: f64 = 0.;
    let mut bdytail: f64 = 0.;
    let mut cdytail: f64 = 0.;
    let mut adztail: f64 = 0.;
    let mut bdztail: f64 = 0.;
    let mut cdztail: f64 = 0.;
    let mut at_blarge: f64 = 0.;
    let mut at_clarge: f64 = 0.;
    let mut bt_clarge: f64 = 0.;
    let mut bt_alarge: f64 = 0.;
    let mut ct_alarge: f64 = 0.;
    let mut ct_blarge: f64 = 0.;
    let mut at_b: [f64; 4] = [0.; 4];
    let mut at_c: [f64; 4] = [0.; 4];
    let mut bt_c: [f64; 4] = [0.; 4];
    let mut bt_a: [f64; 4] = [0.; 4];
    let mut ct_a: [f64; 4] = [0.; 4];
    let mut ct_b: [f64; 4] = [0.; 4];
    let mut at_blen: i32 = 0;
    let mut at_clen: i32 = 0;
    let mut bt_clen: i32 = 0;
    let mut bt_alen: i32 = 0;
    let mut ct_alen: i32 = 0;
    let mut ct_blen: i32 = 0;
    let mut bdxt_cdy1: f64 = 0.;
    let mut cdxt_bdy1: f64 = 0.;
    let mut cdxt_ady1: f64 = 0.;
    let mut adxt_cdy1: f64 = 0.;
    let mut adxt_bdy1: f64 = 0.;
    let mut bdxt_ady1: f64 = 0.;
    let mut bdxt_cdy0: f64 = 0.;
    let mut cdxt_bdy0: f64 = 0.;
    let mut cdxt_ady0: f64 = 0.;
    let mut adxt_cdy0: f64 = 0.;
    let mut adxt_bdy0: f64 = 0.;
    let mut bdxt_ady0: f64 = 0.;
    let mut bdyt_cdx1: f64 = 0.;
    let mut cdyt_bdx1: f64 = 0.;
    let mut cdyt_adx1: f64 = 0.;
    let mut adyt_cdx1: f64 = 0.;
    let mut adyt_bdx1: f64 = 0.;
    let mut bdyt_adx1: f64 = 0.;
    let mut bdyt_cdx0: f64 = 0.;
    let mut cdyt_bdx0: f64 = 0.;
    let mut cdyt_adx0: f64 = 0.;
    let mut adyt_cdx0: f64 = 0.;
    let mut adyt_bdx0: f64 = 0.;
    let mut bdyt_adx0: f64 = 0.;
    let mut bct: [f64; 8] = [0.; 8];
    let mut cat: [f64; 8] = [0.; 8];
    let mut abt: [f64; 8] = [0.; 8];
    let mut bctlen: i32 = 0;
    let mut catlen: i32 = 0;
    let mut abtlen: i32 = 0;
    let mut bdxt_cdyt1: f64 = 0.;
    let mut cdxt_bdyt1: f64 = 0.;
    let mut cdxt_adyt1: f64 = 0.;
    let mut adxt_cdyt1: f64 = 0.;
    let mut adxt_bdyt1: f64 = 0.;
    let mut bdxt_adyt1: f64 = 0.;
    let mut bdxt_cdyt0: f64 = 0.;
    let mut cdxt_bdyt0: f64 = 0.;
    let mut cdxt_adyt0: f64 = 0.;
    let mut adxt_cdyt0: f64 = 0.;
    let mut adxt_bdyt0: f64 = 0.;
    let mut bdxt_adyt0: f64 = 0.;
    let mut u: [f64; 4] = [0.; 4];
    let mut v: [f64; 12] = [0.; 12];
    let mut w: [f64; 16] = [0.; 16];
    let mut u3: f64 = 0.;
    let mut vlength: i32 = 0;
    let mut wlength: i32 = 0;
    let mut negate: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _k: f64 = 0.;
    let mut _0: f64 = 0.;
    adx =
        *pa.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    bdx =
        *pb.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    cdx =
        *pc.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    ady =
        *pa.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    bdy =
        *pb.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    cdy =
        *pc.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    adz =
        *pa.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    bdz =
        *pb.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    cdz =
        *pc.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    Two_Product(bdx, cdy, &mut bdxcdy1, &mut bdxcdy0);
    Two_Product(cdx, bdy, &mut cdxbdy1, &mut cdxbdy0);
    Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc.as_mut_ptr());
    alen =
        scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(), adz,
                                 adet.as_mut_ptr());
    Two_Product(cdx, ady, &mut cdxady1, &mut cdxady0);
    Two_Product(adx, cdy, &mut adxcdy1, &mut adxcdy0);
    Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca.as_mut_ptr());
    blen =
        scale_expansion_zeroelim(4 as i32, ca.as_mut_ptr(), bdz,
                                 bdet.as_mut_ptr());
    Two_Product(adx, bdy, &mut adxbdy1, &mut adxbdy0);
    Two_Product(bdx, ady, &mut bdxady1, &mut bdxady0);
    Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab.as_mut_ptr());
    clen =
        scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(), cdz,
                                 cdet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    finlength =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), clen,
                                    cdet.as_mut_ptr(), fin1.as_mut_ptr());
    det = estimate(finlength, fin1.as_mut_ptr());
    errbound = PARAMS.o3derrboundB * permanent;
    if det >= errbound || -det >= errbound { return det }
    adxtail =
        Two_Diff_Tail(*pa.offset(0 as i32 as isize),
                      *pd.offset(0 as i32 as isize), adx);
    bdxtail =
        Two_Diff_Tail(*pb.offset(0 as i32 as isize),
                      *pd.offset(0 as i32 as isize), bdx);
    cdxtail =
        Two_Diff_Tail(*pc.offset(0 as i32 as isize),
                      *pd.offset(0 as i32 as isize), cdx);
    adytail =
        Two_Diff_Tail(*pa.offset(1 as i32 as isize),
                      *pd.offset(1 as i32 as isize), ady);
    bdytail =
        Two_Diff_Tail(*pb.offset(1 as i32 as isize),
                      *pd.offset(1 as i32 as isize), bdy);
    cdytail =
        Two_Diff_Tail(*pc.offset(1 as i32 as isize),
                      *pd.offset(1 as i32 as isize), cdy);
    adztail =
        Two_Diff_Tail(*pa.offset(2 as i32 as isize),
                      *pd.offset(2 as i32 as isize), adz);
    bdztail =
        Two_Diff_Tail(*pb.offset(2 as i32 as isize),
                      *pd.offset(2 as i32 as isize), bdz);
    cdztail =
        Two_Diff_Tail(*pc.offset(2 as i32 as isize),
                      *pd.offset(2 as i32 as isize), cdz);
    if adxtail == 0.0f64 && bdxtail == 0.0f64 && cdxtail == 0.0f64 &&
           adytail == 0.0f64 && bdytail == 0.0f64 && cdytail == 0.0f64 &&
           adztail == 0.0f64 && bdztail == 0.0f64 && cdztail == 0.0f64 {
        return det
    }
    errbound = PARAMS.o3derrboundC * permanent + PARAMS.resulterrbound * Absolute(det);
    det +=
        adz *
            (bdx * cdytail + cdy * bdxtail - (bdy * cdxtail + cdx * bdytail))
            + adztail * (bdx * cdy - bdy * cdx) +
            (bdz *
                 (cdx * adytail + ady * cdxtail -
                      (cdy * adxtail + adx * cdytail)) +
                 bdztail * (cdx * ady - cdy * adx)) +
            (cdz *
                 (adx * bdytail + bdy * adxtail -
                      (ady * bdxtail + bdx * adytail)) +
                 cdztail * (adx * bdy - ady * bdx));
    if det >= errbound || -det >= errbound { return det }
    finnow = fin1.as_mut_ptr();
    finother = fin2.as_mut_ptr();
    if adxtail == 0.0f64 {
        if adytail == 0.0f64 {
            at_b[0 as i32 as usize] = 0.0f64;
            at_blen = 1 as i32;
            at_c[0 as i32 as usize] = 0.0f64;
            at_clen = 1 as i32
        } else {
            negate = -adytail;
            Two_Product(negate, bdx, &mut at_blarge,
                        &mut *at_b.as_mut_ptr().offset(0 as i32 as
                                                           isize));
            at_b[1 as i32 as usize] = at_blarge;
            at_blen = 2 as i32;
            Two_Product(adytail, cdx, &mut at_clarge,
                        &mut *at_c.as_mut_ptr().offset(0 as i32 as
                                                           isize));
            at_c[1 as i32 as usize] = at_clarge;
            at_clen = 2 as i32
        }
    } else if adytail == 0.0f64 {
        Two_Product(adxtail, bdy, &mut at_blarge,
                    &mut *at_b.as_mut_ptr().offset(0 as i32 as
                                                       isize));
        at_b[1 as i32 as usize] = at_blarge;
        at_blen = 2 as i32;
        negate = -adxtail;
        Two_Product(negate, cdy, &mut at_clarge,
                    &mut *at_c.as_mut_ptr().offset(0 as i32 as
                                                       isize));
        at_c[1 as i32 as usize] = at_clarge;
        at_clen = 2 as i32
    } else {
        Two_Product(adxtail, bdy, &mut adxt_bdy1, &mut adxt_bdy0);
        Two_Product(adytail, bdx, &mut adyt_bdx1, &mut adyt_bdx0);
        Two_Two_Diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0,
                     at_b.as_mut_ptr());
        at_blen = 4 as i32;
        Two_Product(adytail, cdx, &mut adyt_cdx1, &mut adyt_cdx0);
        Two_Product(adxtail, cdy, &mut adxt_cdy1, &mut adxt_cdy0);
        Two_Two_Diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0,
                     at_c.as_mut_ptr());
        at_clen = 4 as i32
    }
    if bdxtail == 0.0f64 {
        if bdytail == 0.0f64 {
            bt_c[0 as i32 as usize] = 0.0f64;
            bt_clen = 1 as i32;
            bt_a[0 as i32 as usize] = 0.0f64;
            bt_alen = 1 as i32
        } else {
            negate = -bdytail;
            Two_Product(negate, cdx, &mut bt_clarge,
                        &mut *bt_c.as_mut_ptr().offset(0 as i32 as
                                                           isize));
            bt_c[1 as i32 as usize] = bt_clarge;
            bt_clen = 2 as i32;
            Two_Product(bdytail, adx, &mut bt_alarge,
                        &mut *bt_a.as_mut_ptr().offset(0 as i32 as
                                                           isize));
            bt_a[1 as i32 as usize] = bt_alarge;
            bt_alen = 2 as i32
        }
    } else if bdytail == 0.0f64 {
        Two_Product(bdxtail, cdy, &mut bt_clarge,
                    &mut *bt_c.as_mut_ptr().offset(0 as i32 as
                                                       isize));
        bt_c[1 as i32 as usize] = bt_clarge;
        bt_clen = 2 as i32;
        negate = -bdxtail;
        Two_Product(negate, ady, &mut bt_alarge,
                    &mut *bt_a.as_mut_ptr().offset(0 as i32 as
                                                       isize));
        bt_a[1 as i32 as usize] = bt_alarge;
        bt_alen = 2 as i32
    } else {
        Two_Product(bdxtail, cdy, &mut bdxt_cdy1, &mut bdxt_cdy0);
        Two_Product(bdytail, cdx, &mut bdyt_cdx1, &mut bdyt_cdx0);
        Two_Two_Diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0,
                     bt_c.as_mut_ptr());
        bt_clen = 4 as i32;
        Two_Product(bdytail, adx, &mut bdyt_adx1, &mut bdyt_adx0);
        Two_Product(bdxtail, ady, &mut bdxt_ady1, &mut bdxt_ady0);
        Two_Two_Diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0,
                     bt_a.as_mut_ptr());
        bt_alen = 4 as i32
    }
    if cdxtail == 0.0f64 {
        if cdytail == 0.0f64 {
            ct_a[0 as i32 as usize] = 0.0f64;
            ct_alen = 1 as i32;
            ct_b[0 as i32 as usize] = 0.0f64;
            ct_blen = 1 as i32
        } else {
            negate = -cdytail;
            Two_Product(negate, adx, &mut ct_alarge,
                        &mut *ct_a.as_mut_ptr().offset(0 as i32 as
                                                           isize));
            ct_a[1 as i32 as usize] = ct_alarge;
            ct_alen = 2 as i32;
            Two_Product(cdytail, bdx, &mut ct_blarge,
                        &mut *ct_b.as_mut_ptr().offset(0 as i32 as
                                                           isize));
            ct_b[1 as i32 as usize] = ct_blarge;
            ct_blen = 2 as i32
        }
    } else if cdytail == 0.0f64 {
        Two_Product(cdxtail, ady, &mut ct_alarge,
                    &mut *ct_a.as_mut_ptr().offset(0 as i32 as
                                                       isize));
        ct_a[1 as i32 as usize] = ct_alarge;
        ct_alen = 2 as i32;
        negate = -cdxtail;
        Two_Product(negate, bdy, &mut ct_blarge,
                    &mut *ct_b.as_mut_ptr().offset(0 as i32 as
                                                       isize));
        ct_b[1 as i32 as usize] = ct_blarge;
        ct_blen = 2 as i32
    } else {
        Two_Product(cdxtail, ady, &mut cdxt_ady1, &mut cdxt_ady0);
        Two_Product(cdytail, adx, &mut cdyt_adx1, &mut cdyt_adx0);
        Two_Two_Diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0,
                     ct_a.as_mut_ptr());
        ct_alen = 4 as i32;
        Two_Product(cdytail, bdx, &mut cdyt_bdx1, &mut cdyt_bdx0);
        Two_Product(cdxtail, bdy, &mut cdxt_bdy1, &mut cdxt_bdy0);
        Two_Two_Diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0,
                     ct_b.as_mut_ptr());
        ct_blen = 4 as i32
    }
    bctlen =
        fast_expansion_sum_zeroelim(bt_clen, bt_c.as_mut_ptr(), ct_blen,
                                    ct_b.as_mut_ptr(), bct.as_mut_ptr());
    wlength =
        scale_expansion_zeroelim(bctlen, bct.as_mut_ptr(), adz,
                                 w.as_mut_ptr());
    finlength =
        fast_expansion_sum_zeroelim(finlength, finnow, wlength,
                                    w.as_mut_ptr(), finother);
    finswap = finnow;
    finnow = finother;
    finother = finswap;
    catlen =
        fast_expansion_sum_zeroelim(ct_alen, ct_a.as_mut_ptr(), at_clen,
                                    at_c.as_mut_ptr(), cat.as_mut_ptr());
    wlength =
        scale_expansion_zeroelim(catlen, cat.as_mut_ptr(), bdz,
                                 w.as_mut_ptr());
    finlength =
        fast_expansion_sum_zeroelim(finlength, finnow, wlength,
                                    w.as_mut_ptr(), finother);
    finswap = finnow;
    finnow = finother;
    finother = finswap;
    abtlen =
        fast_expansion_sum_zeroelim(at_blen, at_b.as_mut_ptr(), bt_alen,
                                    bt_a.as_mut_ptr(), abt.as_mut_ptr());
    wlength =
        scale_expansion_zeroelim(abtlen, abt.as_mut_ptr(), cdz,
                                 w.as_mut_ptr());
    finlength =
        fast_expansion_sum_zeroelim(finlength, finnow, wlength,
                                    w.as_mut_ptr(), finother);
    finswap = finnow;
    finnow = finother;
    finother = finswap;
    if adztail != 0.0f64 {
        vlength =
            scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(),
                                     adztail, v.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, vlength,
                                        v.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if bdztail != 0.0f64 {
        vlength =
            scale_expansion_zeroelim(4 as i32, ca.as_mut_ptr(),
                                     bdztail, v.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, vlength,
                                        v.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if cdztail != 0.0f64 {
        vlength =
            scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(),
                                     cdztail, v.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, vlength,
                                        v.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if adxtail != 0.0f64 {
        if bdytail != 0.0f64 {
            Two_Product(adxtail, bdytail, &mut adxt_bdyt1, &mut adxt_bdyt0);
            Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdz, u.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow,
                                            4 as i32, u.as_mut_ptr(),
                                            finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if cdztail != 0.0f64 {
                Two_One_Product(adxt_bdyt1, adxt_bdyt0, cdztail,
                                u.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow,
                                                4 as i32,
                                                u.as_mut_ptr(), finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
        }
        if cdytail != 0.0f64 {
            negate = -adxtail;
            Two_Product(negate, cdytail, &mut adxt_cdyt1, &mut adxt_cdyt0);
            Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdz, u.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow,
                                            4 as i32, u.as_mut_ptr(),
                                            finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if bdztail != 0.0f64 {
                Two_One_Product(adxt_cdyt1, adxt_cdyt0, bdztail,
                                u.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow,
                                                4 as i32,
                                                u.as_mut_ptr(), finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
        }
    }
    if bdxtail != 0.0f64 {
        if cdytail != 0.0f64 {
            Two_Product(bdxtail, cdytail, &mut bdxt_cdyt1, &mut bdxt_cdyt0);
            Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adz, u.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow,
                                            4 as i32, u.as_mut_ptr(),
                                            finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if adztail != 0.0f64 {
                Two_One_Product(bdxt_cdyt1, bdxt_cdyt0, adztail,
                                u.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow,
                                                4 as i32,
                                                u.as_mut_ptr(), finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
        }
        if adytail != 0.0f64 {
            negate = -bdxtail;
            Two_Product(negate, adytail, &mut bdxt_adyt1, &mut bdxt_adyt0);
            Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdz, u.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow,
                                            4 as i32, u.as_mut_ptr(),
                                            finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if cdztail != 0.0f64 {
                Two_One_Product(bdxt_adyt1, bdxt_adyt0, cdztail,
                                u.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow,
                                                4 as i32,
                                                u.as_mut_ptr(), finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
        }
    }
    if cdxtail != 0.0f64 {
        if adytail != 0.0f64 {
            Two_Product(cdxtail, adytail, &mut cdxt_adyt1, &mut cdxt_adyt0);
            Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdz, u.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow,
                                            4 as i32, u.as_mut_ptr(),
                                            finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if bdztail != 0.0f64 {
                Two_One_Product(cdxt_adyt1, cdxt_adyt0, bdztail,
                                u.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow,
                                                4 as i32,
                                                u.as_mut_ptr(), finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
        }
        if bdytail != 0.0f64 {
            negate = -cdxtail;
            Two_Product(negate, bdytail, &mut cdxt_bdyt1, &mut cdxt_bdyt0);
            Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adz, u.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow,
                                            4 as i32, u.as_mut_ptr(),
                                            finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if adztail != 0.0f64 {
                Two_One_Product(cdxt_bdyt1, cdxt_bdyt0, adztail,
                                u.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow,
                                                4 as i32,
                                                u.as_mut_ptr(), finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
        }
    }
    if adztail != 0.0f64 {
        wlength =
            scale_expansion_zeroelim(bctlen, bct.as_mut_ptr(), adztail,
                                     w.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, wlength,
                                        w.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if bdztail != 0.0f64 {
        wlength =
            scale_expansion_zeroelim(catlen, cat.as_mut_ptr(), bdztail,
                                     w.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, wlength,
                                        w.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if cdztail != 0.0f64 {
        wlength =
            scale_expansion_zeroelim(abtlen, abt.as_mut_ptr(), cdztail,
                                     w.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, wlength,
                                        w.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    return *finnow.offset((finlength - 1 as i32) as isize);
}

pub unsafe fn orient3d(mut pa: *const f64,
                                  mut pb: *const f64,
                                  mut pc: *const f64,
                                  mut pd: *const f64)
 -> f64 {
    let mut adx: f64 = 0.;
    let mut bdx: f64 = 0.;
    let mut cdx: f64 = 0.;
    let mut ady: f64 = 0.;
    let mut bdy: f64 = 0.;
    let mut cdy: f64 = 0.;
    let mut adz: f64 = 0.;
    let mut bdz: f64 = 0.;
    let mut cdz: f64 = 0.;
    let mut bdxcdy: f64 = 0.;
    let mut cdxbdy: f64 = 0.;
    let mut cdxady: f64 = 0.;
    let mut adxcdy: f64 = 0.;
    let mut adxbdy: f64 = 0.;
    let mut bdxady: f64 = 0.;
    let mut det: f64 = 0.;
    let mut permanent: f64 = 0.;
    let mut errbound: f64 = 0.;
    adx =
        *pa.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    bdx =
        *pb.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    cdx =
        *pc.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    ady =
        *pa.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    bdy =
        *pb.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    cdy =
        *pc.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    adz =
        *pa.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    bdz =
        *pb.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    cdz =
        *pc.offset(2 as i32 as isize) -
            *pd.offset(2 as i32 as isize);
    bdxcdy = bdx * cdy;
    cdxbdy = cdx * bdy;
    cdxady = cdx * ady;
    adxcdy = adx * cdy;
    adxbdy = adx * bdy;
    bdxady = bdx * ady;
    det =
        adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) +
            cdz * (adxbdy - bdxady);
    permanent =
        (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(adz) +
            (Absolute(cdxady) + Absolute(adxcdy)) * Absolute(bdz) +
            (Absolute(adxbdy) + Absolute(bdxady)) * Absolute(cdz);
    errbound = PARAMS.o3derrboundA * permanent;
    if det > errbound || -det > errbound { return det }
    return orient3dadapt(pa, pb, pc, pd, permanent);
}
/* ****************************************************************************/
/*                                                                           */
/*  incirclefast()   Approximate 2D incircle test.  Nonrobust.               */
/*  incircleexact()   Exact 2D incircle test.  Robust.                       */
/*  incircleslow()   Another exact 2D incircle test.  Robust.                */
/*  incircle()   Adaptive exact 2D incircle test.  Robust.                   */
/*                                                                           */
/*               Return a positive value if the point pd lies inside the     */
/*               circle passing through pa, pb, and pc; a negative value if  */
/*               it lies outside; and zero if the four points are cocircular.*/
/*               The points pa, pb, and pc must be in counterclockwise       */
/*               order, or the sign of the result will be reversed.          */
/*                                                                           */
/*  Only the first and last routine should be used; the middle two are for   */
/*  timings.                                                                 */
/*                                                                           */
/*  The last three use exact arithmetic to ensure a correct answer.  The     */
/*  result returned is the determinant of a matrix.  In incircle() only,     */
/*  this determinant is computed adaptively, in the sense that exact         */
/*  arithmetic is used only to the degree it is needed to ensure that the    */
/*  returned value has the correct sign.  Hence, incircle() is usually quite */
/*  fast, but will run more slowly when the input points are cocircular or   */
/*  nearly so.                                                               */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn incirclefast(mut pa: *const f64,
                                      mut pb: *const f64,
                                      mut pc: *const f64,
                                      mut pd: *const f64)
 -> f64 {
    let mut adx: f64 = 0.;
    let mut ady: f64 = 0.;
    let mut bdx: f64 = 0.;
    let mut bdy: f64 = 0.;
    let mut cdx: f64 = 0.;
    let mut cdy: f64 = 0.;
    let mut abdet: f64 = 0.;
    let mut bcdet: f64 = 0.;
    let mut cadet: f64 = 0.;
    let mut alift: f64 = 0.;
    let mut blift: f64 = 0.;
    let mut clift: f64 = 0.;
    adx =
        *pa.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    ady =
        *pa.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    bdx =
        *pb.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    bdy =
        *pb.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    cdx =
        *pc.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    cdy =
        *pc.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    abdet = adx * bdy - bdx * ady;
    bcdet = bdx * cdy - cdx * bdy;
    cadet = cdx * ady - adx * cdy;
    alift = adx * adx + ady * ady;
    blift = bdx * bdx + bdy * bdy;
    clift = cdx * cdx + cdy * cdy;
    return alift * bcdet + blift * cadet + clift * abdet;
}

pub unsafe fn incircleexact(mut pa: *mut f64,
                                       mut pb: *mut f64,
                                       mut pc: *mut f64,
                                       mut pd: *mut f64)
 -> f64 {
    let mut axby1: f64 = 0.;
    let mut bxcy1: f64 = 0.;
    let mut cxdy1: f64 = 0.;
    let mut dxay1: f64 = 0.;
    let mut axcy1: f64 = 0.;
    let mut bxdy1: f64 = 0.;
    let mut bxay1: f64 = 0.;
    let mut cxby1: f64 = 0.;
    let mut dxcy1: f64 = 0.;
    let mut axdy1: f64 = 0.;
    let mut cxay1: f64 = 0.;
    let mut dxby1: f64 = 0.;
    let mut axby0: f64 = 0.;
    let mut bxcy0: f64 = 0.;
    let mut cxdy0: f64 = 0.;
    let mut dxay0: f64 = 0.;
    let mut axcy0: f64 = 0.;
    let mut bxdy0: f64 = 0.;
    let mut bxay0: f64 = 0.;
    let mut cxby0: f64 = 0.;
    let mut dxcy0: f64 = 0.;
    let mut axdy0: f64 = 0.;
    let mut cxay0: f64 = 0.;
    let mut dxby0: f64 = 0.;
    let mut ab: [f64; 4] = [0.; 4];
    let mut bc: [f64; 4] = [0.; 4];
    let mut cd: [f64; 4] = [0.; 4];
    let mut da: [f64; 4] = [0.; 4];
    let mut ac: [f64; 4] = [0.; 4];
    let mut bd: [f64; 4] = [0.; 4];
    let mut temp8: [f64; 8] = [0.; 8];
    let mut templen: i32 = 0;
    let mut abc: [f64; 12] = [0.; 12];
    let mut bcd: [f64; 12] = [0.; 12];
    let mut cda: [f64; 12] = [0.; 12];
    let mut dab: [f64; 12] = [0.; 12];
    let mut abclen: i32 = 0;
    let mut bcdlen: i32 = 0;
    let mut cdalen: i32 = 0;
    let mut dablen: i32 = 0;
    let mut det24x: [f64; 24] = [0.; 24];
    let mut det24y: [f64; 24] = [0.; 24];
    let mut det48x: [f64; 48] = [0.; 48];
    let mut det48y: [f64; 48] = [0.; 48];
    let mut xlen: i32 = 0;
    let mut ylen: i32 = 0;
    let mut adet: [f64; 96] = [0.; 96];
    let mut bdet: [f64; 96] = [0.; 96];
    let mut cdet: [f64; 96] = [0.; 96];
    let mut ddet: [f64; 96] = [0.; 96];
    let mut alen: i32 = 0;
    let mut blen: i32 = 0;
    let mut clen: i32 = 0;
    let mut dlen: i32 = 0;
    let mut abdet: [f64; 192] = [0.; 192];
    let mut cddet: [f64; 192] = [0.; 192];
    let mut ablen: i32 = 0;
    let mut cdlen: i32 = 0;
    let mut deter: [f64; 384] = [0.; 384];
    let mut deterlen: i32 = 0;
    let mut i: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    Two_Product(*pa.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut axby1,
                &mut axby0);
    Two_Product(*pb.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut bxay1,
                &mut bxay0);
    Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab.as_mut_ptr());
    Two_Product(*pb.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut bxcy1,
                &mut bxcy0);
    Two_Product(*pc.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut cxby1,
                &mut cxby0);
    Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc.as_mut_ptr());
    Two_Product(*pc.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut cxdy1,
                &mut cxdy0);
    Two_Product(*pd.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut dxcy1,
                &mut dxcy0);
    Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd.as_mut_ptr());
    Two_Product(*pd.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut dxay1,
                &mut dxay0);
    Two_Product(*pa.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut axdy1,
                &mut axdy0);
    Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da.as_mut_ptr());
    Two_Product(*pa.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut axcy1,
                &mut axcy0);
    Two_Product(*pc.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut cxay1,
                &mut cxay0);
    Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac.as_mut_ptr());
    Two_Product(*pb.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut bxdy1,
                &mut bxdy0);
    Two_Product(*pd.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut dxby1,
                &mut dxby0);
    Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd.as_mut_ptr());
    templen =
        fast_expansion_sum_zeroelim(4 as i32, cd.as_mut_ptr(),
                                    4 as i32, da.as_mut_ptr(),
                                    temp8.as_mut_ptr());
    cdalen =
        fast_expansion_sum_zeroelim(templen, temp8.as_mut_ptr(),
                                    4 as i32, ac.as_mut_ptr(),
                                    cda.as_mut_ptr());
    templen =
        fast_expansion_sum_zeroelim(4 as i32, da.as_mut_ptr(),
                                    4 as i32, ab.as_mut_ptr(),
                                    temp8.as_mut_ptr());
    dablen =
        fast_expansion_sum_zeroelim(templen, temp8.as_mut_ptr(),
                                    4 as i32, bd.as_mut_ptr(),
                                    dab.as_mut_ptr());
    i = 0 as i32;
    while i < 4 as i32 {
        bd[i as usize] = -bd[i as usize];
        ac[i as usize] = -ac[i as usize];
        i += 1
    }
    templen =
        fast_expansion_sum_zeroelim(4 as i32, ab.as_mut_ptr(),
                                    4 as i32, bc.as_mut_ptr(),
                                    temp8.as_mut_ptr());
    abclen =
        fast_expansion_sum_zeroelim(templen, temp8.as_mut_ptr(),
                                    4 as i32, ac.as_mut_ptr(),
                                    abc.as_mut_ptr());
    templen =
        fast_expansion_sum_zeroelim(4 as i32, bc.as_mut_ptr(),
                                    4 as i32, cd.as_mut_ptr(),
                                    temp8.as_mut_ptr());
    bcdlen =
        fast_expansion_sum_zeroelim(templen, temp8.as_mut_ptr(),
                                    4 as i32, bd.as_mut_ptr(),
                                    bcd.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(bcdlen, bcd.as_mut_ptr(),
                                 *pa.offset(0 as i32 as isize),
                                 det24x.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, det24x.as_mut_ptr(),
                                 *pa.offset(0 as i32 as isize),
                                 det48x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(bcdlen, bcd.as_mut_ptr(),
                                 *pa.offset(1 as i32 as isize),
                                 det24y.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, det24y.as_mut_ptr(),
                                 *pa.offset(1 as i32 as isize),
                                 det48y.as_mut_ptr());
    alen =
        fast_expansion_sum_zeroelim(xlen, det48x.as_mut_ptr(), ylen,
                                    det48y.as_mut_ptr(), adet.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(cdalen, cda.as_mut_ptr(),
                                 *pb.offset(0 as i32 as isize),
                                 det24x.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, det24x.as_mut_ptr(),
                                 -*pb.offset(0 as i32 as isize),
                                 det48x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(cdalen, cda.as_mut_ptr(),
                                 *pb.offset(1 as i32 as isize),
                                 det24y.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, det24y.as_mut_ptr(),
                                 -*pb.offset(1 as i32 as isize),
                                 det48y.as_mut_ptr());
    blen =
        fast_expansion_sum_zeroelim(xlen, det48x.as_mut_ptr(), ylen,
                                    det48y.as_mut_ptr(), bdet.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(dablen, dab.as_mut_ptr(),
                                 *pc.offset(0 as i32 as isize),
                                 det24x.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, det24x.as_mut_ptr(),
                                 *pc.offset(0 as i32 as isize),
                                 det48x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(dablen, dab.as_mut_ptr(),
                                 *pc.offset(1 as i32 as isize),
                                 det24y.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, det24y.as_mut_ptr(),
                                 *pc.offset(1 as i32 as isize),
                                 det48y.as_mut_ptr());
    clen =
        fast_expansion_sum_zeroelim(xlen, det48x.as_mut_ptr(), ylen,
                                    det48y.as_mut_ptr(), cdet.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(abclen, abc.as_mut_ptr(),
                                 *pd.offset(0 as i32 as isize),
                                 det24x.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, det24x.as_mut_ptr(),
                                 -*pd.offset(0 as i32 as isize),
                                 det48x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(abclen, abc.as_mut_ptr(),
                                 *pd.offset(1 as i32 as isize),
                                 det24y.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, det24y.as_mut_ptr(),
                                 -*pd.offset(1 as i32 as isize),
                                 det48y.as_mut_ptr());
    dlen =
        fast_expansion_sum_zeroelim(xlen, det48x.as_mut_ptr(), ylen,
                                    det48y.as_mut_ptr(), ddet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    cdlen =
        fast_expansion_sum_zeroelim(clen, cdet.as_mut_ptr(), dlen,
                                    ddet.as_mut_ptr(), cddet.as_mut_ptr());
    deterlen =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), cdlen,
                                    cddet.as_mut_ptr(), deter.as_mut_ptr());
    return deter[(deterlen - 1 as i32) as usize];
}

pub unsafe fn incircleslow(mut pa: *const f64,
                                      mut pb: *const f64,
                                      mut pc: *const f64,
                                      mut pd: *const f64)
 -> f64 {
    let mut adx: f64 = 0.;
    let mut bdx: f64 = 0.;
    let mut cdx: f64 = 0.;
    let mut ady: f64 = 0.;
    let mut bdy: f64 = 0.;
    let mut cdy: f64 = 0.;
    let mut adxtail: f64 = 0.;
    let mut bdxtail: f64 = 0.;
    let mut cdxtail: f64 = 0.;
    let mut adytail: f64 = 0.;
    let mut bdytail: f64 = 0.;
    let mut cdytail: f64 = 0.;
    let mut negate: f64 = 0.;
    let mut negatetail: f64 = 0.;
    let mut axby7: f64 = 0.;
    let mut bxcy7: f64 = 0.;
    let mut axcy7: f64 = 0.;
    let mut bxay7: f64 = 0.;
    let mut cxby7: f64 = 0.;
    let mut cxay7: f64 = 0.;
    let mut axby: [f64; 8] = [0.; 8];
    let mut bxcy: [f64; 8] = [0.; 8];
    let mut axcy: [f64; 8] = [0.; 8];
    let mut bxay: [f64; 8] = [0.; 8];
    let mut cxby: [f64; 8] = [0.; 8];
    let mut cxay: [f64; 8] = [0.; 8];
    let mut temp16: [f64; 16] = [0.; 16];
    let mut temp16len: i32 = 0;
    let mut detx: [f64; 32] = [0.; 32];
    let mut detxx: [f64; 64] = [0.; 64];
    let mut detxt: [f64; 32] = [0.; 32];
    let mut detxxt: [f64; 64] = [0.; 64];
    let mut detxtxt: [f64; 64] = [0.; 64];
    let mut xlen: i32 = 0;
    let mut xxlen: i32 = 0;
    let mut xtlen: i32 = 0;
    let mut xxtlen: i32 = 0;
    let mut xtxtlen: i32 = 0;
    let mut x1: [f64; 128] = [0.; 128];
    let mut x2: [f64; 192] = [0.; 192];
    let mut x1len: i32 = 0;
    let mut x2len: i32 = 0;
    let mut dety: [f64; 32] = [0.; 32];
    let mut detyy: [f64; 64] = [0.; 64];
    let mut detyt: [f64; 32] = [0.; 32];
    let mut detyyt: [f64; 64] = [0.; 64];
    let mut detytyt: [f64; 64] = [0.; 64];
    let mut ylen: i32 = 0;
    let mut yylen: i32 = 0;
    let mut ytlen: i32 = 0;
    let mut yytlen: i32 = 0;
    let mut ytytlen: i32 = 0;
    let mut y1: [f64; 128] = [0.; 128];
    let mut y2: [f64; 192] = [0.; 192];
    let mut y1len: i32 = 0;
    let mut y2len: i32 = 0;
    let mut adet: [f64; 384] = [0.; 384];
    let mut bdet: [f64; 384] = [0.; 384];
    let mut cdet: [f64; 384] = [0.; 384];
    let mut abdet: [f64; 768] = [0.; 768];
    let mut deter: [f64; 1152] = [0.; 1152];
    let mut alen: i32 = 0;
    let mut blen: i32 = 0;
    let mut clen: i32 = 0;
    let mut ablen: i32 = 0;
    let mut deterlen: i32 = 0;
    let mut i: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut a0hi: f64 = 0.;
    let mut a0lo: f64 = 0.;
    let mut a1hi: f64 = 0.;
    let mut a1lo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _k: f64 = 0.;
    let mut _l: f64 = 0.;
    let mut _m: f64 = 0.;
    let mut _n: f64 = 0.;
    let mut _0: f64 = 0.;
    let mut _1: f64 = 0.;
    let mut _2: f64 = 0.;
    Two_Diff(*pa.offset(0 as i32 as isize),
             *pd.offset(0 as i32 as isize), &mut adx, &mut adxtail);
    Two_Diff(*pa.offset(1 as i32 as isize),
             *pd.offset(1 as i32 as isize), &mut ady, &mut adytail);
    Two_Diff(*pb.offset(0 as i32 as isize),
             *pd.offset(0 as i32 as isize), &mut bdx, &mut bdxtail);
    Two_Diff(*pb.offset(1 as i32 as isize),
             *pd.offset(1 as i32 as isize), &mut bdy, &mut bdytail);
    Two_Diff(*pc.offset(0 as i32 as isize),
             *pd.offset(0 as i32 as isize), &mut cdx, &mut cdxtail);
    Two_Diff(*pc.offset(1 as i32 as isize),
             *pd.offset(1 as i32 as isize), &mut cdy, &mut cdytail);
    Two_Two_Product(adx, adxtail, bdy, bdytail, axby.as_mut_ptr());
    negate = -ady;
    negatetail = -adytail;
    Two_Two_Product(bdx, bdxtail, negate, negatetail, bxay.as_mut_ptr());
    Two_Two_Product(bdx, bdxtail, cdy, cdytail, bxcy.as_mut_ptr());
    negate = -bdy;
    negatetail = -bdytail;
    Two_Two_Product(cdx, cdxtail, negate, negatetail, cxby.as_mut_ptr());
    Two_Two_Product(cdx, cdxtail, ady, adytail, cxay.as_mut_ptr());
    negate = -cdy;
    negatetail = -cdytail;
    Two_Two_Product(adx, adxtail, negate, negatetail, axcy.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(8 as i32, bxcy.as_mut_ptr(),
                                    8 as i32, cxby.as_mut_ptr(),
                                    temp16.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), adx,
                                 detx.as_mut_ptr());
    xxlen =
        scale_expansion_zeroelim(xlen, detx.as_mut_ptr(), adx,
                                 detxx.as_mut_ptr());
    xtlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), adxtail,
                                 detxt.as_mut_ptr());
    xxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), adx,
                                 detxxt.as_mut_ptr());
    i = 0 as i32;
    while i < xxtlen { detxxt[i as usize] *= 2.0f64; i += 1 }
    xtxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), adxtail,
                                 detxtxt.as_mut_ptr());
    x1len =
        fast_expansion_sum_zeroelim(xxlen, detxx.as_mut_ptr(), xxtlen,
                                    detxxt.as_mut_ptr(), x1.as_mut_ptr());
    x2len =
        fast_expansion_sum_zeroelim(x1len, x1.as_mut_ptr(), xtxtlen,
                                    detxtxt.as_mut_ptr(), x2.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), ady,
                                 dety.as_mut_ptr());
    yylen =
        scale_expansion_zeroelim(ylen, dety.as_mut_ptr(), ady,
                                 detyy.as_mut_ptr());
    ytlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), adytail,
                                 detyt.as_mut_ptr());
    yytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), ady,
                                 detyyt.as_mut_ptr());
    i = 0 as i32;
    while i < yytlen { detyyt[i as usize] *= 2.0f64; i += 1 }
    ytytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), adytail,
                                 detytyt.as_mut_ptr());
    y1len =
        fast_expansion_sum_zeroelim(yylen, detyy.as_mut_ptr(), yytlen,
                                    detyyt.as_mut_ptr(), y1.as_mut_ptr());
    y2len =
        fast_expansion_sum_zeroelim(y1len, y1.as_mut_ptr(), ytytlen,
                                    detytyt.as_mut_ptr(), y2.as_mut_ptr());
    alen =
        fast_expansion_sum_zeroelim(x2len, x2.as_mut_ptr(), y2len,
                                    y2.as_mut_ptr(), adet.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(8 as i32, cxay.as_mut_ptr(),
                                    8 as i32, axcy.as_mut_ptr(),
                                    temp16.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), bdx,
                                 detx.as_mut_ptr());
    xxlen =
        scale_expansion_zeroelim(xlen, detx.as_mut_ptr(), bdx,
                                 detxx.as_mut_ptr());
    xtlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), bdxtail,
                                 detxt.as_mut_ptr());
    xxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), bdx,
                                 detxxt.as_mut_ptr());
    i = 0 as i32;
    while i < xxtlen { detxxt[i as usize] *= 2.0f64; i += 1 }
    xtxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), bdxtail,
                                 detxtxt.as_mut_ptr());
    x1len =
        fast_expansion_sum_zeroelim(xxlen, detxx.as_mut_ptr(), xxtlen,
                                    detxxt.as_mut_ptr(), x1.as_mut_ptr());
    x2len =
        fast_expansion_sum_zeroelim(x1len, x1.as_mut_ptr(), xtxtlen,
                                    detxtxt.as_mut_ptr(), x2.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), bdy,
                                 dety.as_mut_ptr());
    yylen =
        scale_expansion_zeroelim(ylen, dety.as_mut_ptr(), bdy,
                                 detyy.as_mut_ptr());
    ytlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), bdytail,
                                 detyt.as_mut_ptr());
    yytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), bdy,
                                 detyyt.as_mut_ptr());
    i = 0 as i32;
    while i < yytlen { detyyt[i as usize] *= 2.0f64; i += 1 }
    ytytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), bdytail,
                                 detytyt.as_mut_ptr());
    y1len =
        fast_expansion_sum_zeroelim(yylen, detyy.as_mut_ptr(), yytlen,
                                    detyyt.as_mut_ptr(), y1.as_mut_ptr());
    y2len =
        fast_expansion_sum_zeroelim(y1len, y1.as_mut_ptr(), ytytlen,
                                    detytyt.as_mut_ptr(), y2.as_mut_ptr());
    blen =
        fast_expansion_sum_zeroelim(x2len, x2.as_mut_ptr(), y2len,
                                    y2.as_mut_ptr(), bdet.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(8 as i32, axby.as_mut_ptr(),
                                    8 as i32, bxay.as_mut_ptr(),
                                    temp16.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), cdx,
                                 detx.as_mut_ptr());
    xxlen =
        scale_expansion_zeroelim(xlen, detx.as_mut_ptr(), cdx,
                                 detxx.as_mut_ptr());
    xtlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), cdxtail,
                                 detxt.as_mut_ptr());
    xxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), cdx,
                                 detxxt.as_mut_ptr());
    i = 0 as i32;
    while i < xxtlen { detxxt[i as usize] *= 2.0f64; i += 1 }
    xtxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), cdxtail,
                                 detxtxt.as_mut_ptr());
    x1len =
        fast_expansion_sum_zeroelim(xxlen, detxx.as_mut_ptr(), xxtlen,
                                    detxxt.as_mut_ptr(), x1.as_mut_ptr());
    x2len =
        fast_expansion_sum_zeroelim(x1len, x1.as_mut_ptr(), xtxtlen,
                                    detxtxt.as_mut_ptr(), x2.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), cdy,
                                 dety.as_mut_ptr());
    yylen =
        scale_expansion_zeroelim(ylen, dety.as_mut_ptr(), cdy,
                                 detyy.as_mut_ptr());
    ytlen =
        scale_expansion_zeroelim(temp16len, temp16.as_mut_ptr(), cdytail,
                                 detyt.as_mut_ptr());
    yytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), cdy,
                                 detyyt.as_mut_ptr());
    i = 0 as i32;
    while i < yytlen { detyyt[i as usize] *= 2.0f64; i += 1 }
    ytytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), cdytail,
                                 detytyt.as_mut_ptr());
    y1len =
        fast_expansion_sum_zeroelim(yylen, detyy.as_mut_ptr(), yytlen,
                                    detyyt.as_mut_ptr(), y1.as_mut_ptr());
    y2len =
        fast_expansion_sum_zeroelim(y1len, y1.as_mut_ptr(), ytytlen,
                                    detytyt.as_mut_ptr(), y2.as_mut_ptr());
    clen =
        fast_expansion_sum_zeroelim(x2len, x2.as_mut_ptr(), y2len,
                                    y2.as_mut_ptr(), cdet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    deterlen =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), clen,
                                    cdet.as_mut_ptr(), deter.as_mut_ptr());
    return deter[(deterlen - 1 as i32) as usize];
}

pub unsafe fn incircleadapt(mut pa: *const f64,
                                       mut pb: *const f64,
                                       mut pc: *const f64,
                                       mut pd: *const f64,
                                       mut permanent: f64)
 -> f64 {
    let mut adx: f64 = 0.;
    let mut bdx: f64 = 0.;
    let mut cdx: f64 = 0.;
    let mut ady: f64 = 0.;
    let mut bdy: f64 = 0.;
    let mut cdy: f64 = 0.;
    let mut det: f64 = 0.;
    let mut errbound: f64 = 0.;
    let mut bdxcdy1: f64 = 0.;
    let mut cdxbdy1: f64 = 0.;
    let mut cdxady1: f64 = 0.;
    let mut adxcdy1: f64 = 0.;
    let mut adxbdy1: f64 = 0.;
    let mut bdxady1: f64 = 0.;
    let mut bdxcdy0: f64 = 0.;
    let mut cdxbdy0: f64 = 0.;
    let mut cdxady0: f64 = 0.;
    let mut adxcdy0: f64 = 0.;
    let mut adxbdy0: f64 = 0.;
    let mut bdxady0: f64 = 0.;
    let mut bc: [f64; 4] = [0.; 4];
    let mut ca: [f64; 4] = [0.; 4];
    let mut ab: [f64; 4] = [0.; 4];
    let mut bc3: f64 = 0.;
    let mut ca3: f64 = 0.;
    let mut ab3: f64 = 0.;
    let mut axbc: [f64; 8] = [0.; 8];
    let mut axxbc: [f64; 16] = [0.; 16];
    let mut aybc: [f64; 8] = [0.; 8];
    let mut ayybc: [f64; 16] = [0.; 16];
    let mut adet: [f64; 32] = [0.; 32];
    let mut axbclen: i32 = 0;
    let mut axxbclen: i32 = 0;
    let mut aybclen: i32 = 0;
    let mut ayybclen: i32 = 0;
    let mut alen: i32 = 0;
    let mut bxca: [f64; 8] = [0.; 8];
    let mut bxxca: [f64; 16] = [0.; 16];
    let mut byca: [f64; 8] = [0.; 8];
    let mut byyca: [f64; 16] = [0.; 16];
    let mut bdet: [f64; 32] = [0.; 32];
    let mut bxcalen: i32 = 0;
    let mut bxxcalen: i32 = 0;
    let mut bycalen: i32 = 0;
    let mut byycalen: i32 = 0;
    let mut blen: i32 = 0;
    let mut cxab: [f64; 8] = [0.; 8];
    let mut cxxab: [f64; 16] = [0.; 16];
    let mut cyab: [f64; 8] = [0.; 8];
    let mut cyyab: [f64; 16] = [0.; 16];
    let mut cdet: [f64; 32] = [0.; 32];
    let mut cxablen: i32 = 0;
    let mut cxxablen: i32 = 0;
    let mut cyablen: i32 = 0;
    let mut cyyablen: i32 = 0;
    let mut clen: i32 = 0;
    let mut abdet: [f64; 64] = [0.; 64];
    let mut ablen: i32 = 0;
    let mut fin1: [f64; 1152] = [0.; 1152];
    let mut fin2: [f64; 1152] = [0.; 1152];
    let mut finnow: *mut f64 = 0 as *mut f64;
    let mut finother: *mut f64 = 0 as *mut f64;
    let mut finswap: *mut f64 = 0 as *mut f64;
    let mut finlength: i32 = 0;
    let mut adxtail: f64 = 0.;
    let mut bdxtail: f64 = 0.;
    let mut cdxtail: f64 = 0.;
    let mut adytail: f64 = 0.;
    let mut bdytail: f64 = 0.;
    let mut cdytail: f64 = 0.;
    let mut adxadx1: f64 = 0.;
    let mut adyady1: f64 = 0.;
    let mut bdxbdx1: f64 = 0.;
    let mut bdybdy1: f64 = 0.;
    let mut cdxcdx1: f64 = 0.;
    let mut cdycdy1: f64 = 0.;
    let mut adxadx0: f64 = 0.;
    let mut adyady0: f64 = 0.;
    let mut bdxbdx0: f64 = 0.;
    let mut bdybdy0: f64 = 0.;
    let mut cdxcdx0: f64 = 0.;
    let mut cdycdy0: f64 = 0.;
    let mut aa: [f64; 4] = [0.; 4];
    let mut bb: [f64; 4] = [0.; 4];
    let mut cc: [f64; 4] = [0.; 4];
    let mut aa3: f64 = 0.;
    let mut bb3: f64 = 0.;
    let mut cc3: f64 = 0.;
    let mut ti1: f64 = 0.;
    let mut tj1: f64 = 0.;
    let mut ti0: f64 = 0.;
    let mut tj0: f64 = 0.;
    let mut u: [f64; 4] = [0.; 4];
    let mut v: [f64; 4] = [0.; 4];
    let mut u3: f64 = 0.;
    let mut v3: f64 = 0.;
    let mut temp8: [f64; 8] = [0.; 8];
    let mut temp16a: [f64; 16] = [0.; 16];
    let mut temp16b: [f64; 16] = [0.; 16];
    let mut temp16c: [f64; 16] = [0.; 16];
    let mut temp32a: [f64; 32] = [0.; 32];
    let mut temp32b: [f64; 32] = [0.; 32];
    let mut temp48: [f64; 48] = [0.; 48];
    let mut temp64: [f64; 64] = [0.; 64];
    let mut temp8len: i32 = 0;
    let mut temp16alen: i32 = 0;
    let mut temp16blen: i32 = 0;
    let mut temp16clen: i32 = 0;
    let mut temp32alen: i32 = 0;
    let mut temp32blen: i32 = 0;
    let mut temp48len: i32 = 0;
    let mut temp64len: i32 = 0;
    let mut axtbb: [f64; 8] = [0.; 8];
    let mut axtcc: [f64; 8] = [0.; 8];
    let mut aytbb: [f64; 8] = [0.; 8];
    let mut aytcc: [f64; 8] = [0.; 8];
    let mut axtbblen: i32 = 0;
    let mut axtcclen: i32 = 0;
    let mut aytbblen: i32 = 0;
    let mut aytcclen: i32 = 0;
    let mut bxtaa: [f64; 8] = [0.; 8];
    let mut bxtcc: [f64; 8] = [0.; 8];
    let mut bytaa: [f64; 8] = [0.; 8];
    let mut bytcc: [f64; 8] = [0.; 8];
    let mut bxtaalen: i32 = 0;
    let mut bxtcclen: i32 = 0;
    let mut bytaalen: i32 = 0;
    let mut bytcclen: i32 = 0;
    let mut cxtaa: [f64; 8] = [0.; 8];
    let mut cxtbb: [f64; 8] = [0.; 8];
    let mut cytaa: [f64; 8] = [0.; 8];
    let mut cytbb: [f64; 8] = [0.; 8];
    let mut cxtaalen: i32 = 0;
    let mut cxtbblen: i32 = 0;
    let mut cytaalen: i32 = 0;
    let mut cytbblen: i32 = 0;
    let mut axtbc: [f64; 8] = [0.; 8];
    let mut aytbc: [f64; 8] = [0.; 8];
    let mut bxtca: [f64; 8] = [0.; 8];
    let mut bytca: [f64; 8] = [0.; 8];
    let mut cxtab: [f64; 8] = [0.; 8];
    let mut cytab: [f64; 8] = [0.; 8];
    let mut axtbclen: i32 = 0;
    let mut aytbclen: i32 = 0;
    let mut bxtcalen: i32 = 0;
    let mut bytcalen: i32 = 0;
    let mut cxtablen: i32 = 0;
    let mut cytablen: i32 = 0;
    let mut axtbct: [f64; 16] = [0.; 16];
    let mut aytbct: [f64; 16] = [0.; 16];
    let mut bxtcat: [f64; 16] = [0.; 16];
    let mut bytcat: [f64; 16] = [0.; 16];
    let mut cxtabt: [f64; 16] = [0.; 16];
    let mut cytabt: [f64; 16] = [0.; 16];
    let mut axtbctlen: i32 = 0;
    let mut aytbctlen: i32 = 0;
    let mut bxtcatlen: i32 = 0;
    let mut bytcatlen: i32 = 0;
    let mut cxtabtlen: i32 = 0;
    let mut cytabtlen: i32 = 0;
    let mut axtbctt: [f64; 8] = [0.; 8];
    let mut aytbctt: [f64; 8] = [0.; 8];
    let mut bxtcatt: [f64; 8] = [0.; 8];
    let mut bytcatt: [f64; 8] = [0.; 8];
    let mut cxtabtt: [f64; 8] = [0.; 8];
    let mut cytabtt: [f64; 8] = [0.; 8];
    let mut axtbcttlen: i32 = 0;
    let mut aytbcttlen: i32 = 0;
    let mut bxtcattlen: i32 = 0;
    let mut bytcattlen: i32 = 0;
    let mut cxtabttlen: i32 = 0;
    let mut cytabttlen: i32 = 0;
    let mut abt: [f64; 8] = [0.; 8];
    let mut bct: [f64; 8] = [0.; 8];
    let mut cat: [f64; 8] = [0.; 8];
    let mut abtlen: i32 = 0;
    let mut bctlen: i32 = 0;
    let mut catlen: i32 = 0;
    let mut abtt: [f64; 4] = [0.; 4];
    let mut bctt: [f64; 4] = [0.; 4];
    let mut catt: [f64; 4] = [0.; 4];
    let mut abttlen: i32 = 0;
    let mut bcttlen: i32 = 0;
    let mut cattlen: i32 = 0;
    let mut abtt3: f64 = 0.;
    let mut bctt3: f64 = 0.;
    let mut catt3: f64 = 0.;
    let mut negate: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    adx =
        *pa.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    bdx =
        *pb.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    cdx =
        *pc.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    ady =
        *pa.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    bdy =
        *pb.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    cdy =
        *pc.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    Two_Product(bdx, cdy, &mut bdxcdy1, &mut bdxcdy0);
    Two_Product(cdx, bdy, &mut cdxbdy1, &mut cdxbdy0);
    Two_Two_Diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0, bc.as_mut_ptr());
    axbclen =
        scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(), adx,
                                 axbc.as_mut_ptr());
    axxbclen =
        scale_expansion_zeroelim(axbclen, axbc.as_mut_ptr(), adx,
                                 axxbc.as_mut_ptr());
    aybclen =
        scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(), ady,
                                 aybc.as_mut_ptr());
    ayybclen =
        scale_expansion_zeroelim(aybclen, aybc.as_mut_ptr(), ady,
                                 ayybc.as_mut_ptr());
    alen =
        fast_expansion_sum_zeroelim(axxbclen, axxbc.as_mut_ptr(), ayybclen,
                                    ayybc.as_mut_ptr(), adet.as_mut_ptr());
    Two_Product(cdx, ady, &mut cdxady1, &mut cdxady0);
    Two_Product(adx, cdy, &mut adxcdy1, &mut adxcdy0);
    Two_Two_Diff(cdxady1, cdxady0, adxcdy1, adxcdy0, ca.as_mut_ptr());
    bxcalen =
        scale_expansion_zeroelim(4 as i32, ca.as_mut_ptr(), bdx,
                                 bxca.as_mut_ptr());
    bxxcalen =
        scale_expansion_zeroelim(bxcalen, bxca.as_mut_ptr(), bdx,
                                 bxxca.as_mut_ptr());
    bycalen =
        scale_expansion_zeroelim(4 as i32, ca.as_mut_ptr(), bdy,
                                 byca.as_mut_ptr());
    byycalen =
        scale_expansion_zeroelim(bycalen, byca.as_mut_ptr(), bdy,
                                 byyca.as_mut_ptr());
    blen =
        fast_expansion_sum_zeroelim(bxxcalen, bxxca.as_mut_ptr(), byycalen,
                                    byyca.as_mut_ptr(), bdet.as_mut_ptr());
    Two_Product(adx, bdy, &mut adxbdy1, &mut adxbdy0);
    Two_Product(bdx, ady, &mut bdxady1, &mut bdxady0);
    Two_Two_Diff(adxbdy1, adxbdy0, bdxady1, bdxady0, ab.as_mut_ptr());
    cxablen =
        scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(), cdx,
                                 cxab.as_mut_ptr());
    cxxablen =
        scale_expansion_zeroelim(cxablen, cxab.as_mut_ptr(), cdx,
                                 cxxab.as_mut_ptr());
    cyablen =
        scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(), cdy,
                                 cyab.as_mut_ptr());
    cyyablen =
        scale_expansion_zeroelim(cyablen, cyab.as_mut_ptr(), cdy,
                                 cyyab.as_mut_ptr());
    clen =
        fast_expansion_sum_zeroelim(cxxablen, cxxab.as_mut_ptr(), cyyablen,
                                    cyyab.as_mut_ptr(), cdet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    finlength =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), clen,
                                    cdet.as_mut_ptr(), fin1.as_mut_ptr());
    det = estimate(finlength, fin1.as_mut_ptr());
    errbound = PARAMS.iccerrboundB * permanent;
    if det >= errbound || -det >= errbound { return det }
    adxtail =
        Two_Diff_Tail(*pa.offset(0 as i32 as isize),
                      *pd.offset(0 as i32 as isize), adx);
    adytail =
        Two_Diff_Tail(*pa.offset(1 as i32 as isize),
                      *pd.offset(1 as i32 as isize), ady);
    bdxtail =
        Two_Diff_Tail(*pb.offset(0 as i32 as isize),
                      *pd.offset(0 as i32 as isize), bdx);
    bdytail =
        Two_Diff_Tail(*pb.offset(1 as i32 as isize),
                      *pd.offset(1 as i32 as isize), bdy);
    cdxtail =
        Two_Diff_Tail(*pc.offset(0 as i32 as isize),
                      *pd.offset(0 as i32 as isize), cdx);
    cdytail =
        Two_Diff_Tail(*pc.offset(1 as i32 as isize),
                      *pd.offset(1 as i32 as isize), cdy);
    if adxtail == 0.0f64 && bdxtail == 0.0f64 && cdxtail == 0.0f64 &&
           adytail == 0.0f64 && bdytail == 0.0f64 && cdytail == 0.0f64 {
        return det
    }
    errbound = PARAMS.iccerrboundC * permanent + PARAMS.resulterrbound * Absolute(det);
    det +=
        (adx * adx + ady * ady) *
            (bdx * cdytail + cdy * bdxtail - (bdy * cdxtail + cdx * bdytail))
            +
            2.0f64 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)
            +
            ((bdx * bdx + bdy * bdy) *
                 (cdx * adytail + ady * cdxtail -
                      (cdy * adxtail + adx * cdytail)) +
                 2.0f64 * (bdx * bdxtail + bdy * bdytail) *
                     (cdx * ady - cdy * adx)) +
            ((cdx * cdx + cdy * cdy) *
                 (adx * bdytail + bdy * adxtail -
                      (ady * bdxtail + bdx * adytail)) +
                 2.0f64 * (cdx * cdxtail + cdy * cdytail) *
                     (adx * bdy - ady * bdx));
    if det >= errbound || -det >= errbound { return det }
    finnow = fin1.as_mut_ptr();
    finother = fin2.as_mut_ptr();
    if bdxtail != 0.0f64 || bdytail != 0.0f64 || cdxtail != 0.0f64 ||
           cdytail != 0.0f64 {
        Square(adx, &mut adxadx1, &mut adxadx0);
        Square(ady, &mut adyady1, &mut adyady0);
        Two_Two_Sum(adxadx1, adxadx0, adyady1, adyady0, aa.as_mut_ptr());
    }
    if cdxtail != 0.0f64 || cdytail != 0.0f64 || adxtail != 0.0f64 ||
           adytail != 0.0f64 {
        Square(bdx, &mut bdxbdx1, &mut bdxbdx0);
        Square(bdy, &mut bdybdy1, &mut bdybdy0);
        Two_Two_Sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0, bb.as_mut_ptr());
    }
    if adxtail != 0.0f64 || adytail != 0.0f64 || bdxtail != 0.0f64 ||
           bdytail != 0.0f64 {
        Square(cdx, &mut cdxcdx1, &mut cdxcdx0);
        Square(cdy, &mut cdycdy1, &mut cdycdy0);
        Two_Two_Sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0, cc.as_mut_ptr());
    }
    if adxtail != 0.0f64 {
        axtbclen =
            scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(),
                                     adxtail, axtbc.as_mut_ptr());
        temp16alen =
            scale_expansion_zeroelim(axtbclen, axtbc.as_mut_ptr(),
                                     2.0f64 * adx, temp16a.as_mut_ptr());
        axtcclen =
            scale_expansion_zeroelim(4 as i32, cc.as_mut_ptr(),
                                     adxtail, axtcc.as_mut_ptr());
        temp16blen =
            scale_expansion_zeroelim(axtcclen, axtcc.as_mut_ptr(), bdy,
                                     temp16b.as_mut_ptr());
        axtbblen =
            scale_expansion_zeroelim(4 as i32, bb.as_mut_ptr(),
                                     adxtail, axtbb.as_mut_ptr());
        temp16clen =
            scale_expansion_zeroelim(axtbblen, axtbb.as_mut_ptr(), -cdy,
                                     temp16c.as_mut_ptr());
        temp32alen =
            fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                        temp16blen, temp16b.as_mut_ptr(),
                                        temp32a.as_mut_ptr());
        temp48len =
            fast_expansion_sum_zeroelim(temp16clen, temp16c.as_mut_ptr(),
                                        temp32alen, temp32a.as_mut_ptr(),
                                        temp48.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                        temp48.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if adytail != 0.0f64 {
        aytbclen =
            scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(),
                                     adytail, aytbc.as_mut_ptr());
        temp16alen =
            scale_expansion_zeroelim(aytbclen, aytbc.as_mut_ptr(),
                                     2.0f64 * ady, temp16a.as_mut_ptr());
        aytbblen =
            scale_expansion_zeroelim(4 as i32, bb.as_mut_ptr(),
                                     adytail, aytbb.as_mut_ptr());
        temp16blen =
            scale_expansion_zeroelim(aytbblen, aytbb.as_mut_ptr(), cdx,
                                     temp16b.as_mut_ptr());
        aytcclen =
            scale_expansion_zeroelim(4 as i32, cc.as_mut_ptr(),
                                     adytail, aytcc.as_mut_ptr());
        temp16clen =
            scale_expansion_zeroelim(aytcclen, aytcc.as_mut_ptr(), -bdx,
                                     temp16c.as_mut_ptr());
        temp32alen =
            fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                        temp16blen, temp16b.as_mut_ptr(),
                                        temp32a.as_mut_ptr());
        temp48len =
            fast_expansion_sum_zeroelim(temp16clen, temp16c.as_mut_ptr(),
                                        temp32alen, temp32a.as_mut_ptr(),
                                        temp48.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                        temp48.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if bdxtail != 0.0f64 {
        bxtcalen =
            scale_expansion_zeroelim(4 as i32, ca.as_mut_ptr(),
                                     bdxtail, bxtca.as_mut_ptr());
        temp16alen =
            scale_expansion_zeroelim(bxtcalen, bxtca.as_mut_ptr(),
                                     2.0f64 * bdx, temp16a.as_mut_ptr());
        bxtaalen =
            scale_expansion_zeroelim(4 as i32, aa.as_mut_ptr(),
                                     bdxtail, bxtaa.as_mut_ptr());
        temp16blen =
            scale_expansion_zeroelim(bxtaalen, bxtaa.as_mut_ptr(), cdy,
                                     temp16b.as_mut_ptr());
        bxtcclen =
            scale_expansion_zeroelim(4 as i32, cc.as_mut_ptr(),
                                     bdxtail, bxtcc.as_mut_ptr());
        temp16clen =
            scale_expansion_zeroelim(bxtcclen, bxtcc.as_mut_ptr(), -ady,
                                     temp16c.as_mut_ptr());
        temp32alen =
            fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                        temp16blen, temp16b.as_mut_ptr(),
                                        temp32a.as_mut_ptr());
        temp48len =
            fast_expansion_sum_zeroelim(temp16clen, temp16c.as_mut_ptr(),
                                        temp32alen, temp32a.as_mut_ptr(),
                                        temp48.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                        temp48.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if bdytail != 0.0f64 {
        bytcalen =
            scale_expansion_zeroelim(4 as i32, ca.as_mut_ptr(),
                                     bdytail, bytca.as_mut_ptr());
        temp16alen =
            scale_expansion_zeroelim(bytcalen, bytca.as_mut_ptr(),
                                     2.0f64 * bdy, temp16a.as_mut_ptr());
        bytcclen =
            scale_expansion_zeroelim(4 as i32, cc.as_mut_ptr(),
                                     bdytail, bytcc.as_mut_ptr());
        temp16blen =
            scale_expansion_zeroelim(bytcclen, bytcc.as_mut_ptr(), adx,
                                     temp16b.as_mut_ptr());
        bytaalen =
            scale_expansion_zeroelim(4 as i32, aa.as_mut_ptr(),
                                     bdytail, bytaa.as_mut_ptr());
        temp16clen =
            scale_expansion_zeroelim(bytaalen, bytaa.as_mut_ptr(), -cdx,
                                     temp16c.as_mut_ptr());
        temp32alen =
            fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                        temp16blen, temp16b.as_mut_ptr(),
                                        temp32a.as_mut_ptr());
        temp48len =
            fast_expansion_sum_zeroelim(temp16clen, temp16c.as_mut_ptr(),
                                        temp32alen, temp32a.as_mut_ptr(),
                                        temp48.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                        temp48.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if cdxtail != 0.0f64 {
        cxtablen =
            scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(),
                                     cdxtail, cxtab.as_mut_ptr());
        temp16alen =
            scale_expansion_zeroelim(cxtablen, cxtab.as_mut_ptr(),
                                     2.0f64 * cdx, temp16a.as_mut_ptr());
        cxtbblen =
            scale_expansion_zeroelim(4 as i32, bb.as_mut_ptr(),
                                     cdxtail, cxtbb.as_mut_ptr());
        temp16blen =
            scale_expansion_zeroelim(cxtbblen, cxtbb.as_mut_ptr(), ady,
                                     temp16b.as_mut_ptr());
        cxtaalen =
            scale_expansion_zeroelim(4 as i32, aa.as_mut_ptr(),
                                     cdxtail, cxtaa.as_mut_ptr());
        temp16clen =
            scale_expansion_zeroelim(cxtaalen, cxtaa.as_mut_ptr(), -bdy,
                                     temp16c.as_mut_ptr());
        temp32alen =
            fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                        temp16blen, temp16b.as_mut_ptr(),
                                        temp32a.as_mut_ptr());
        temp48len =
            fast_expansion_sum_zeroelim(temp16clen, temp16c.as_mut_ptr(),
                                        temp32alen, temp32a.as_mut_ptr(),
                                        temp48.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                        temp48.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if cdytail != 0.0f64 {
        cytablen =
            scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(),
                                     cdytail, cytab.as_mut_ptr());
        temp16alen =
            scale_expansion_zeroelim(cytablen, cytab.as_mut_ptr(),
                                     2.0f64 * cdy, temp16a.as_mut_ptr());
        cytaalen =
            scale_expansion_zeroelim(4 as i32, aa.as_mut_ptr(),
                                     cdytail, cytaa.as_mut_ptr());
        temp16blen =
            scale_expansion_zeroelim(cytaalen, cytaa.as_mut_ptr(), bdx,
                                     temp16b.as_mut_ptr());
        cytbblen =
            scale_expansion_zeroelim(4 as i32, bb.as_mut_ptr(),
                                     cdytail, cytbb.as_mut_ptr());
        temp16clen =
            scale_expansion_zeroelim(cytbblen, cytbb.as_mut_ptr(), -adx,
                                     temp16c.as_mut_ptr());
        temp32alen =
            fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                        temp16blen, temp16b.as_mut_ptr(),
                                        temp32a.as_mut_ptr());
        temp48len =
            fast_expansion_sum_zeroelim(temp16clen, temp16c.as_mut_ptr(),
                                        temp32alen, temp32a.as_mut_ptr(),
                                        temp48.as_mut_ptr());
        finlength =
            fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                        temp48.as_mut_ptr(), finother);
        finswap = finnow;
        finnow = finother;
        finother = finswap
    }
    if adxtail != 0.0f64 || adytail != 0.0f64 {
        if bdxtail != 0.0f64 || bdytail != 0.0f64 || cdxtail != 0.0f64 ||
               cdytail != 0.0f64 {
            Two_Product(bdxtail, cdy, &mut ti1, &mut ti0);
            Two_Product(bdx, cdytail, &mut tj1, &mut tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, u.as_mut_ptr());
            negate = -bdy;
            Two_Product(cdxtail, negate, &mut ti1, &mut ti0);
            negate = -bdytail;
            Two_Product(cdx, negate, &mut tj1, &mut tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, v.as_mut_ptr());
            bctlen =
                fast_expansion_sum_zeroelim(4 as i32, u.as_mut_ptr(),
                                            4 as i32, v.as_mut_ptr(),
                                            bct.as_mut_ptr());
            Two_Product(bdxtail, cdytail, &mut ti1, &mut ti0);
            Two_Product(cdxtail, bdytail, &mut tj1, &mut tj0);
            Two_Two_Diff(ti1, ti0, tj1, tj0, bctt.as_mut_ptr());
            bcttlen = 4 as i32
        } else {
            bct[0 as i32 as usize] = 0.0f64;
            bctlen = 1 as i32;
            bctt[0 as i32 as usize] = 0.0f64;
            bcttlen = 1 as i32
        }
        if adxtail != 0.0f64 {
            temp16alen =
                scale_expansion_zeroelim(axtbclen, axtbc.as_mut_ptr(),
                                         adxtail, temp16a.as_mut_ptr());
            axtbctlen =
                scale_expansion_zeroelim(bctlen, bct.as_mut_ptr(), adxtail,
                                         axtbct.as_mut_ptr());
            temp32alen =
                scale_expansion_zeroelim(axtbctlen, axtbct.as_mut_ptr(),
                                         2.0f64 * adx, temp32a.as_mut_ptr());
            temp48len =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp32alen, temp32a.as_mut_ptr(),
                                            temp48.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if bdytail != 0.0f64 {
                temp8len =
                    scale_expansion_zeroelim(4 as i32,
                                             cc.as_mut_ptr(), adxtail,
                                             temp8.as_mut_ptr());
                temp16alen =
                    scale_expansion_zeroelim(temp8len, temp8.as_mut_ptr(),
                                             bdytail, temp16a.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a.as_mut_ptr(),
                                                finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
            if cdytail != 0.0f64 {
                temp8len =
                    scale_expansion_zeroelim(4 as i32,
                                             bb.as_mut_ptr(), -adxtail,
                                             temp8.as_mut_ptr());
                temp16alen =
                    scale_expansion_zeroelim(temp8len, temp8.as_mut_ptr(),
                                             cdytail, temp16a.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a.as_mut_ptr(),
                                                finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
            temp32alen =
                scale_expansion_zeroelim(axtbctlen, axtbct.as_mut_ptr(),
                                         adxtail, temp32a.as_mut_ptr());
            axtbcttlen =
                scale_expansion_zeroelim(bcttlen, bctt.as_mut_ptr(), adxtail,
                                         axtbctt.as_mut_ptr());
            temp16alen =
                scale_expansion_zeroelim(axtbcttlen, axtbctt.as_mut_ptr(),
                                         2.0f64 * adx, temp16a.as_mut_ptr());
            temp16blen =
                scale_expansion_zeroelim(axtbcttlen, axtbctt.as_mut_ptr(),
                                         adxtail, temp16b.as_mut_ptr());
            temp32blen =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp16blen, temp16b.as_mut_ptr(),
                                            temp32b.as_mut_ptr());
            temp64len =
                fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                            temp32blen, temp32b.as_mut_ptr(),
                                            temp64.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                            temp64.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap
        }
        if adytail != 0.0f64 {
            temp16alen =
                scale_expansion_zeroelim(aytbclen, aytbc.as_mut_ptr(),
                                         adytail, temp16a.as_mut_ptr());
            aytbctlen =
                scale_expansion_zeroelim(bctlen, bct.as_mut_ptr(), adytail,
                                         aytbct.as_mut_ptr());
            temp32alen =
                scale_expansion_zeroelim(aytbctlen, aytbct.as_mut_ptr(),
                                         2.0f64 * ady, temp32a.as_mut_ptr());
            temp48len =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp32alen, temp32a.as_mut_ptr(),
                                            temp48.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            temp32alen =
                scale_expansion_zeroelim(aytbctlen, aytbct.as_mut_ptr(),
                                         adytail, temp32a.as_mut_ptr());
            aytbcttlen =
                scale_expansion_zeroelim(bcttlen, bctt.as_mut_ptr(), adytail,
                                         aytbctt.as_mut_ptr());
            temp16alen =
                scale_expansion_zeroelim(aytbcttlen, aytbctt.as_mut_ptr(),
                                         2.0f64 * ady, temp16a.as_mut_ptr());
            temp16blen =
                scale_expansion_zeroelim(aytbcttlen, aytbctt.as_mut_ptr(),
                                         adytail, temp16b.as_mut_ptr());
            temp32blen =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp16blen, temp16b.as_mut_ptr(),
                                            temp32b.as_mut_ptr());
            temp64len =
                fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                            temp32blen, temp32b.as_mut_ptr(),
                                            temp64.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                            temp64.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap
        }
    }
    if bdxtail != 0.0f64 || bdytail != 0.0f64 {
        if cdxtail != 0.0f64 || cdytail != 0.0f64 || adxtail != 0.0f64 ||
               adytail != 0.0f64 {
            Two_Product(cdxtail, ady, &mut ti1, &mut ti0);
            Two_Product(cdx, adytail, &mut tj1, &mut tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, u.as_mut_ptr());
            negate = -cdy;
            Two_Product(adxtail, negate, &mut ti1, &mut ti0);
            negate = -cdytail;
            Two_Product(adx, negate, &mut tj1, &mut tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, v.as_mut_ptr());
            catlen =
                fast_expansion_sum_zeroelim(4 as i32, u.as_mut_ptr(),
                                            4 as i32, v.as_mut_ptr(),
                                            cat.as_mut_ptr());
            Two_Product(cdxtail, adytail, &mut ti1, &mut ti0);
            Two_Product(adxtail, cdytail, &mut tj1, &mut tj0);
            Two_Two_Diff(ti1, ti0, tj1, tj0, catt.as_mut_ptr());
            cattlen = 4 as i32
        } else {
            cat[0 as i32 as usize] = 0.0f64;
            catlen = 1 as i32;
            catt[0 as i32 as usize] = 0.0f64;
            cattlen = 1 as i32
        }
        if bdxtail != 0.0f64 {
            temp16alen =
                scale_expansion_zeroelim(bxtcalen, bxtca.as_mut_ptr(),
                                         bdxtail, temp16a.as_mut_ptr());
            bxtcatlen =
                scale_expansion_zeroelim(catlen, cat.as_mut_ptr(), bdxtail,
                                         bxtcat.as_mut_ptr());
            temp32alen =
                scale_expansion_zeroelim(bxtcatlen, bxtcat.as_mut_ptr(),
                                         2.0f64 * bdx, temp32a.as_mut_ptr());
            temp48len =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp32alen, temp32a.as_mut_ptr(),
                                            temp48.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if cdytail != 0.0f64 {
                temp8len =
                    scale_expansion_zeroelim(4 as i32,
                                             aa.as_mut_ptr(), bdxtail,
                                             temp8.as_mut_ptr());
                temp16alen =
                    scale_expansion_zeroelim(temp8len, temp8.as_mut_ptr(),
                                             cdytail, temp16a.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a.as_mut_ptr(),
                                                finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
            if adytail != 0.0f64 {
                temp8len =
                    scale_expansion_zeroelim(4 as i32,
                                             cc.as_mut_ptr(), -bdxtail,
                                             temp8.as_mut_ptr());
                temp16alen =
                    scale_expansion_zeroelim(temp8len, temp8.as_mut_ptr(),
                                             adytail, temp16a.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a.as_mut_ptr(),
                                                finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
            temp32alen =
                scale_expansion_zeroelim(bxtcatlen, bxtcat.as_mut_ptr(),
                                         bdxtail, temp32a.as_mut_ptr());
            bxtcattlen =
                scale_expansion_zeroelim(cattlen, catt.as_mut_ptr(), bdxtail,
                                         bxtcatt.as_mut_ptr());
            temp16alen =
                scale_expansion_zeroelim(bxtcattlen, bxtcatt.as_mut_ptr(),
                                         2.0f64 * bdx, temp16a.as_mut_ptr());
            temp16blen =
                scale_expansion_zeroelim(bxtcattlen, bxtcatt.as_mut_ptr(),
                                         bdxtail, temp16b.as_mut_ptr());
            temp32blen =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp16blen, temp16b.as_mut_ptr(),
                                            temp32b.as_mut_ptr());
            temp64len =
                fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                            temp32blen, temp32b.as_mut_ptr(),
                                            temp64.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                            temp64.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap
        }
        if bdytail != 0.0f64 {
            temp16alen =
                scale_expansion_zeroelim(bytcalen, bytca.as_mut_ptr(),
                                         bdytail, temp16a.as_mut_ptr());
            bytcatlen =
                scale_expansion_zeroelim(catlen, cat.as_mut_ptr(), bdytail,
                                         bytcat.as_mut_ptr());
            temp32alen =
                scale_expansion_zeroelim(bytcatlen, bytcat.as_mut_ptr(),
                                         2.0f64 * bdy, temp32a.as_mut_ptr());
            temp48len =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp32alen, temp32a.as_mut_ptr(),
                                            temp48.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            temp32alen =
                scale_expansion_zeroelim(bytcatlen, bytcat.as_mut_ptr(),
                                         bdytail, temp32a.as_mut_ptr());
            bytcattlen =
                scale_expansion_zeroelim(cattlen, catt.as_mut_ptr(), bdytail,
                                         bytcatt.as_mut_ptr());
            temp16alen =
                scale_expansion_zeroelim(bytcattlen, bytcatt.as_mut_ptr(),
                                         2.0f64 * bdy, temp16a.as_mut_ptr());
            temp16blen =
                scale_expansion_zeroelim(bytcattlen, bytcatt.as_mut_ptr(),
                                         bdytail, temp16b.as_mut_ptr());
            temp32blen =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp16blen, temp16b.as_mut_ptr(),
                                            temp32b.as_mut_ptr());
            temp64len =
                fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                            temp32blen, temp32b.as_mut_ptr(),
                                            temp64.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                            temp64.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap
        }
    }
    if cdxtail != 0.0f64 || cdytail != 0.0f64 {
        if adxtail != 0.0f64 || adytail != 0.0f64 || bdxtail != 0.0f64 ||
               bdytail != 0.0f64 {
            Two_Product(adxtail, bdy, &mut ti1, &mut ti0);
            Two_Product(adx, bdytail, &mut tj1, &mut tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, u.as_mut_ptr());
            negate = -ady;
            Two_Product(bdxtail, negate, &mut ti1, &mut ti0);
            negate = -adytail;
            Two_Product(bdx, negate, &mut tj1, &mut tj0);
            Two_Two_Sum(ti1, ti0, tj1, tj0, v.as_mut_ptr());
            abtlen =
                fast_expansion_sum_zeroelim(4 as i32, u.as_mut_ptr(),
                                            4 as i32, v.as_mut_ptr(),
                                            abt.as_mut_ptr());
            Two_Product(adxtail, bdytail, &mut ti1, &mut ti0);
            Two_Product(bdxtail, adytail, &mut tj1, &mut tj0);
            Two_Two_Diff(ti1, ti0, tj1, tj0, abtt.as_mut_ptr());
            abttlen = 4 as i32
        } else {
            abt[0 as i32 as usize] = 0.0f64;
            abtlen = 1 as i32;
            abtt[0 as i32 as usize] = 0.0f64;
            abttlen = 1 as i32
        }
        if cdxtail != 0.0f64 {
            temp16alen =
                scale_expansion_zeroelim(cxtablen, cxtab.as_mut_ptr(),
                                         cdxtail, temp16a.as_mut_ptr());
            cxtabtlen =
                scale_expansion_zeroelim(abtlen, abt.as_mut_ptr(), cdxtail,
                                         cxtabt.as_mut_ptr());
            temp32alen =
                scale_expansion_zeroelim(cxtabtlen, cxtabt.as_mut_ptr(),
                                         2.0f64 * cdx, temp32a.as_mut_ptr());
            temp48len =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp32alen, temp32a.as_mut_ptr(),
                                            temp48.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            if adytail != 0.0f64 {
                temp8len =
                    scale_expansion_zeroelim(4 as i32,
                                             bb.as_mut_ptr(), cdxtail,
                                             temp8.as_mut_ptr());
                temp16alen =
                    scale_expansion_zeroelim(temp8len, temp8.as_mut_ptr(),
                                             adytail, temp16a.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a.as_mut_ptr(),
                                                finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
            if bdytail != 0.0f64 {
                temp8len =
                    scale_expansion_zeroelim(4 as i32,
                                             aa.as_mut_ptr(), -cdxtail,
                                             temp8.as_mut_ptr());
                temp16alen =
                    scale_expansion_zeroelim(temp8len, temp8.as_mut_ptr(),
                                             bdytail, temp16a.as_mut_ptr());
                finlength =
                    fast_expansion_sum_zeroelim(finlength, finnow, temp16alen,
                                                temp16a.as_mut_ptr(),
                                                finother);
                finswap = finnow;
                finnow = finother;
                finother = finswap
            }
            temp32alen =
                scale_expansion_zeroelim(cxtabtlen, cxtabt.as_mut_ptr(),
                                         cdxtail, temp32a.as_mut_ptr());
            cxtabttlen =
                scale_expansion_zeroelim(abttlen, abtt.as_mut_ptr(), cdxtail,
                                         cxtabtt.as_mut_ptr());
            temp16alen =
                scale_expansion_zeroelim(cxtabttlen, cxtabtt.as_mut_ptr(),
                                         2.0f64 * cdx, temp16a.as_mut_ptr());
            temp16blen =
                scale_expansion_zeroelim(cxtabttlen, cxtabtt.as_mut_ptr(),
                                         cdxtail, temp16b.as_mut_ptr());
            temp32blen =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp16blen, temp16b.as_mut_ptr(),
                                            temp32b.as_mut_ptr());
            temp64len =
                fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                            temp32blen, temp32b.as_mut_ptr(),
                                            temp64.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                            temp64.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap
        }
        if cdytail != 0.0f64 {
            temp16alen =
                scale_expansion_zeroelim(cytablen, cytab.as_mut_ptr(),
                                         cdytail, temp16a.as_mut_ptr());
            cytabtlen =
                scale_expansion_zeroelim(abtlen, abt.as_mut_ptr(), cdytail,
                                         cytabt.as_mut_ptr());
            temp32alen =
                scale_expansion_zeroelim(cytabtlen, cytabt.as_mut_ptr(),
                                         2.0f64 * cdy, temp32a.as_mut_ptr());
            temp48len =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp32alen, temp32a.as_mut_ptr(),
                                            temp48.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp48len,
                                            temp48.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap;
            temp32alen =
                scale_expansion_zeroelim(cytabtlen, cytabt.as_mut_ptr(),
                                         cdytail, temp32a.as_mut_ptr());
            cytabttlen =
                scale_expansion_zeroelim(abttlen, abtt.as_mut_ptr(), cdytail,
                                         cytabtt.as_mut_ptr());
            temp16alen =
                scale_expansion_zeroelim(cytabttlen, cytabtt.as_mut_ptr(),
                                         2.0f64 * cdy, temp16a.as_mut_ptr());
            temp16blen =
                scale_expansion_zeroelim(cytabttlen, cytabtt.as_mut_ptr(),
                                         cdytail, temp16b.as_mut_ptr());
            temp32blen =
                fast_expansion_sum_zeroelim(temp16alen, temp16a.as_mut_ptr(),
                                            temp16blen, temp16b.as_mut_ptr(),
                                            temp32b.as_mut_ptr());
            temp64len =
                fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                            temp32blen, temp32b.as_mut_ptr(),
                                            temp64.as_mut_ptr());
            finlength =
                fast_expansion_sum_zeroelim(finlength, finnow, temp64len,
                                            temp64.as_mut_ptr(), finother);
            finswap = finnow;
            finnow = finother;
            finother = finswap
        }
    }
    return *finnow.offset((finlength - 1 as i32) as isize);
}

pub unsafe fn incircle(mut pa: *const f64,
                                  mut pb: *const f64,
                                  mut pc: *const f64,
                                  mut pd: *const f64)
 -> f64 {
    let mut adx: f64 = 0.;
    let mut bdx: f64 = 0.;
    let mut cdx: f64 = 0.;
    let mut ady: f64 = 0.;
    let mut bdy: f64 = 0.;
    let mut cdy: f64 = 0.;
    let mut bdxcdy: f64 = 0.;
    let mut cdxbdy: f64 = 0.;
    let mut cdxady: f64 = 0.;
    let mut adxcdy: f64 = 0.;
    let mut adxbdy: f64 = 0.;
    let mut bdxady: f64 = 0.;
    let mut alift: f64 = 0.;
    let mut blift: f64 = 0.;
    let mut clift: f64 = 0.;
    let mut det: f64 = 0.;
    let mut permanent: f64 = 0.;
    let mut errbound: f64 = 0.;
    adx =
        *pa.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    bdx =
        *pb.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    cdx =
        *pc.offset(0 as i32 as isize) -
            *pd.offset(0 as i32 as isize);
    ady =
        *pa.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    bdy =
        *pb.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    cdy =
        *pc.offset(1 as i32 as isize) -
            *pd.offset(1 as i32 as isize);
    bdxcdy = bdx * cdy;
    cdxbdy = cdx * bdy;
    alift = adx * adx + ady * ady;
    cdxady = cdx * ady;
    adxcdy = adx * cdy;
    blift = bdx * bdx + bdy * bdy;
    adxbdy = adx * bdy;
    bdxady = bdx * ady;
    clift = cdx * cdx + cdy * cdy;
    det =
        alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) +
            clift * (adxbdy - bdxady);
    permanent =
        (Absolute(bdxcdy) + Absolute(cdxbdy)) * alift +
            (Absolute(cdxady) + Absolute(adxcdy)) * blift +
            (Absolute(adxbdy) + Absolute(bdxady)) * clift;
    errbound = PARAMS.iccerrboundA * permanent;
    if det > errbound || -det > errbound { return det }
    return incircleadapt(pa, pb, pc, pd, permanent);
}
/* ****************************************************************************/
/*                                                                           */
/*  inspherefast()   Approximate 3D insphere test.  Nonrobust.               */
/*  insphereexact()   Exact 3D insphere test.  Robust.                       */
/*  insphereslow()   Another exact 3D insphere test.  Robust.                */
/*  insphere()   Adaptive exact 3D insphere test.  Robust.                   */
/*                                                                           */
/*               Return a positive value if the point pe lies inside the     */
/*               sphere passing through pa, pb, pc, and pd; a negative value */
/*               if it lies outside; and zero if the five points are         */
/*               cospherical.  The points pa, pb, pc, and pd must be ordered */
/*               so that they have a positive orientation (as defined by     */
/*               orient3d()), or the sign of the result will be reversed.    */
/*                                                                           */
/*  Only the first and last routine should be used; the middle two are for   */
/*  timings.                                                                 */
/*                                                                           */
/*  The last three use exact arithmetic to ensure a correct answer.  The     */
/*  result returned is the determinant of a matrix.  In insphere() only,     */
/*  this determinant is computed adaptively, in the sense that exact         */
/*  arithmetic is used only to the degree it is needed to ensure that the    */
/*  returned value has the correct sign.  Hence, insphere() is usually quite */
/*  fast, but will run more slowly when the input points are cospherical or  */
/*  nearly so.                                                               */
/*                                                                           */
/* ****************************************************************************/

pub unsafe fn inspherefast(mut pa: *const f64,
                                      mut pb: *const f64,
                                      mut pc: *const f64,
                                      mut pd: *const f64,
                                      mut pe: *const f64)
 -> f64 {
    let mut aex: f64 = 0.;
    let mut bex: f64 = 0.;
    let mut cex: f64 = 0.;
    let mut dex: f64 = 0.;
    let mut aey: f64 = 0.;
    let mut bey: f64 = 0.;
    let mut cey: f64 = 0.;
    let mut dey: f64 = 0.;
    let mut aez: f64 = 0.;
    let mut bez: f64 = 0.;
    let mut cez: f64 = 0.;
    let mut dez: f64 = 0.;
    let mut alift: f64 = 0.;
    let mut blift: f64 = 0.;
    let mut clift: f64 = 0.;
    let mut dlift: f64 = 0.;
    let mut ab: f64 = 0.;
    let mut bc: f64 = 0.;
    let mut cd: f64 = 0.;
    let mut da: f64 = 0.;
    let mut ac: f64 = 0.;
    let mut bd: f64 = 0.;
    let mut abc: f64 = 0.;
    let mut bcd: f64 = 0.;
    let mut cda: f64 = 0.;
    let mut dab: f64 = 0.;
    aex =
        *pa.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    bex =
        *pb.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    cex =
        *pc.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    dex =
        *pd.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    aey =
        *pa.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    bey =
        *pb.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    cey =
        *pc.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    dey =
        *pd.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    aez =
        *pa.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    bez =
        *pb.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    cez =
        *pc.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    dez =
        *pd.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    ab = aex * bey - bex * aey;
    bc = bex * cey - cex * bey;
    cd = cex * dey - dex * cey;
    da = dex * aey - aex * dey;
    ac = aex * cey - cex * aey;
    bd = bex * dey - dex * bey;
    abc = aez * bc - bez * ac + cez * ab;
    bcd = bez * cd - cez * bd + dez * bc;
    cda = cez * da + dez * ac + aez * cd;
    dab = dez * ab + aez * bd + bez * da;
    alift = aex * aex + aey * aey + aez * aez;
    blift = bex * bex + bey * bey + bez * bez;
    clift = cex * cex + cey * cey + cez * cez;
    dlift = dex * dex + dey * dey + dez * dez;
    return dlift * abc - clift * dab + (blift * cda - alift * bcd);
}

pub unsafe fn insphereexact(mut pa: *const f64,
                                       mut pb: *const f64,
                                       mut pc: *const f64,
                                       mut pd: *const f64,
                                       mut pe: *const f64)
 -> f64 {
    let mut axby1: f64 = 0.;
    let mut bxcy1: f64 = 0.;
    let mut cxdy1: f64 = 0.;
    let mut dxey1: f64 = 0.;
    let mut exay1: f64 = 0.;
    let mut bxay1: f64 = 0.;
    let mut cxby1: f64 = 0.;
    let mut dxcy1: f64 = 0.;
    let mut exdy1: f64 = 0.;
    let mut axey1: f64 = 0.;
    let mut axcy1: f64 = 0.;
    let mut bxdy1: f64 = 0.;
    let mut cxey1: f64 = 0.;
    let mut dxay1: f64 = 0.;
    let mut exby1: f64 = 0.;
    let mut cxay1: f64 = 0.;
    let mut dxby1: f64 = 0.;
    let mut excy1: f64 = 0.;
    let mut axdy1: f64 = 0.;
    let mut bxey1: f64 = 0.;
    let mut axby0: f64 = 0.;
    let mut bxcy0: f64 = 0.;
    let mut cxdy0: f64 = 0.;
    let mut dxey0: f64 = 0.;
    let mut exay0: f64 = 0.;
    let mut bxay0: f64 = 0.;
    let mut cxby0: f64 = 0.;
    let mut dxcy0: f64 = 0.;
    let mut exdy0: f64 = 0.;
    let mut axey0: f64 = 0.;
    let mut axcy0: f64 = 0.;
    let mut bxdy0: f64 = 0.;
    let mut cxey0: f64 = 0.;
    let mut dxay0: f64 = 0.;
    let mut exby0: f64 = 0.;
    let mut cxay0: f64 = 0.;
    let mut dxby0: f64 = 0.;
    let mut excy0: f64 = 0.;
    let mut axdy0: f64 = 0.;
    let mut bxey0: f64 = 0.;
    let mut ab: [f64; 4] = [0.; 4];
    let mut bc: [f64; 4] = [0.; 4];
    let mut cd: [f64; 4] = [0.; 4];
    let mut de: [f64; 4] = [0.; 4];
    let mut ea: [f64; 4] = [0.; 4];
    let mut ac: [f64; 4] = [0.; 4];
    let mut bd: [f64; 4] = [0.; 4];
    let mut ce: [f64; 4] = [0.; 4];
    let mut da: [f64; 4] = [0.; 4];
    let mut eb: [f64; 4] = [0.; 4];
    let mut temp8a: [f64; 8] = [0.; 8];
    let mut temp8b: [f64; 8] = [0.; 8];
    let mut temp16: [f64; 16] = [0.; 16];
    let mut temp8alen: i32 = 0;
    let mut temp8blen: i32 = 0;
    let mut temp16len: i32 = 0;
    let mut abc: [f64; 24] = [0.; 24];
    let mut bcd: [f64; 24] = [0.; 24];
    let mut cde: [f64; 24] = [0.; 24];
    let mut dea: [f64; 24] = [0.; 24];
    let mut eab: [f64; 24] = [0.; 24];
    let mut abd: [f64; 24] = [0.; 24];
    let mut bce: [f64; 24] = [0.; 24];
    let mut cda: [f64; 24] = [0.; 24];
    let mut deb: [f64; 24] = [0.; 24];
    let mut eac: [f64; 24] = [0.; 24];
    let mut abclen: i32 = 0;
    let mut bcdlen: i32 = 0;
    let mut cdelen: i32 = 0;
    let mut dealen: i32 = 0;
    let mut eablen: i32 = 0;
    let mut abdlen: i32 = 0;
    let mut bcelen: i32 = 0;
    let mut cdalen: i32 = 0;
    let mut deblen: i32 = 0;
    let mut eaclen: i32 = 0;
    let mut temp48a: [f64; 48] = [0.; 48];
    let mut temp48b: [f64; 48] = [0.; 48];
    let mut temp48alen: i32 = 0;
    let mut temp48blen: i32 = 0;
    let mut abcd: [f64; 96] = [0.; 96];
    let mut bcde: [f64; 96] = [0.; 96];
    let mut cdea: [f64; 96] = [0.; 96];
    let mut deab: [f64; 96] = [0.; 96];
    let mut eabc: [f64; 96] = [0.; 96];
    let mut abcdlen: i32 = 0;
    let mut bcdelen: i32 = 0;
    let mut cdealen: i32 = 0;
    let mut deablen: i32 = 0;
    let mut eabclen: i32 = 0;
    let mut temp192: [f64; 192] = [0.; 192];
    let mut det384x: [f64; 384] = [0.; 384];
    let mut det384y: [f64; 384] = [0.; 384];
    let mut det384z: [f64; 384] = [0.; 384];
    let mut xlen: i32 = 0;
    let mut ylen: i32 = 0;
    let mut zlen: i32 = 0;
    let mut detxy: [f64; 768] = [0.; 768];
    let mut xylen: i32 = 0;
    let mut adet: [f64; 1152] = [0.; 1152];
    let mut bdet: [f64; 1152] = [0.; 1152];
    let mut cdet: [f64; 1152] = [0.; 1152];
    let mut ddet: [f64; 1152] = [0.; 1152];
    let mut edet: [f64; 1152] = [0.; 1152];
    let mut alen: i32 = 0;
    let mut blen: i32 = 0;
    let mut clen: i32 = 0;
    let mut dlen: i32 = 0;
    let mut elen: i32 = 0;
    let mut abdet: [f64; 2304] = [0.; 2304];
    let mut cddet: [f64; 2304] = [0.; 2304];
    let mut cdedet: [f64; 3456] = [0.; 3456];
    let mut ablen: i32 = 0;
    let mut cdlen: i32 = 0;
    let mut deter: [f64; 5760] = [0.; 5760];
    let mut deterlen: i32 = 0;
    let mut i: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    Two_Product(*pa.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut axby1,
                &mut axby0);
    Two_Product(*pb.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut bxay1,
                &mut bxay0);
    Two_Two_Diff(axby1, axby0, bxay1, bxay0, ab.as_mut_ptr());
    Two_Product(*pb.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut bxcy1,
                &mut bxcy0);
    Two_Product(*pc.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut cxby1,
                &mut cxby0);
    Two_Two_Diff(bxcy1, bxcy0, cxby1, cxby0, bc.as_mut_ptr());
    Two_Product(*pc.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut cxdy1,
                &mut cxdy0);
    Two_Product(*pd.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut dxcy1,
                &mut dxcy0);
    Two_Two_Diff(cxdy1, cxdy0, dxcy1, dxcy0, cd.as_mut_ptr());
    Two_Product(*pd.offset(0 as i32 as isize),
                *pe.offset(1 as i32 as isize), &mut dxey1,
                &mut dxey0);
    Two_Product(*pe.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut exdy1,
                &mut exdy0);
    Two_Two_Diff(dxey1, dxey0, exdy1, exdy0, de.as_mut_ptr());
    Two_Product(*pe.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut exay1,
                &mut exay0);
    Two_Product(*pa.offset(0 as i32 as isize),
                *pe.offset(1 as i32 as isize), &mut axey1,
                &mut axey0);
    Two_Two_Diff(exay1, exay0, axey1, axey0, ea.as_mut_ptr());
    Two_Product(*pa.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut axcy1,
                &mut axcy0);
    Two_Product(*pc.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut cxay1,
                &mut cxay0);
    Two_Two_Diff(axcy1, axcy0, cxay1, cxay0, ac.as_mut_ptr());
    Two_Product(*pb.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut bxdy1,
                &mut bxdy0);
    Two_Product(*pd.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut dxby1,
                &mut dxby0);
    Two_Two_Diff(bxdy1, bxdy0, dxby1, dxby0, bd.as_mut_ptr());
    Two_Product(*pc.offset(0 as i32 as isize),
                *pe.offset(1 as i32 as isize), &mut cxey1,
                &mut cxey0);
    Two_Product(*pe.offset(0 as i32 as isize),
                *pc.offset(1 as i32 as isize), &mut excy1,
                &mut excy0);
    Two_Two_Diff(cxey1, cxey0, excy1, excy0, ce.as_mut_ptr());
    Two_Product(*pd.offset(0 as i32 as isize),
                *pa.offset(1 as i32 as isize), &mut dxay1,
                &mut dxay0);
    Two_Product(*pa.offset(0 as i32 as isize),
                *pd.offset(1 as i32 as isize), &mut axdy1,
                &mut axdy0);
    Two_Two_Diff(dxay1, dxay0, axdy1, axdy0, da.as_mut_ptr());
    Two_Product(*pe.offset(0 as i32 as isize),
                *pb.offset(1 as i32 as isize), &mut exby1,
                &mut exby0);
    Two_Product(*pb.offset(0 as i32 as isize),
                *pe.offset(1 as i32 as isize), &mut bxey1,
                &mut bxey0);
    Two_Two_Diff(exby1, exby0, bxey1, bxey0, eb.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(),
                                 *pa.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, ac.as_mut_ptr(),
                                 -*pb.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(),
                                 *pc.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    abclen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), abc.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, cd.as_mut_ptr(),
                                 *pb.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, bd.as_mut_ptr(),
                                 -*pc.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(),
                                 *pd.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    bcdlen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), bcd.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, de.as_mut_ptr(),
                                 *pc.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, ce.as_mut_ptr(),
                                 -*pd.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, cd.as_mut_ptr(),
                                 *pe.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    cdelen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), cde.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ea.as_mut_ptr(),
                                 *pd.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, da.as_mut_ptr(),
                                 -*pe.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, de.as_mut_ptr(),
                                 *pa.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    dealen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), dea.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(),
                                 *pe.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, eb.as_mut_ptr(),
                                 -*pa.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ea.as_mut_ptr(),
                                 *pb.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    eablen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), eab.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, bd.as_mut_ptr(),
                                 *pa.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, da.as_mut_ptr(),
                                 *pb.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(),
                                 *pd.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    abdlen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), abd.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ce.as_mut_ptr(),
                                 *pb.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, eb.as_mut_ptr(),
                                 *pc.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(),
                                 *pe.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    bcelen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), bce.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, da.as_mut_ptr(),
                                 *pc.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, ac.as_mut_ptr(),
                                 *pd.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, cd.as_mut_ptr(),
                                 *pa.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    cdalen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), cda.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, eb.as_mut_ptr(),
                                 *pd.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, bd.as_mut_ptr(),
                                 *pe.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, de.as_mut_ptr(),
                                 *pb.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    deblen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), deb.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ac.as_mut_ptr(),
                                 *pe.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, ce.as_mut_ptr(),
                                 *pa.offset(2 as i32 as isize),
                                 temp8b.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ea.as_mut_ptr(),
                                 *pc.offset(2 as i32 as isize),
                                 temp8a.as_mut_ptr());
    eaclen =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), eac.as_mut_ptr());
    temp48alen =
        fast_expansion_sum_zeroelim(cdelen, cde.as_mut_ptr(), bcelen,
                                    bce.as_mut_ptr(), temp48a.as_mut_ptr());
    temp48blen =
        fast_expansion_sum_zeroelim(deblen, deb.as_mut_ptr(), bcdlen,
                                    bcd.as_mut_ptr(), temp48b.as_mut_ptr());
    i = 0 as i32;
    while i < temp48blen {
        temp48b[i as usize] = -temp48b[i as usize];
        i += 1
    }
    bcdelen =
        fast_expansion_sum_zeroelim(temp48alen, temp48a.as_mut_ptr(),
                                    temp48blen, temp48b.as_mut_ptr(),
                                    bcde.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(bcdelen, bcde.as_mut_ptr(),
                                 *pa.offset(0 as i32 as isize),
                                 temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, temp192.as_mut_ptr(),
                                 *pa.offset(0 as i32 as isize),
                                 det384x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(bcdelen, bcde.as_mut_ptr(),
                                 *pa.offset(1 as i32 as isize),
                                 temp192.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, temp192.as_mut_ptr(),
                                 *pa.offset(1 as i32 as isize),
                                 det384y.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(bcdelen, bcde.as_mut_ptr(),
                                 *pa.offset(2 as i32 as isize),
                                 temp192.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(zlen, temp192.as_mut_ptr(),
                                 *pa.offset(2 as i32 as isize),
                                 det384z.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, det384x.as_mut_ptr(), ylen,
                                    det384y.as_mut_ptr(), detxy.as_mut_ptr());
    alen =
        fast_expansion_sum_zeroelim(xylen, detxy.as_mut_ptr(), zlen,
                                    det384z.as_mut_ptr(), adet.as_mut_ptr());
    temp48alen =
        fast_expansion_sum_zeroelim(dealen, dea.as_mut_ptr(), cdalen,
                                    cda.as_mut_ptr(), temp48a.as_mut_ptr());
    temp48blen =
        fast_expansion_sum_zeroelim(eaclen, eac.as_mut_ptr(), cdelen,
                                    cde.as_mut_ptr(), temp48b.as_mut_ptr());
    i = 0 as i32;
    while i < temp48blen {
        temp48b[i as usize] = -temp48b[i as usize];
        i += 1
    }
    cdealen =
        fast_expansion_sum_zeroelim(temp48alen, temp48a.as_mut_ptr(),
                                    temp48blen, temp48b.as_mut_ptr(),
                                    cdea.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(cdealen, cdea.as_mut_ptr(),
                                 *pb.offset(0 as i32 as isize),
                                 temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, temp192.as_mut_ptr(),
                                 *pb.offset(0 as i32 as isize),
                                 det384x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(cdealen, cdea.as_mut_ptr(),
                                 *pb.offset(1 as i32 as isize),
                                 temp192.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, temp192.as_mut_ptr(),
                                 *pb.offset(1 as i32 as isize),
                                 det384y.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(cdealen, cdea.as_mut_ptr(),
                                 *pb.offset(2 as i32 as isize),
                                 temp192.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(zlen, temp192.as_mut_ptr(),
                                 *pb.offset(2 as i32 as isize),
                                 det384z.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, det384x.as_mut_ptr(), ylen,
                                    det384y.as_mut_ptr(), detxy.as_mut_ptr());
    blen =
        fast_expansion_sum_zeroelim(xylen, detxy.as_mut_ptr(), zlen,
                                    det384z.as_mut_ptr(), bdet.as_mut_ptr());
    temp48alen =
        fast_expansion_sum_zeroelim(eablen, eab.as_mut_ptr(), deblen,
                                    deb.as_mut_ptr(), temp48a.as_mut_ptr());
    temp48blen =
        fast_expansion_sum_zeroelim(abdlen, abd.as_mut_ptr(), dealen,
                                    dea.as_mut_ptr(), temp48b.as_mut_ptr());
    i = 0 as i32;
    while i < temp48blen {
        temp48b[i as usize] = -temp48b[i as usize];
        i += 1
    }
    deablen =
        fast_expansion_sum_zeroelim(temp48alen, temp48a.as_mut_ptr(),
                                    temp48blen, temp48b.as_mut_ptr(),
                                    deab.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(deablen, deab.as_mut_ptr(),
                                 *pc.offset(0 as i32 as isize),
                                 temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, temp192.as_mut_ptr(),
                                 *pc.offset(0 as i32 as isize),
                                 det384x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(deablen, deab.as_mut_ptr(),
                                 *pc.offset(1 as i32 as isize),
                                 temp192.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, temp192.as_mut_ptr(),
                                 *pc.offset(1 as i32 as isize),
                                 det384y.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(deablen, deab.as_mut_ptr(),
                                 *pc.offset(2 as i32 as isize),
                                 temp192.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(zlen, temp192.as_mut_ptr(),
                                 *pc.offset(2 as i32 as isize),
                                 det384z.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, det384x.as_mut_ptr(), ylen,
                                    det384y.as_mut_ptr(), detxy.as_mut_ptr());
    clen =
        fast_expansion_sum_zeroelim(xylen, detxy.as_mut_ptr(), zlen,
                                    det384z.as_mut_ptr(), cdet.as_mut_ptr());
    temp48alen =
        fast_expansion_sum_zeroelim(abclen, abc.as_mut_ptr(), eaclen,
                                    eac.as_mut_ptr(), temp48a.as_mut_ptr());
    temp48blen =
        fast_expansion_sum_zeroelim(bcelen, bce.as_mut_ptr(), eablen,
                                    eab.as_mut_ptr(), temp48b.as_mut_ptr());
    i = 0 as i32;
    while i < temp48blen {
        temp48b[i as usize] = -temp48b[i as usize];
        i += 1
    }
    eabclen =
        fast_expansion_sum_zeroelim(temp48alen, temp48a.as_mut_ptr(),
                                    temp48blen, temp48b.as_mut_ptr(),
                                    eabc.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(eabclen, eabc.as_mut_ptr(),
                                 *pd.offset(0 as i32 as isize),
                                 temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, temp192.as_mut_ptr(),
                                 *pd.offset(0 as i32 as isize),
                                 det384x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(eabclen, eabc.as_mut_ptr(),
                                 *pd.offset(1 as i32 as isize),
                                 temp192.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, temp192.as_mut_ptr(),
                                 *pd.offset(1 as i32 as isize),
                                 det384y.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(eabclen, eabc.as_mut_ptr(),
                                 *pd.offset(2 as i32 as isize),
                                 temp192.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(zlen, temp192.as_mut_ptr(),
                                 *pd.offset(2 as i32 as isize),
                                 det384z.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, det384x.as_mut_ptr(), ylen,
                                    det384y.as_mut_ptr(), detxy.as_mut_ptr());
    dlen =
        fast_expansion_sum_zeroelim(xylen, detxy.as_mut_ptr(), zlen,
                                    det384z.as_mut_ptr(), ddet.as_mut_ptr());
    temp48alen =
        fast_expansion_sum_zeroelim(bcdlen, bcd.as_mut_ptr(), abdlen,
                                    abd.as_mut_ptr(), temp48a.as_mut_ptr());
    temp48blen =
        fast_expansion_sum_zeroelim(cdalen, cda.as_mut_ptr(), abclen,
                                    abc.as_mut_ptr(), temp48b.as_mut_ptr());
    i = 0 as i32;
    while i < temp48blen {
        temp48b[i as usize] = -temp48b[i as usize];
        i += 1
    }
    abcdlen =
        fast_expansion_sum_zeroelim(temp48alen, temp48a.as_mut_ptr(),
                                    temp48blen, temp48b.as_mut_ptr(),
                                    abcd.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(abcdlen, abcd.as_mut_ptr(),
                                 *pe.offset(0 as i32 as isize),
                                 temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(xlen, temp192.as_mut_ptr(),
                                 *pe.offset(0 as i32 as isize),
                                 det384x.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(abcdlen, abcd.as_mut_ptr(),
                                 *pe.offset(1 as i32 as isize),
                                 temp192.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(ylen, temp192.as_mut_ptr(),
                                 *pe.offset(1 as i32 as isize),
                                 det384y.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(abcdlen, abcd.as_mut_ptr(),
                                 *pe.offset(2 as i32 as isize),
                                 temp192.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(zlen, temp192.as_mut_ptr(),
                                 *pe.offset(2 as i32 as isize),
                                 det384z.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, det384x.as_mut_ptr(), ylen,
                                    det384y.as_mut_ptr(), detxy.as_mut_ptr());
    elen =
        fast_expansion_sum_zeroelim(xylen, detxy.as_mut_ptr(), zlen,
                                    det384z.as_mut_ptr(), edet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    cdlen =
        fast_expansion_sum_zeroelim(clen, cdet.as_mut_ptr(), dlen,
                                    ddet.as_mut_ptr(), cddet.as_mut_ptr());
    cdelen =
        fast_expansion_sum_zeroelim(cdlen, cddet.as_mut_ptr(), elen,
                                    edet.as_mut_ptr(), cdedet.as_mut_ptr());
    deterlen =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), cdelen,
                                    cdedet.as_mut_ptr(), deter.as_mut_ptr());
    return deter[(deterlen - 1 as i32) as usize];
}

pub unsafe fn insphereslow(mut pa: *const f64,
                                      mut pb: *const f64,
                                      mut pc: *const f64,
                                      mut pd: *const f64,
                                      mut pe: *const f64)
 -> f64 {
    let mut aex: f64 = 0.;
    let mut bex: f64 = 0.;
    let mut cex: f64 = 0.;
    let mut dex: f64 = 0.;
    let mut aey: f64 = 0.;
    let mut bey: f64 = 0.;
    let mut cey: f64 = 0.;
    let mut dey: f64 = 0.;
    let mut aez: f64 = 0.;
    let mut bez: f64 = 0.;
    let mut cez: f64 = 0.;
    let mut dez: f64 = 0.;
    let mut aextail: f64 = 0.;
    let mut bextail: f64 = 0.;
    let mut cextail: f64 = 0.;
    let mut dextail: f64 = 0.;
    let mut aeytail: f64 = 0.;
    let mut beytail: f64 = 0.;
    let mut ceytail: f64 = 0.;
    let mut deytail: f64 = 0.;
    let mut aeztail: f64 = 0.;
    let mut beztail: f64 = 0.;
    let mut ceztail: f64 = 0.;
    let mut deztail: f64 = 0.;
    let mut negate: f64 = 0.;
    let mut negatetail: f64 = 0.;
    let mut axby7: f64 = 0.;
    let mut bxcy7: f64 = 0.;
    let mut cxdy7: f64 = 0.;
    let mut dxay7: f64 = 0.;
    let mut axcy7: f64 = 0.;
    let mut bxdy7: f64 = 0.;
    let mut bxay7: f64 = 0.;
    let mut cxby7: f64 = 0.;
    let mut dxcy7: f64 = 0.;
    let mut axdy7: f64 = 0.;
    let mut cxay7: f64 = 0.;
    let mut dxby7: f64 = 0.;
    let mut axby: [f64; 8] = [0.; 8];
    let mut bxcy: [f64; 8] = [0.; 8];
    let mut cxdy: [f64; 8] = [0.; 8];
    let mut dxay: [f64; 8] = [0.; 8];
    let mut axcy: [f64; 8] = [0.; 8];
    let mut bxdy: [f64; 8] = [0.; 8];
    let mut bxay: [f64; 8] = [0.; 8];
    let mut cxby: [f64; 8] = [0.; 8];
    let mut dxcy: [f64; 8] = [0.; 8];
    let mut axdy: [f64; 8] = [0.; 8];
    let mut cxay: [f64; 8] = [0.; 8];
    let mut dxby: [f64; 8] = [0.; 8];
    let mut ab: [f64; 16] = [0.; 16];
    let mut bc: [f64; 16] = [0.; 16];
    let mut cd: [f64; 16] = [0.; 16];
    let mut da: [f64; 16] = [0.; 16];
    let mut ac: [f64; 16] = [0.; 16];
    let mut bd: [f64; 16] = [0.; 16];
    let mut ablen: i32 = 0;
    let mut bclen: i32 = 0;
    let mut cdlen: i32 = 0;
    let mut dalen: i32 = 0;
    let mut aclen: i32 = 0;
    let mut bdlen: i32 = 0;
    let mut temp32a: [f64; 32] = [0.; 32];
    let mut temp32b: [f64; 32] = [0.; 32];
    let mut temp64a: [f64; 64] = [0.; 64];
    let mut temp64b: [f64; 64] = [0.; 64];
    let mut temp64c: [f64; 64] = [0.; 64];
    let mut temp32alen: i32 = 0;
    let mut temp32blen: i32 = 0;
    let mut temp64alen: i32 = 0;
    let mut temp64blen: i32 = 0;
    let mut temp64clen: i32 = 0;
    let mut temp128: [f64; 128] = [0.; 128];
    let mut temp192: [f64; 192] = [0.; 192];
    let mut temp128len: i32 = 0;
    let mut temp192len: i32 = 0;
    let mut detx: [f64; 384] = [0.; 384];
    let mut detxx: [f64; 768] = [0.; 768];
    let mut detxt: [f64; 384] = [0.; 384];
    let mut detxxt: [f64; 768] = [0.; 768];
    let mut detxtxt: [f64; 768] = [0.; 768];
    let mut xlen: i32 = 0;
    let mut xxlen: i32 = 0;
    let mut xtlen: i32 = 0;
    let mut xxtlen: i32 = 0;
    let mut xtxtlen: i32 = 0;
    let mut x1: [f64; 1536] = [0.; 1536];
    let mut x2: [f64; 2304] = [0.; 2304];
    let mut x1len: i32 = 0;
    let mut x2len: i32 = 0;
    let mut dety: [f64; 384] = [0.; 384];
    let mut detyy: [f64; 768] = [0.; 768];
    let mut detyt: [f64; 384] = [0.; 384];
    let mut detyyt: [f64; 768] = [0.; 768];
    let mut detytyt: [f64; 768] = [0.; 768];
    let mut ylen: i32 = 0;
    let mut yylen: i32 = 0;
    let mut ytlen: i32 = 0;
    let mut yytlen: i32 = 0;
    let mut ytytlen: i32 = 0;
    let mut y1: [f64; 1536] = [0.; 1536];
    let mut y2: [f64; 2304] = [0.; 2304];
    let mut y1len: i32 = 0;
    let mut y2len: i32 = 0;
    let mut detz: [f64; 384] = [0.; 384];
    let mut detzz: [f64; 768] = [0.; 768];
    let mut detzt: [f64; 384] = [0.; 384];
    let mut detzzt: [f64; 768] = [0.; 768];
    let mut detztzt: [f64; 768] = [0.; 768];
    let mut zlen: i32 = 0;
    let mut zzlen: i32 = 0;
    let mut ztlen: i32 = 0;
    let mut zztlen: i32 = 0;
    let mut ztztlen: i32 = 0;
    let mut z1: [f64; 1536] = [0.; 1536];
    let mut z2: [f64; 2304] = [0.; 2304];
    let mut z1len: i32 = 0;
    let mut z2len: i32 = 0;
    let mut detxy: [f64; 4608] = [0.; 4608];
    let mut xylen: i32 = 0;
    let mut adet: [f64; 6912] = [0.; 6912];
    let mut bdet: [f64; 6912] = [0.; 6912];
    let mut cdet: [f64; 6912] = [0.; 6912];
    let mut ddet: [f64; 6912] = [0.; 6912];
    let mut alen: i32 = 0;
    let mut blen: i32 = 0;
    let mut clen: i32 = 0;
    let mut dlen: i32 = 0;
    let mut abdet: [f64; 13824] = [0.; 13824];
    let mut cddet: [f64; 13824] = [0.; 13824];
    let mut deter: [f64; 27648] = [0.; 27648];
    let mut deterlen: i32 = 0;
    let mut i: i32 = 0;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut a0hi: f64 = 0.;
    let mut a0lo: f64 = 0.;
    let mut a1hi: f64 = 0.;
    let mut a1lo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _k: f64 = 0.;
    let mut _l: f64 = 0.;
    let mut _m: f64 = 0.;
    let mut _n: f64 = 0.;
    let mut _0: f64 = 0.;
    let mut _1: f64 = 0.;
    let mut _2: f64 = 0.;
    Two_Diff(*pa.offset(0 as i32 as isize),
             *pe.offset(0 as i32 as isize), &mut aex, &mut aextail);
    Two_Diff(*pa.offset(1 as i32 as isize),
             *pe.offset(1 as i32 as isize), &mut aey, &mut aeytail);
    Two_Diff(*pa.offset(2 as i32 as isize),
             *pe.offset(2 as i32 as isize), &mut aez, &mut aeztail);
    Two_Diff(*pb.offset(0 as i32 as isize),
             *pe.offset(0 as i32 as isize), &mut bex, &mut bextail);
    Two_Diff(*pb.offset(1 as i32 as isize),
             *pe.offset(1 as i32 as isize), &mut bey, &mut beytail);
    Two_Diff(*pb.offset(2 as i32 as isize),
             *pe.offset(2 as i32 as isize), &mut bez, &mut beztail);
    Two_Diff(*pc.offset(0 as i32 as isize),
             *pe.offset(0 as i32 as isize), &mut cex, &mut cextail);
    Two_Diff(*pc.offset(1 as i32 as isize),
             *pe.offset(1 as i32 as isize), &mut cey, &mut ceytail);
    Two_Diff(*pc.offset(2 as i32 as isize),
             *pe.offset(2 as i32 as isize), &mut cez, &mut ceztail);
    Two_Diff(*pd.offset(0 as i32 as isize),
             *pe.offset(0 as i32 as isize), &mut dex, &mut dextail);
    Two_Diff(*pd.offset(1 as i32 as isize),
             *pe.offset(1 as i32 as isize), &mut dey, &mut deytail);
    Two_Diff(*pd.offset(2 as i32 as isize),
             *pe.offset(2 as i32 as isize), &mut dez, &mut deztail);
    Two_Two_Product(aex, aextail, bey, beytail, axby.as_mut_ptr());
    negate = -aey;
    negatetail = -aeytail;
    Two_Two_Product(bex, bextail, negate, negatetail, bxay.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(8 as i32, axby.as_mut_ptr(),
                                    8 as i32, bxay.as_mut_ptr(),
                                    ab.as_mut_ptr());
    Two_Two_Product(bex, bextail, cey, ceytail, bxcy.as_mut_ptr());
    negate = -bey;
    negatetail = -beytail;
    Two_Two_Product(cex, cextail, negate, negatetail, cxby.as_mut_ptr());
    bclen =
        fast_expansion_sum_zeroelim(8 as i32, bxcy.as_mut_ptr(),
                                    8 as i32, cxby.as_mut_ptr(),
                                    bc.as_mut_ptr());
    Two_Two_Product(cex, cextail, dey, deytail, cxdy.as_mut_ptr());
    negate = -cey;
    negatetail = -ceytail;
    Two_Two_Product(dex, dextail, negate, negatetail, dxcy.as_mut_ptr());
    cdlen =
        fast_expansion_sum_zeroelim(8 as i32, cxdy.as_mut_ptr(),
                                    8 as i32, dxcy.as_mut_ptr(),
                                    cd.as_mut_ptr());
    Two_Two_Product(dex, dextail, aey, aeytail, dxay.as_mut_ptr());
    negate = -dey;
    negatetail = -deytail;
    Two_Two_Product(aex, aextail, negate, negatetail, axdy.as_mut_ptr());
    dalen =
        fast_expansion_sum_zeroelim(8 as i32, dxay.as_mut_ptr(),
                                    8 as i32, axdy.as_mut_ptr(),
                                    da.as_mut_ptr());
    Two_Two_Product(aex, aextail, cey, ceytail, axcy.as_mut_ptr());
    negate = -aey;
    negatetail = -aeytail;
    Two_Two_Product(cex, cextail, negate, negatetail, cxay.as_mut_ptr());
    aclen =
        fast_expansion_sum_zeroelim(8 as i32, axcy.as_mut_ptr(),
                                    8 as i32, cxay.as_mut_ptr(),
                                    ac.as_mut_ptr());
    Two_Two_Product(bex, bextail, dey, deytail, bxdy.as_mut_ptr());
    negate = -bey;
    negatetail = -beytail;
    Two_Two_Product(dex, dextail, negate, negatetail, dxby.as_mut_ptr());
    bdlen =
        fast_expansion_sum_zeroelim(8 as i32, bxdy.as_mut_ptr(),
                                    8 as i32, dxby.as_mut_ptr(),
                                    bd.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(cdlen, cd.as_mut_ptr(), -bez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(cdlen, cd.as_mut_ptr(), -beztail,
                                 temp32b.as_mut_ptr());
    temp64alen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64a.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(bdlen, bd.as_mut_ptr(), cez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(bdlen, bd.as_mut_ptr(), ceztail,
                                 temp32b.as_mut_ptr());
    temp64blen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64b.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(bclen, bc.as_mut_ptr(), -dez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(bclen, bc.as_mut_ptr(), -deztail,
                                 temp32b.as_mut_ptr());
    temp64clen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64c.as_mut_ptr());
    temp128len =
        fast_expansion_sum_zeroelim(temp64alen, temp64a.as_mut_ptr(),
                                    temp64blen, temp64b.as_mut_ptr(),
                                    temp128.as_mut_ptr());
    temp192len =
        fast_expansion_sum_zeroelim(temp64clen, temp64c.as_mut_ptr(),
                                    temp128len, temp128.as_mut_ptr(),
                                    temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), aex,
                                 detx.as_mut_ptr());
    xxlen =
        scale_expansion_zeroelim(xlen, detx.as_mut_ptr(), aex,
                                 detxx.as_mut_ptr());
    xtlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), aextail,
                                 detxt.as_mut_ptr());
    xxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), aex,
                                 detxxt.as_mut_ptr());
    i = 0 as i32;
    while i < xxtlen { detxxt[i as usize] *= 2.0f64; i += 1 }
    xtxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), aextail,
                                 detxtxt.as_mut_ptr());
    x1len =
        fast_expansion_sum_zeroelim(xxlen, detxx.as_mut_ptr(), xxtlen,
                                    detxxt.as_mut_ptr(), x1.as_mut_ptr());
    x2len =
        fast_expansion_sum_zeroelim(x1len, x1.as_mut_ptr(), xtxtlen,
                                    detxtxt.as_mut_ptr(), x2.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), aey,
                                 dety.as_mut_ptr());
    yylen =
        scale_expansion_zeroelim(ylen, dety.as_mut_ptr(), aey,
                                 detyy.as_mut_ptr());
    ytlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), aeytail,
                                 detyt.as_mut_ptr());
    yytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), aey,
                                 detyyt.as_mut_ptr());
    i = 0 as i32;
    while i < yytlen { detyyt[i as usize] *= 2.0f64; i += 1 }
    ytytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), aeytail,
                                 detytyt.as_mut_ptr());
    y1len =
        fast_expansion_sum_zeroelim(yylen, detyy.as_mut_ptr(), yytlen,
                                    detyyt.as_mut_ptr(), y1.as_mut_ptr());
    y2len =
        fast_expansion_sum_zeroelim(y1len, y1.as_mut_ptr(), ytytlen,
                                    detytyt.as_mut_ptr(), y2.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), aez,
                                 detz.as_mut_ptr());
    zzlen =
        scale_expansion_zeroelim(zlen, detz.as_mut_ptr(), aez,
                                 detzz.as_mut_ptr());
    ztlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), aeztail,
                                 detzt.as_mut_ptr());
    zztlen =
        scale_expansion_zeroelim(ztlen, detzt.as_mut_ptr(), aez,
                                 detzzt.as_mut_ptr());
    i = 0 as i32;
    while i < zztlen { detzzt[i as usize] *= 2.0f64; i += 1 }
    ztztlen =
        scale_expansion_zeroelim(ztlen, detzt.as_mut_ptr(), aeztail,
                                 detztzt.as_mut_ptr());
    z1len =
        fast_expansion_sum_zeroelim(zzlen, detzz.as_mut_ptr(), zztlen,
                                    detzzt.as_mut_ptr(), z1.as_mut_ptr());
    z2len =
        fast_expansion_sum_zeroelim(z1len, z1.as_mut_ptr(), ztztlen,
                                    detztzt.as_mut_ptr(), z2.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(x2len, x2.as_mut_ptr(), y2len,
                                    y2.as_mut_ptr(), detxy.as_mut_ptr());
    alen =
        fast_expansion_sum_zeroelim(z2len, z2.as_mut_ptr(), xylen,
                                    detxy.as_mut_ptr(), adet.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(dalen, da.as_mut_ptr(), cez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(dalen, da.as_mut_ptr(), ceztail,
                                 temp32b.as_mut_ptr());
    temp64alen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64a.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(aclen, ac.as_mut_ptr(), dez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(aclen, ac.as_mut_ptr(), deztail,
                                 temp32b.as_mut_ptr());
    temp64blen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64b.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(cdlen, cd.as_mut_ptr(), aez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(cdlen, cd.as_mut_ptr(), aeztail,
                                 temp32b.as_mut_ptr());
    temp64clen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64c.as_mut_ptr());
    temp128len =
        fast_expansion_sum_zeroelim(temp64alen, temp64a.as_mut_ptr(),
                                    temp64blen, temp64b.as_mut_ptr(),
                                    temp128.as_mut_ptr());
    temp192len =
        fast_expansion_sum_zeroelim(temp64clen, temp64c.as_mut_ptr(),
                                    temp128len, temp128.as_mut_ptr(),
                                    temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), bex,
                                 detx.as_mut_ptr());
    xxlen =
        scale_expansion_zeroelim(xlen, detx.as_mut_ptr(), bex,
                                 detxx.as_mut_ptr());
    xtlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), bextail,
                                 detxt.as_mut_ptr());
    xxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), bex,
                                 detxxt.as_mut_ptr());
    i = 0 as i32;
    while i < xxtlen { detxxt[i as usize] *= 2.0f64; i += 1 }
    xtxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), bextail,
                                 detxtxt.as_mut_ptr());
    x1len =
        fast_expansion_sum_zeroelim(xxlen, detxx.as_mut_ptr(), xxtlen,
                                    detxxt.as_mut_ptr(), x1.as_mut_ptr());
    x2len =
        fast_expansion_sum_zeroelim(x1len, x1.as_mut_ptr(), xtxtlen,
                                    detxtxt.as_mut_ptr(), x2.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), bey,
                                 dety.as_mut_ptr());
    yylen =
        scale_expansion_zeroelim(ylen, dety.as_mut_ptr(), bey,
                                 detyy.as_mut_ptr());
    ytlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), beytail,
                                 detyt.as_mut_ptr());
    yytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), bey,
                                 detyyt.as_mut_ptr());
    i = 0 as i32;
    while i < yytlen { detyyt[i as usize] *= 2.0f64; i += 1 }
    ytytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), beytail,
                                 detytyt.as_mut_ptr());
    y1len =
        fast_expansion_sum_zeroelim(yylen, detyy.as_mut_ptr(), yytlen,
                                    detyyt.as_mut_ptr(), y1.as_mut_ptr());
    y2len =
        fast_expansion_sum_zeroelim(y1len, y1.as_mut_ptr(), ytytlen,
                                    detytyt.as_mut_ptr(), y2.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), bez,
                                 detz.as_mut_ptr());
    zzlen =
        scale_expansion_zeroelim(zlen, detz.as_mut_ptr(), bez,
                                 detzz.as_mut_ptr());
    ztlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), beztail,
                                 detzt.as_mut_ptr());
    zztlen =
        scale_expansion_zeroelim(ztlen, detzt.as_mut_ptr(), bez,
                                 detzzt.as_mut_ptr());
    i = 0 as i32;
    while i < zztlen { detzzt[i as usize] *= 2.0f64; i += 1 }
    ztztlen =
        scale_expansion_zeroelim(ztlen, detzt.as_mut_ptr(), beztail,
                                 detztzt.as_mut_ptr());
    z1len =
        fast_expansion_sum_zeroelim(zzlen, detzz.as_mut_ptr(), zztlen,
                                    detzzt.as_mut_ptr(), z1.as_mut_ptr());
    z2len =
        fast_expansion_sum_zeroelim(z1len, z1.as_mut_ptr(), ztztlen,
                                    detztzt.as_mut_ptr(), z2.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(x2len, x2.as_mut_ptr(), y2len,
                                    y2.as_mut_ptr(), detxy.as_mut_ptr());
    blen =
        fast_expansion_sum_zeroelim(z2len, z2.as_mut_ptr(), xylen,
                                    detxy.as_mut_ptr(), bdet.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(ablen, ab.as_mut_ptr(), -dez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(ablen, ab.as_mut_ptr(), -deztail,
                                 temp32b.as_mut_ptr());
    temp64alen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64a.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(bdlen, bd.as_mut_ptr(), -aez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(bdlen, bd.as_mut_ptr(), -aeztail,
                                 temp32b.as_mut_ptr());
    temp64blen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64b.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(dalen, da.as_mut_ptr(), -bez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(dalen, da.as_mut_ptr(), -beztail,
                                 temp32b.as_mut_ptr());
    temp64clen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64c.as_mut_ptr());
    temp128len =
        fast_expansion_sum_zeroelim(temp64alen, temp64a.as_mut_ptr(),
                                    temp64blen, temp64b.as_mut_ptr(),
                                    temp128.as_mut_ptr());
    temp192len =
        fast_expansion_sum_zeroelim(temp64clen, temp64c.as_mut_ptr(),
                                    temp128len, temp128.as_mut_ptr(),
                                    temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), cex,
                                 detx.as_mut_ptr());
    xxlen =
        scale_expansion_zeroelim(xlen, detx.as_mut_ptr(), cex,
                                 detxx.as_mut_ptr());
    xtlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), cextail,
                                 detxt.as_mut_ptr());
    xxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), cex,
                                 detxxt.as_mut_ptr());
    i = 0 as i32;
    while i < xxtlen { detxxt[i as usize] *= 2.0f64; i += 1 }
    xtxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), cextail,
                                 detxtxt.as_mut_ptr());
    x1len =
        fast_expansion_sum_zeroelim(xxlen, detxx.as_mut_ptr(), xxtlen,
                                    detxxt.as_mut_ptr(), x1.as_mut_ptr());
    x2len =
        fast_expansion_sum_zeroelim(x1len, x1.as_mut_ptr(), xtxtlen,
                                    detxtxt.as_mut_ptr(), x2.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), cey,
                                 dety.as_mut_ptr());
    yylen =
        scale_expansion_zeroelim(ylen, dety.as_mut_ptr(), cey,
                                 detyy.as_mut_ptr());
    ytlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), ceytail,
                                 detyt.as_mut_ptr());
    yytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), cey,
                                 detyyt.as_mut_ptr());
    i = 0 as i32;
    while i < yytlen { detyyt[i as usize] *= 2.0f64; i += 1 }
    ytytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), ceytail,
                                 detytyt.as_mut_ptr());
    y1len =
        fast_expansion_sum_zeroelim(yylen, detyy.as_mut_ptr(), yytlen,
                                    detyyt.as_mut_ptr(), y1.as_mut_ptr());
    y2len =
        fast_expansion_sum_zeroelim(y1len, y1.as_mut_ptr(), ytytlen,
                                    detytyt.as_mut_ptr(), y2.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), cez,
                                 detz.as_mut_ptr());
    zzlen =
        scale_expansion_zeroelim(zlen, detz.as_mut_ptr(), cez,
                                 detzz.as_mut_ptr());
    ztlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), ceztail,
                                 detzt.as_mut_ptr());
    zztlen =
        scale_expansion_zeroelim(ztlen, detzt.as_mut_ptr(), cez,
                                 detzzt.as_mut_ptr());
    i = 0 as i32;
    while i < zztlen { detzzt[i as usize] *= 2.0f64; i += 1 }
    ztztlen =
        scale_expansion_zeroelim(ztlen, detzt.as_mut_ptr(), ceztail,
                                 detztzt.as_mut_ptr());
    z1len =
        fast_expansion_sum_zeroelim(zzlen, detzz.as_mut_ptr(), zztlen,
                                    detzzt.as_mut_ptr(), z1.as_mut_ptr());
    z2len =
        fast_expansion_sum_zeroelim(z1len, z1.as_mut_ptr(), ztztlen,
                                    detztzt.as_mut_ptr(), z2.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(x2len, x2.as_mut_ptr(), y2len,
                                    y2.as_mut_ptr(), detxy.as_mut_ptr());
    clen =
        fast_expansion_sum_zeroelim(z2len, z2.as_mut_ptr(), xylen,
                                    detxy.as_mut_ptr(), cdet.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(bclen, bc.as_mut_ptr(), aez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(bclen, bc.as_mut_ptr(), aeztail,
                                 temp32b.as_mut_ptr());
    temp64alen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64a.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(aclen, ac.as_mut_ptr(), -bez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(aclen, ac.as_mut_ptr(), -beztail,
                                 temp32b.as_mut_ptr());
    temp64blen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64b.as_mut_ptr());
    temp32alen =
        scale_expansion_zeroelim(ablen, ab.as_mut_ptr(), cez,
                                 temp32a.as_mut_ptr());
    temp32blen =
        scale_expansion_zeroelim(ablen, ab.as_mut_ptr(), ceztail,
                                 temp32b.as_mut_ptr());
    temp64clen =
        fast_expansion_sum_zeroelim(temp32alen, temp32a.as_mut_ptr(),
                                    temp32blen, temp32b.as_mut_ptr(),
                                    temp64c.as_mut_ptr());
    temp128len =
        fast_expansion_sum_zeroelim(temp64alen, temp64a.as_mut_ptr(),
                                    temp64blen, temp64b.as_mut_ptr(),
                                    temp128.as_mut_ptr());
    temp192len =
        fast_expansion_sum_zeroelim(temp64clen, temp64c.as_mut_ptr(),
                                    temp128len, temp128.as_mut_ptr(),
                                    temp192.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), dex,
                                 detx.as_mut_ptr());
    xxlen =
        scale_expansion_zeroelim(xlen, detx.as_mut_ptr(), dex,
                                 detxx.as_mut_ptr());
    xtlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), dextail,
                                 detxt.as_mut_ptr());
    xxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), dex,
                                 detxxt.as_mut_ptr());
    i = 0 as i32;
    while i < xxtlen { detxxt[i as usize] *= 2.0f64; i += 1 }
    xtxtlen =
        scale_expansion_zeroelim(xtlen, detxt.as_mut_ptr(), dextail,
                                 detxtxt.as_mut_ptr());
    x1len =
        fast_expansion_sum_zeroelim(xxlen, detxx.as_mut_ptr(), xxtlen,
                                    detxxt.as_mut_ptr(), x1.as_mut_ptr());
    x2len =
        fast_expansion_sum_zeroelim(x1len, x1.as_mut_ptr(), xtxtlen,
                                    detxtxt.as_mut_ptr(), x2.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), dey,
                                 dety.as_mut_ptr());
    yylen =
        scale_expansion_zeroelim(ylen, dety.as_mut_ptr(), dey,
                                 detyy.as_mut_ptr());
    ytlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), deytail,
                                 detyt.as_mut_ptr());
    yytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), dey,
                                 detyyt.as_mut_ptr());
    i = 0 as i32;
    while i < yytlen { detyyt[i as usize] *= 2.0f64; i += 1 }
    ytytlen =
        scale_expansion_zeroelim(ytlen, detyt.as_mut_ptr(), deytail,
                                 detytyt.as_mut_ptr());
    y1len =
        fast_expansion_sum_zeroelim(yylen, detyy.as_mut_ptr(), yytlen,
                                    detyyt.as_mut_ptr(), y1.as_mut_ptr());
    y2len =
        fast_expansion_sum_zeroelim(y1len, y1.as_mut_ptr(), ytytlen,
                                    detytyt.as_mut_ptr(), y2.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), dez,
                                 detz.as_mut_ptr());
    zzlen =
        scale_expansion_zeroelim(zlen, detz.as_mut_ptr(), dez,
                                 detzz.as_mut_ptr());
    ztlen =
        scale_expansion_zeroelim(temp192len, temp192.as_mut_ptr(), deztail,
                                 detzt.as_mut_ptr());
    zztlen =
        scale_expansion_zeroelim(ztlen, detzt.as_mut_ptr(), dez,
                                 detzzt.as_mut_ptr());
    i = 0 as i32;
    while i < zztlen { detzzt[i as usize] *= 2.0f64; i += 1 }
    ztztlen =
        scale_expansion_zeroelim(ztlen, detzt.as_mut_ptr(), deztail,
                                 detztzt.as_mut_ptr());
    z1len =
        fast_expansion_sum_zeroelim(zzlen, detzz.as_mut_ptr(), zztlen,
                                    detzzt.as_mut_ptr(), z1.as_mut_ptr());
    z2len =
        fast_expansion_sum_zeroelim(z1len, z1.as_mut_ptr(), ztztlen,
                                    detztzt.as_mut_ptr(), z2.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(x2len, x2.as_mut_ptr(), y2len,
                                    y2.as_mut_ptr(), detxy.as_mut_ptr());
    dlen =
        fast_expansion_sum_zeroelim(z2len, z2.as_mut_ptr(), xylen,
                                    detxy.as_mut_ptr(), ddet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    cdlen =
        fast_expansion_sum_zeroelim(clen, cdet.as_mut_ptr(), dlen,
                                    ddet.as_mut_ptr(), cddet.as_mut_ptr());
    deterlen =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), cdlen,
                                    cddet.as_mut_ptr(), deter.as_mut_ptr());
    return deter[(deterlen - 1 as i32) as usize];
}

pub unsafe fn insphereadapt(mut pa: *const f64,
                                       mut pb: *const f64,
                                       mut pc: *const f64,
                                       mut pd: *const f64,
                                       mut pe: *const f64,
                                       mut permanent: f64)
 -> f64 {
    let mut aex: f64 = 0.;
    let mut bex: f64 = 0.;
    let mut cex: f64 = 0.;
    let mut dex: f64 = 0.;
    let mut aey: f64 = 0.;
    let mut bey: f64 = 0.;
    let mut cey: f64 = 0.;
    let mut dey: f64 = 0.;
    let mut aez: f64 = 0.;
    let mut bez: f64 = 0.;
    let mut cez: f64 = 0.;
    let mut dez: f64 = 0.;
    let mut det: f64 = 0.;
    let mut errbound: f64 = 0.;
    let mut aexbey1: f64 = 0.;
    let mut bexaey1: f64 = 0.;
    let mut bexcey1: f64 = 0.;
    let mut cexbey1: f64 = 0.;
    let mut cexdey1: f64 = 0.;
    let mut dexcey1: f64 = 0.;
    let mut dexaey1: f64 = 0.;
    let mut aexdey1: f64 = 0.;
    let mut aexcey1: f64 = 0.;
    let mut cexaey1: f64 = 0.;
    let mut bexdey1: f64 = 0.;
    let mut dexbey1: f64 = 0.;
    let mut aexbey0: f64 = 0.;
    let mut bexaey0: f64 = 0.;
    let mut bexcey0: f64 = 0.;
    let mut cexbey0: f64 = 0.;
    let mut cexdey0: f64 = 0.;
    let mut dexcey0: f64 = 0.;
    let mut dexaey0: f64 = 0.;
    let mut aexdey0: f64 = 0.;
    let mut aexcey0: f64 = 0.;
    let mut cexaey0: f64 = 0.;
    let mut bexdey0: f64 = 0.;
    let mut dexbey0: f64 = 0.;
    let mut ab: [f64; 4] = [0.; 4];
    let mut bc: [f64; 4] = [0.; 4];
    let mut cd: [f64; 4] = [0.; 4];
    let mut da: [f64; 4] = [0.; 4];
    let mut ac: [f64; 4] = [0.; 4];
    let mut bd: [f64; 4] = [0.; 4];
    let mut ab3: f64 = 0.;
    let mut bc3: f64 = 0.;
    let mut cd3: f64 = 0.;
    let mut da3: f64 = 0.;
    let mut ac3: f64 = 0.;
    let mut bd3: f64 = 0.;
    let mut abeps: f64 = 0.;
    let mut bceps: f64 = 0.;
    let mut cdeps: f64 = 0.;
    let mut daeps: f64 = 0.;
    let mut aceps: f64 = 0.;
    let mut bdeps: f64 = 0.;
    let mut temp8a: [f64; 8] = [0.; 8];
    let mut temp8b: [f64; 8] = [0.; 8];
    let mut temp8c: [f64; 8] = [0.; 8];
    let mut temp16: [f64; 16] = [0.; 16];
    let mut temp24: [f64; 24] = [0.; 24];
    let mut temp48: [f64; 48] = [0.; 48];
    let mut temp8alen: i32 = 0;
    let mut temp8blen: i32 = 0;
    let mut temp8clen: i32 = 0;
    let mut temp16len: i32 = 0;
    let mut temp24len: i32 = 0;
    let mut temp48len: i32 = 0;
    let mut xdet: [f64; 96] = [0.; 96];
    let mut ydet: [f64; 96] = [0.; 96];
    let mut zdet: [f64; 96] = [0.; 96];
    let mut xydet: [f64; 192] = [0.; 192];
    let mut xlen: i32 = 0;
    let mut ylen: i32 = 0;
    let mut zlen: i32 = 0;
    let mut xylen: i32 = 0;
    let mut adet: [f64; 288] = [0.; 288];
    let mut bdet: [f64; 288] = [0.; 288];
    let mut cdet: [f64; 288] = [0.; 288];
    let mut ddet: [f64; 288] = [0.; 288];
    let mut alen: i32 = 0;
    let mut blen: i32 = 0;
    let mut clen: i32 = 0;
    let mut dlen: i32 = 0;
    let mut abdet: [f64; 576] = [0.; 576];
    let mut cddet: [f64; 576] = [0.; 576];
    let mut ablen: i32 = 0;
    let mut cdlen: i32 = 0;
    let mut fin1: [f64; 1152] = [0.; 1152];
    let mut finlength: i32 = 0;
    let mut aextail: f64 = 0.;
    let mut bextail: f64 = 0.;
    let mut cextail: f64 = 0.;
    let mut dextail: f64 = 0.;
    let mut aeytail: f64 = 0.;
    let mut beytail: f64 = 0.;
    let mut ceytail: f64 = 0.;
    let mut deytail: f64 = 0.;
    let mut aeztail: f64 = 0.;
    let mut beztail: f64 = 0.;
    let mut ceztail: f64 = 0.;
    let mut deztail: f64 = 0.;
    let mut bvirt: f64 = 0.;
    let mut avirt: f64 = 0.;
    let mut bround: f64 = 0.;
    let mut around: f64 = 0.;
    let mut c: f64 = 0.;
    let mut abig: f64 = 0.;
    let mut ahi: f64 = 0.;
    let mut alo: f64 = 0.;
    let mut bhi: f64 = 0.;
    let mut blo: f64 = 0.;
    let mut err1: f64 = 0.;
    let mut err2: f64 = 0.;
    let mut err3: f64 = 0.;
    let mut _i: f64 = 0.;
    let mut _j: f64 = 0.;
    let mut _0: f64 = 0.;
    aex =
        *pa.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    bex =
        *pb.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    cex =
        *pc.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    dex =
        *pd.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    aey =
        *pa.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    bey =
        *pb.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    cey =
        *pc.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    dey =
        *pd.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    aez =
        *pa.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    bez =
        *pb.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    cez =
        *pc.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    dez =
        *pd.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    Two_Product(aex, bey, &mut aexbey1, &mut aexbey0);
    Two_Product(bex, aey, &mut bexaey1, &mut bexaey0);
    Two_Two_Diff(aexbey1, aexbey0, bexaey1, bexaey0, ab.as_mut_ptr());
    Two_Product(bex, cey, &mut bexcey1, &mut bexcey0);
    Two_Product(cex, bey, &mut cexbey1, &mut cexbey0);
    Two_Two_Diff(bexcey1, bexcey0, cexbey1, cexbey0, bc.as_mut_ptr());
    Two_Product(cex, dey, &mut cexdey1, &mut cexdey0);
    Two_Product(dex, cey, &mut dexcey1, &mut dexcey0);
    Two_Two_Diff(cexdey1, cexdey0, dexcey1, dexcey0, cd.as_mut_ptr());
    Two_Product(dex, aey, &mut dexaey1, &mut dexaey0);
    Two_Product(aex, dey, &mut aexdey1, &mut aexdey0);
    Two_Two_Diff(dexaey1, dexaey0, aexdey1, aexdey0, da.as_mut_ptr());
    Two_Product(aex, cey, &mut aexcey1, &mut aexcey0);
    Two_Product(cex, aey, &mut cexaey1, &mut cexaey0);
    Two_Two_Diff(aexcey1, aexcey0, cexaey1, cexaey0, ac.as_mut_ptr());
    Two_Product(bex, dey, &mut bexdey1, &mut bexdey0);
    Two_Product(dex, bey, &mut dexbey1, &mut dexbey0);
    Two_Two_Diff(bexdey1, bexdey0, dexbey1, dexbey0, bd.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, cd.as_mut_ptr(), bez,
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, bd.as_mut_ptr(), -cez,
                                 temp8b.as_mut_ptr());
    temp8clen =
        scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(), dez,
                                 temp8c.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp24len =
        fast_expansion_sum_zeroelim(temp8clen, temp8c.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), temp24.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), aex,
                                 temp48.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), -aex,
                                 xdet.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), aey,
                                 temp48.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), -aey,
                                 ydet.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), aez,
                                 temp48.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), -aez,
                                 zdet.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, xdet.as_mut_ptr(), ylen,
                                    ydet.as_mut_ptr(), xydet.as_mut_ptr());
    alen =
        fast_expansion_sum_zeroelim(xylen, xydet.as_mut_ptr(), zlen,
                                    zdet.as_mut_ptr(), adet.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, da.as_mut_ptr(), cez,
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, ac.as_mut_ptr(), dez,
                                 temp8b.as_mut_ptr());
    temp8clen =
        scale_expansion_zeroelim(4 as i32, cd.as_mut_ptr(), aez,
                                 temp8c.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp24len =
        fast_expansion_sum_zeroelim(temp8clen, temp8c.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), temp24.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), bex,
                                 temp48.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), bex,
                                 xdet.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), bey,
                                 temp48.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), bey,
                                 ydet.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), bez,
                                 temp48.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), bez,
                                 zdet.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, xdet.as_mut_ptr(), ylen,
                                    ydet.as_mut_ptr(), xydet.as_mut_ptr());
    blen =
        fast_expansion_sum_zeroelim(xylen, xydet.as_mut_ptr(), zlen,
                                    zdet.as_mut_ptr(), bdet.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(), dez,
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, bd.as_mut_ptr(), aez,
                                 temp8b.as_mut_ptr());
    temp8clen =
        scale_expansion_zeroelim(4 as i32, da.as_mut_ptr(), bez,
                                 temp8c.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp24len =
        fast_expansion_sum_zeroelim(temp8clen, temp8c.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), temp24.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), cex,
                                 temp48.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), -cex,
                                 xdet.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), cey,
                                 temp48.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), -cey,
                                 ydet.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), cez,
                                 temp48.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), -cez,
                                 zdet.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, xdet.as_mut_ptr(), ylen,
                                    ydet.as_mut_ptr(), xydet.as_mut_ptr());
    clen =
        fast_expansion_sum_zeroelim(xylen, xydet.as_mut_ptr(), zlen,
                                    zdet.as_mut_ptr(), cdet.as_mut_ptr());
    temp8alen =
        scale_expansion_zeroelim(4 as i32, bc.as_mut_ptr(), aez,
                                 temp8a.as_mut_ptr());
    temp8blen =
        scale_expansion_zeroelim(4 as i32, ac.as_mut_ptr(), -bez,
                                 temp8b.as_mut_ptr());
    temp8clen =
        scale_expansion_zeroelim(4 as i32, ab.as_mut_ptr(), cez,
                                 temp8c.as_mut_ptr());
    temp16len =
        fast_expansion_sum_zeroelim(temp8alen, temp8a.as_mut_ptr(), temp8blen,
                                    temp8b.as_mut_ptr(), temp16.as_mut_ptr());
    temp24len =
        fast_expansion_sum_zeroelim(temp8clen, temp8c.as_mut_ptr(), temp16len,
                                    temp16.as_mut_ptr(), temp24.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), dex,
                                 temp48.as_mut_ptr());
    xlen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), dex,
                                 xdet.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), dey,
                                 temp48.as_mut_ptr());
    ylen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), dey,
                                 ydet.as_mut_ptr());
    temp48len =
        scale_expansion_zeroelim(temp24len, temp24.as_mut_ptr(), dez,
                                 temp48.as_mut_ptr());
    zlen =
        scale_expansion_zeroelim(temp48len, temp48.as_mut_ptr(), dez,
                                 zdet.as_mut_ptr());
    xylen =
        fast_expansion_sum_zeroelim(xlen, xdet.as_mut_ptr(), ylen,
                                    ydet.as_mut_ptr(), xydet.as_mut_ptr());
    dlen =
        fast_expansion_sum_zeroelim(xylen, xydet.as_mut_ptr(), zlen,
                                    zdet.as_mut_ptr(), ddet.as_mut_ptr());
    ablen =
        fast_expansion_sum_zeroelim(alen, adet.as_mut_ptr(), blen,
                                    bdet.as_mut_ptr(), abdet.as_mut_ptr());
    cdlen =
        fast_expansion_sum_zeroelim(clen, cdet.as_mut_ptr(), dlen,
                                    ddet.as_mut_ptr(), cddet.as_mut_ptr());
    finlength =
        fast_expansion_sum_zeroelim(ablen, abdet.as_mut_ptr(), cdlen,
                                    cddet.as_mut_ptr(), fin1.as_mut_ptr());
    det = estimate(finlength, fin1.as_mut_ptr());
    errbound = PARAMS.isperrboundB * permanent;
    if det >= errbound || -det >= errbound { return det }
    aextail =
        Two_Diff_Tail(*pa.offset(0 as i32 as isize),
                      *pe.offset(0 as i32 as isize), aex);
    aeytail =
        Two_Diff_Tail(*pa.offset(1 as i32 as isize),
                      *pe.offset(1 as i32 as isize), aey);
    aeztail =
        Two_Diff_Tail(*pa.offset(2 as i32 as isize),
                      *pe.offset(2 as i32 as isize), aez);
    bextail =
        Two_Diff_Tail(*pb.offset(0 as i32 as isize),
                      *pe.offset(0 as i32 as isize), bex);
    beytail =
        Two_Diff_Tail(*pb.offset(1 as i32 as isize),
                      *pe.offset(1 as i32 as isize), bey);
    beztail =
        Two_Diff_Tail(*pb.offset(2 as i32 as isize),
                      *pe.offset(2 as i32 as isize), bez);
    cextail =
        Two_Diff_Tail(*pc.offset(0 as i32 as isize),
                      *pe.offset(0 as i32 as isize), cex);
    ceytail =
        Two_Diff_Tail(*pc.offset(1 as i32 as isize),
                      *pe.offset(1 as i32 as isize), cey);
    ceztail =
        Two_Diff_Tail(*pc.offset(2 as i32 as isize),
                      *pe.offset(2 as i32 as isize), cez);
    dextail =
        Two_Diff_Tail(*pd.offset(0 as i32 as isize),
                      *pe.offset(0 as i32 as isize), dex);
    deytail =
        Two_Diff_Tail(*pd.offset(1 as i32 as isize),
                      *pe.offset(1 as i32 as isize), dey);
    deztail =
        Two_Diff_Tail(*pd.offset(2 as i32 as isize),
                      *pe.offset(2 as i32 as isize), dez);
    if aextail == 0.0f64 && aeytail == 0.0f64 && aeztail == 0.0f64 &&
           bextail == 0.0f64 && beytail == 0.0f64 && beztail == 0.0f64 &&
           cextail == 0.0f64 && ceytail == 0.0f64 && ceztail == 0.0f64 &&
           dextail == 0.0f64 && deytail == 0.0f64 && deztail == 0.0f64 {
        return det
    }
    errbound = PARAMS.isperrboundC * permanent + PARAMS.resulterrbound * Absolute(det);
    abeps = aex * beytail + bey * aextail - (aey * bextail + bex * aeytail);
    bceps = bex * ceytail + cey * bextail - (bey * cextail + cex * beytail);
    cdeps = cex * deytail + dey * cextail - (cey * dextail + dex * ceytail);
    daeps = dex * aeytail + aey * dextail - (dey * aextail + aex * deytail);
    aceps = aex * ceytail + cey * aextail - (aey * cextail + cex * aeytail);
    bdeps = bex * deytail + dey * bextail - (bey * dextail + dex * beytail);
    det +=
        (bex * bex + bey * bey + bez * bez) *
            (cez * daeps + dez * aceps + aez * cdeps +
                 (ceztail * da3 + deztail * ac3 + aeztail * cd3)) +
            (dex * dex + dey * dey + dez * dez) *
                (aez * bceps - bez * aceps + cez * abeps +
                     (aeztail * bc3 - beztail * ac3 + ceztail * ab3)) -
            ((aex * aex + aey * aey + aez * aez) *
                 (bez * cdeps - cez * bdeps + dez * bceps +
                      (beztail * cd3 - ceztail * bd3 + deztail * bc3)) +
                 (cex * cex + cey * cey + cez * cez) *
                     (dez * abeps + aez * bdeps + bez * daeps +
                          (deztail * ab3 + aeztail * bd3 + beztail * da3))) +
            2.0f64 *
                ((bex * bextail + bey * beytail + bez * beztail) *
                     (cez * da3 + dez * ac3 + aez * cd3) +
                     (dex * dextail + dey * deytail + dez * deztail) *
                         (aez * bc3 - bez * ac3 + cez * ab3) -
                     ((aex * aextail + aey * aeytail + aez * aeztail) *
                          (bez * cd3 - cez * bd3 + dez * bc3) +
                          (cex * cextail + cey * ceytail + cez * ceztail) *
                              (dez * ab3 + aez * bd3 + bez * da3)));
    if det >= errbound || -det >= errbound { return det }
    return insphereexact(pa, pb, pc, pd, pe);
}

pub unsafe fn insphere(mut pa: *const f64,
                                  mut pb: *const f64,
                                  mut pc: *const f64,
                                  mut pd: *const f64,
                                  mut pe: *const f64)
 -> f64 {
    let mut aex: f64 = 0.;
    let mut bex: f64 = 0.;
    let mut cex: f64 = 0.;
    let mut dex: f64 = 0.;
    let mut aey: f64 = 0.;
    let mut bey: f64 = 0.;
    let mut cey: f64 = 0.;
    let mut dey: f64 = 0.;
    let mut aez: f64 = 0.;
    let mut bez: f64 = 0.;
    let mut cez: f64 = 0.;
    let mut dez: f64 = 0.;
    let mut aexbey: f64 = 0.;
    let mut bexaey: f64 = 0.;
    let mut bexcey: f64 = 0.;
    let mut cexbey: f64 = 0.;
    let mut cexdey: f64 = 0.;
    let mut dexcey: f64 = 0.;
    let mut dexaey: f64 = 0.;
    let mut aexdey: f64 = 0.;
    let mut aexcey: f64 = 0.;
    let mut cexaey: f64 = 0.;
    let mut bexdey: f64 = 0.;
    let mut dexbey: f64 = 0.;
    let mut alift: f64 = 0.;
    let mut blift: f64 = 0.;
    let mut clift: f64 = 0.;
    let mut dlift: f64 = 0.;
    let mut ab: f64 = 0.;
    let mut bc: f64 = 0.;
    let mut cd: f64 = 0.;
    let mut da: f64 = 0.;
    let mut ac: f64 = 0.;
    let mut bd: f64 = 0.;
    let mut abc: f64 = 0.;
    let mut bcd: f64 = 0.;
    let mut cda: f64 = 0.;
    let mut dab: f64 = 0.;
    let mut aezplus: f64 = 0.;
    let mut bezplus: f64 = 0.;
    let mut cezplus: f64 = 0.;
    let mut dezplus: f64 = 0.;
    let mut aexbeyplus: f64 = 0.;
    let mut bexaeyplus: f64 = 0.;
    let mut bexceyplus: f64 = 0.;
    let mut cexbeyplus: f64 = 0.;
    let mut cexdeyplus: f64 = 0.;
    let mut dexceyplus: f64 = 0.;
    let mut dexaeyplus: f64 = 0.;
    let mut aexdeyplus: f64 = 0.;
    let mut aexceyplus: f64 = 0.;
    let mut cexaeyplus: f64 = 0.;
    let mut bexdeyplus: f64 = 0.;
    let mut dexbeyplus: f64 = 0.;
    let mut det: f64 = 0.;
    let mut permanent: f64 = 0.;
    let mut errbound: f64 = 0.;
    aex =
        *pa.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    bex =
        *pb.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    cex =
        *pc.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    dex =
        *pd.offset(0 as i32 as isize) -
            *pe.offset(0 as i32 as isize);
    aey =
        *pa.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    bey =
        *pb.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    cey =
        *pc.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    dey =
        *pd.offset(1 as i32 as isize) -
            *pe.offset(1 as i32 as isize);
    aez =
        *pa.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    bez =
        *pb.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    cez =
        *pc.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    dez =
        *pd.offset(2 as i32 as isize) -
            *pe.offset(2 as i32 as isize);
    aexbey = aex * bey;
    bexaey = bex * aey;
    ab = aexbey - bexaey;
    bexcey = bex * cey;
    cexbey = cex * bey;
    bc = bexcey - cexbey;
    cexdey = cex * dey;
    dexcey = dex * cey;
    cd = cexdey - dexcey;
    dexaey = dex * aey;
    aexdey = aex * dey;
    da = dexaey - aexdey;
    aexcey = aex * cey;
    cexaey = cex * aey;
    ac = aexcey - cexaey;
    bexdey = bex * dey;
    dexbey = dex * bey;
    bd = bexdey - dexbey;
    abc = aez * bc - bez * ac + cez * ab;
    bcd = bez * cd - cez * bd + dez * bc;
    cda = cez * da + dez * ac + aez * cd;
    dab = dez * ab + aez * bd + bez * da;
    alift = aex * aex + aey * aey + aez * aez;
    blift = bex * bex + bey * bey + bez * bez;
    clift = cex * cex + cey * cey + cez * cez;
    dlift = dex * dex + dey * dey + dez * dez;
    det = dlift * abc - clift * dab + (blift * cda - alift * bcd);
    aezplus = Absolute(aez);
    bezplus = Absolute(bez);
    cezplus = Absolute(cez);
    dezplus = Absolute(dez);
    aexbeyplus = Absolute(aexbey);
    bexaeyplus = Absolute(bexaey);
    bexceyplus = Absolute(bexcey);
    cexbeyplus = Absolute(cexbey);
    cexdeyplus = Absolute(cexdey);
    dexceyplus = Absolute(dexcey);
    dexaeyplus = Absolute(dexaey);
    aexdeyplus = Absolute(aexdey);
    aexceyplus = Absolute(aexcey);
    cexaeyplus = Absolute(cexaey);
    bexdeyplus = Absolute(bexdey);
    dexbeyplus = Absolute(dexbey);
    permanent =
        ((cexdeyplus + dexceyplus) * bezplus +
             (dexbeyplus + bexdeyplus) * cezplus +
             (bexceyplus + cexbeyplus) * dezplus) * alift +
            ((dexaeyplus + aexdeyplus) * cezplus +
                 (aexceyplus + cexaeyplus) * dezplus +
                 (cexdeyplus + dexceyplus) * aezplus) * blift +
            ((aexbeyplus + bexaeyplus) * dezplus +
                 (bexdeyplus + dexbeyplus) * aezplus +
                 (dexaeyplus + aexdeyplus) * bezplus) * clift +
            ((bexceyplus + cexbeyplus) * aezplus +
                 (cexaeyplus + aexceyplus) * bezplus +
                 (aexbeyplus + bexaeyplus) * cezplus) * dlift;
    errbound = PARAMS.isperrboundA * permanent;
    if det > errbound || -det > errbound { return det }
    return insphereadapt(pa, pb, pc, pd, pe, permanent);
}
