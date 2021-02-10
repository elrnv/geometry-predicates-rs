/// Returns the absolute value of the given number.
///
/// This function exists since [`std::f64::abs`](std::f64::abs) is not available in core.
/// See [#50145](https://github.com/rust-lang/rust/issues/50145)
///
/// This implementation is identical to [`std::f64::abs`](std::f64::abs) on x86 but not on ARM at the time of this writing.
#[inline]
pub fn abs(a: f64) -> f64 {
    f64::from_bits(a.to_bits() & 0x7FFF_FFFF_FFFF_FFFF)
}

#[derive(Debug)]
struct PredicateParams {
    // Used to split floats in half.
    splitter: f64, // = 2^ceiling(p / 2) + 1.
    /* A set of coefficients used to calculate maximum roundoff errors.          */
    resulterrbound: f64,
    ccwerrbound_a: f64,
    ccwerrbound_b: f64,
    ccwerrbound_c: f64,
    o3derrbound_a: f64,
    o3derrbound_b: f64,
    o3derrbound_c: f64,
    iccerrbound_a: f64,
    iccerrbound_b: f64,
    iccerrbound_c: f64,
    isperrbound_a: f64,
    isperrbound_b: f64,
    isperrbound_c: f64,
}

// EPSILON and PARAMS.slitter were pregenerated using exactinit on a machine with IEEE 754 floats.
// See `exactinit` function below for details.

/// The largest power of two such that 1.0 + epsilon = 1.0 in floating-point
/// arithmetic.
///
/// This number bounds the relative roundoff error. It is used for
/// floating-point error analysis.
const EPSILON: f64 = 0.000_000_000_000_000_111_022_302_462_515_65;

///  Constants used in exact arithmetic.
///
///  See exactinit() for the function used to generate these values.
const PARAMS: PredicateParams = PredicateParams {
    ///  Used to split floating-point numbers into two half-length significands
    ///  for exact multiplication.
    splitter: 134_217_729f64,
    resulterrbound: (3.0 + 8.0 * EPSILON) * EPSILON,
    ccwerrbound_a: (3.0 + 16.0 * EPSILON) * EPSILON,
    ccwerrbound_b: (2.0 + 12.0 * EPSILON) * EPSILON,
    ccwerrbound_c: (9.0 + 64.0 * EPSILON) * EPSILON * EPSILON,
    o3derrbound_a: (7.0f64 + 56.0f64 * EPSILON) * EPSILON,
    o3derrbound_b: (3.0f64 + 28.0f64 * EPSILON) * EPSILON,
    o3derrbound_c: (26.0f64 + 288.0f64 * EPSILON) * EPSILON * EPSILON,
    iccerrbound_a: (10.0 + 96.0 * EPSILON) * EPSILON,
    iccerrbound_b: (4.0 + 48.0 * EPSILON) * EPSILON,
    iccerrbound_c: (44.0 + 576.0 * EPSILON) * EPSILON * EPSILON,
    isperrbound_a: (16.0f64 + 224.0f64 * EPSILON) * EPSILON,
    isperrbound_b: (5.0f64 + 72.0f64 * EPSILON) * EPSILON,
    isperrbound_c: (71.0f64 + 1408.0f64 * EPSILON) * EPSILON * EPSILON,
};

/* ****************************************************************************/
/*                                                                           */
/*  exactinit()   Initialize the variables used for exact arithmetic.        */
/*                                                                           */
/*  `epsilon' is the largest power of two such that 1.0 + epsilon = 1.0 in   */
/*  floating-point arithmetic.  `epsilon' bounds the relative roundoff       */
/*  error.  It is used for floating-point error analysis.                    */
/*                                                                           */
/*  `splitter' is used to split floating-point numbers into two half-        */
/*  length significands for exact multiplication.                            */
/*                                                                           */
/*  I imagine that a highly optimizing compiler might be too smart for its   */
/*  own good, and somehow cause this routine to fail, if it pretends that    */
/*  floating-point arithmetic is too much like real arithmetic.              */
/*                                                                           */
/*  Don't change this routine unless you fully understand it.                */
/*                                                                           */
/* ****************************************************************************/
#[allow(dead_code)] // This function is for reference only.
fn exactinit() -> PredicateParams {
    let mut check = 1.0_f64;
    let mut lastcheck;
    let mut every_other = 1_i32;
    let mut epsilon = 1.0f64;
    let mut splitter = 1.0f64;
    loop {
        /* Repeatedly divide `epsilon' by two until it is too small to add to    */
        /*   one without causing roundoff.  (Also check if the sum is equal to   */
        /*   the previous sum, for machines that round up instead of using exact */
        /*   rounding.  Not that this library will work on such machines anyway. */
        lastcheck = check;
        epsilon *= 0.5;
        if every_other != 0 {
            splitter *= 2.0f64
        }
        every_other = (every_other == 0) as i32;
        check = 1.0f64 + epsilon;
        if !(check != 1.0f64 && check != lastcheck) {
            break;
        }
    }
    splitter += 1.0f64;
    PredicateParams {
        splitter,
        /* Error bounds for orientation and incircle tests. */
        resulterrbound: (3.0f64 + 8.0f64 * epsilon) * epsilon,
        ccwerrbound_a: (3.0f64 + 16.0f64 * epsilon) * epsilon,
        ccwerrbound_b: (2.0f64 + 12.0f64 * epsilon) * epsilon,
        ccwerrbound_c: (9.0f64 + 64.0f64 * epsilon) * epsilon * epsilon,
        o3derrbound_a: (7.0f64 + 56.0f64 * epsilon) * epsilon,
        o3derrbound_b: (3.0f64 + 28.0f64 * epsilon) * epsilon,
        o3derrbound_c: (26.0f64 + 288.0f64 * epsilon) * epsilon * epsilon,
        iccerrbound_a: (10.0f64 + 96.0f64 * epsilon) * epsilon,
        iccerrbound_b: (4.0f64 + 48.0f64 * epsilon) * epsilon,
        iccerrbound_c: (44.0f64 + 576.0f64 * epsilon) * epsilon * epsilon,
        isperrbound_a: (16.0f64 + 224.0f64 * epsilon) * epsilon,
        isperrbound_b: (5.0f64 + 72.0f64 * epsilon) * epsilon,
        isperrbound_c: (71.0f64 + 1408.0f64 * epsilon) * epsilon * epsilon,
    }
}

/* Many of the operations are broken up into two pieces, a main part that    */
/*   performs an approximate operation, and a "tail" that computes the       */
/*   roundoff error of that operation.                                       */

#[inline]
pub fn fast_two_sum_tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt: f64 = x - a;
    b - bvirt
}

#[inline]
pub fn fast_two_sum(a: f64, b: f64) -> [f64; 2] {
    let x: f64 = a + b;
    [fast_two_sum_tail(a, b, x), x]
}

#[inline]
pub fn fast_two_diff_tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt: f64 = a - x;
    return bvirt - b;
}

#[inline]
pub fn fast_two_diff(a: f64, b: f64) -> [f64; 2] {
    let x: f64 = a - b;
    [fast_two_diff_tail(a, b, x), x]
}

#[inline]
pub fn two_sum_tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt: f64 = x - a;
    let avirt: f64 = x - bvirt;
    let bround: f64 = b - bvirt;
    let around: f64 = a - avirt;
    around + bround
}

#[inline]
pub fn two_sum(a: f64, b: f64) -> [f64; 2] {
    let x: f64 = a + b;
    [two_sum_tail(a, b, x), x]
}

#[inline]
pub fn two_diff_tail(a: f64, b: f64, x: f64) -> f64 {
    let bvirt: f64 = a - x;
    let avirt: f64 = x + bvirt;
    let bround: f64 = bvirt - b;
    let around: f64 = a - avirt;
    around + bround
}

#[inline]
pub fn two_diff(a: f64, b: f64) -> [f64; 2] {
    let x: f64 = a - b;
    [two_diff_tail(a, b, x), x]
}

#[inline]
pub fn split(a: f64) -> [f64; 2] {
    let c: f64 = PARAMS.splitter * a;
    let abig: f64 = c - a;
    let ahi = c - abig;
    let alo = a - ahi;
    [alo, ahi]
}

#[inline]
pub fn two_product_tail(a: f64, b: f64, x: f64) -> f64 {
    let [alo, ahi] = split(a);
    let [blo, bhi] = split(b);
    let err1: f64 = x - ahi * bhi;
    let err2: f64 = err1 - alo * bhi;
    let err3: f64 = err2 - ahi * blo;
    alo * blo - err3
}

#[inline]
pub fn two_product(a: f64, b: f64) -> [f64; 2] {
    let x = a * b;
    [two_product_tail(a, b, x), x]
}

/// Same as [`two_product`] where one of the inputs has
/// already been split.
///
/// Avoids redundant splitting.
#[inline]
pub fn two_product_presplit(a: f64, b: f64, bhi: f64, blo: f64) -> [f64; 2] {
    let x = a * b;
    let [alo, ahi] = split(a);
    let err1: f64 = x - ahi * bhi;
    let err2: f64 = err1 - alo * bhi;
    let err3: f64 = err2 - ahi * blo;
    [alo * blo - err3, x]
}

/// Same as [`two_product`] where both of the inputs have
/// already been split.
///
/// Avoids redundant splitting.
#[inline]
pub fn two_product_2presplit(a: f64, ahi: f64, alo: f64, b: f64, bhi: f64, blo: f64) -> [f64; 2] {
    let x = a * b;
    let err1: f64 = x - ahi * bhi;
    let err2: f64 = err1 - alo * bhi;
    let err3: f64 = err2 - ahi * blo;
    [alo * blo - err3, x]
}

#[inline]
pub fn square_tail(a: f64, x: f64) -> f64 {
    let [alo, ahi] = split(a);
    let err1: f64 = x - ahi * ahi;
    let err3: f64 = err1 - (ahi + ahi) * alo;
    alo * alo - err3
}

/// Squaring can be done more quickly than [`two_product`].
#[inline]
pub fn square(a: f64) -> [f64; 2] {
    let x = a * a;
    [square_tail(a, x), x]
}

// Macros for summing expansions of various fixed lengths.  These are all
// unrolled versions of expansion_sum().

#[inline]
pub fn two_one_sum(a1: f64, a0: f64, b: f64) -> [f64; 3] {
    let [x0, _i] = two_sum(a0, b);
    let [x1, x2] = two_sum(a1, _i);
    [x0, x1, x2]
}

#[inline]
pub fn two_one_diff(a1: f64, a0: f64, b: f64) -> [f64; 3] {
    let [x0, _i] = two_diff(a0, b);
    let [x1, x2] = two_sum(a1, _i);
    [x2, x1, x0]
}

#[inline]
pub fn two_two_sum(a1: f64, a0: f64, b1: f64, b0: f64) -> [f64; 4] {
    let [x0, _0, _j] = two_one_sum(a1, a0, b0);
    let [x1, x2, x3] = two_one_sum(_j, _0, b1);
    [x0, x1, x2, x3]
}

#[inline]
pub fn two_two_diff(a1: f64, a0: f64, b1: f64, b0: f64) -> [f64; 4] {
    let [_j, _0, x0] = two_one_diff(a1, a0, b0);
    let [x3, x2, x1] = two_one_diff(_j, _0, b1);
    [x0, x1, x2, x3]
}

#[inline]
pub fn four_one_sum(a3: f64, a2: f64, a1: f64, a0: f64, b: f64) -> [f64; 5] {
    let [x0, x1, _j] = two_one_sum(a1, a0, b);
    let [x2, x3, x4] = two_one_sum(a3, a2, _j);
    [x0, x1, x2, x3, x4]
}

#[inline]
pub fn four_two_sum(a3: f64, a2: f64, a1: f64, a0: f64, b1: f64, b0: f64) -> [f64; 6] {
    let [x0, _0, _1, _2, _k] = four_one_sum(a3, a2, a1, a0, b0);
    let [x1, x2, x3, x4, x5] = four_one_sum(_k, _2, _1, _0, b1);
    [x0, x1, x2, x3, x4, x5]
}

#[inline]
pub fn four_four_sum(
    a3: f64,
    a2: f64,
    a1: f64,
    a0: f64,
    b4: f64,
    b3: f64,
    b1: f64,
    b0: f64,
) -> [f64; 8] {
    let [x0, x1, _0, _1, _2, _l] = four_two_sum(a3, a2, a1, a0, b1, b0);
    let [x2, x3, x4, x5, x6, x7] = four_two_sum(_l, _2, _1, _0, b4, b3);
    [x7, x6, x5, x4, x3, x2, x1, x0]
}

#[inline]
pub fn eight_one_sum(
    a7: f64,
    a6: f64,
    a5: f64,
    a4: f64,
    a3: f64,
    a2: f64,
    a1: f64,
    a0: f64,
    b: f64,
) -> [f64; 9] {
    let [x0, x1, x2, x3, _j] = four_one_sum(a3, a2, a1, a0, b);
    let [x4, x5, x6, x7, x8] = four_one_sum(a7, a6, a5, a4, _j);
    [x0, x1, x2, x3, x4, x5, x6, x7, x8]
}

#[inline]
pub fn eight_two_sum(
    a7: f64,
    a6: f64,
    a5: f64,
    a4: f64,
    a3: f64,
    a2: f64,
    a1: f64,
    a0: f64,
    b1: f64,
    b0: f64,
) -> [f64; 10] {
    let [x0, _0, _1, _2, _3, _4, _5, _6, _k] = eight_one_sum(a7, a6, a5, a4, a3, a2, a1, a0, b0);
    let [x1, x2, x3, x4, x5, x6, x7, x8, x9] = eight_one_sum(_k, _6, _5, _4, _3, _2, _1, _0, b1);
    [x0, x1, x2, x3, x4, x5, x6, x7, x8, x9]
}

#[inline]
pub fn eight_four_sum(
    a7: f64,
    a6: f64,
    a5: f64,
    a4: f64,
    a3: f64,
    a2: f64,
    a1: f64,
    a0: f64,
    b4: f64,
    b3: f64,
    b1: f64,
    b0: f64,
) -> [f64; 12] {
    let [x0, x1, _0, _1, _2, _3, _4, _5, _6, _l] =
        eight_two_sum(a7, a6, a5, a4, a3, a2, a1, a0, b1, b0);
    let [x2, x3, x4, x5, x6, x7, x8, x9, x10, x11] =
        eight_two_sum(_l, _6, _5, _4, _3, _2, _1, _0, b4, b3);
    [x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11]
}

// Macros for multiplying expansions of various fixed lengths.

#[inline]
pub fn two_one_product(a1: f64, a0: f64, b: f64) -> [f64; 4] {
    let [blo, bhi] = split(b);
    let [x0, _i] = two_product_presplit(a0, b, bhi, blo);
    let [_0, _j] = two_product_presplit(a1, b, bhi, blo);
    let [x1, _k] = two_sum(_i, _0);
    let [x2, x3] = fast_two_sum(_j, _k);
    [x0, x1, x2, x3]
}

#[inline]
pub fn four_one_product(a3: f64, a2: f64, a1: f64, a0: f64, b: f64) -> [f64; 8] {
    let [blo, bhi] = split(b);
    let [x0, _i] = two_product_presplit(a0, b, bhi, blo);
    let [_0, _j] = two_product_presplit(a1, b, bhi, blo);
    let [x1, _k] = two_sum(_i, _0);
    let [x2, _i] = fast_two_sum(_j, _k);
    let [_0, _j] = two_product_presplit(a2, b, bhi, blo);
    let [x3, _k] = two_sum(_i, _0);
    let [x4, _i] = fast_two_sum(_j, _k);
    let [_0, _j] = two_product_presplit(a3, b, bhi, blo);
    let [x5, _k] = two_sum(_i, _0);
    let [x6, x7] = fast_two_sum(_j, _k);
    [x0, x1, x2, x3, x4, x5, x6, x7]
}

#[inline]
pub fn two_two_product(a1: f64, a0: f64, b1: f64, b0: f64) -> [f64; 8] {
    let [a0lo, a0hi] = split(a0);
    let [blo, bhi] = split(b0);
    let [x0, _i] = two_product_2presplit(a0, a0hi, a0lo, b0, bhi, blo);
    let [a1lo, a1hi] = split(a1);
    let [_0, _j] = two_product_2presplit(a1, a1hi, a1lo, b0, bhi, blo);
    let [_1, _k] = two_sum(_i, _0);
    let [_2, _l] = fast_two_sum(_j, _k);
    let [blo, bhi] = split(b1);
    let [_0, _i] = two_product_2presplit(a0, a0hi, a0lo, b1, bhi, blo);
    let [x1, _k] = two_sum(_1, _0);
    let [_1, _j] = two_sum(_2, _k);
    let [_2, _m] = two_sum(_l, _j);
    let [_0, _j] = two_product_2presplit(a1, a1hi, a1lo, b1, bhi, blo);
    let [_0, _n] = two_sum(_i, _0);
    let [x2, _i] = two_sum(_1, _0);
    let [_1, _k] = two_sum(_2, _i);
    let [_2, _l] = two_sum(_m, _k);
    let [_0, _k] = two_sum(_j, _n);
    let [x3, _j] = two_sum(_1, _0);
    let [_1, _i] = two_sum(_2, _j);
    let [_2, _m] = two_sum(_l, _i);
    let [x4, _i] = two_sum(_1, _k);
    let [x5, _k] = two_sum(_2, _i);
    let [x6, x7] = two_sum(_m, _k);
    [x0, x1, x2, x3, x4, x5, x6, x7]
}

// An expansion of length two can be squared more quickly than finding the
// product of two different expansions of length two, and the result is
// guaranteed to have no more than six (rather than eight) components.

#[inline]
pub fn two_square(a1: f64, a0: f64) -> [f64; 6] {
    let [x0, _j] = square(a0);
    let _0: f64 = a0 + a0;
    let [_1, _k] = two_product(a1, _0);
    let [x1, _2, _l] = two_one_sum(_k, _1, _j);
    let [_1, _j] = square(a1);
    let [x2, x3, x4, x5] = two_two_sum(_j, _1, _l, _2);
    [x0, x1, x2, x3, x4, x5]
}

///  Adds a scalar to an expansion.
///
///  Sets `h = e + b`.  See [the paper](http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf) for details.
///
///  Maintains the nonoverlapping property.  If round-to-even is used (as
///  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent
///  properties as well.  (That is, if `e` has one of these properties, so
///  will `h`.)
#[inline]
pub fn grow_expansion(e: &[f64], b: f64, h: &mut [f64]) -> usize {
    let mut q = b;
    let mut eindex = 0;
    while eindex < e.len() {
        let [hnew, q_new] = two_sum(q, e[eindex]);
        q = q_new;
        h[eindex] = hnew;
        eindex += 1;
    }
    h[eindex] = q;
    eindex + 1
}

///  Adds a scalar to an expansion, eliminating zero components from the output
///  expansion.
///                                                                        
///  Sets `h = e + b`. See [the paper](http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf) for details.        
///                                                                        
///  Maintains the nonoverlapping property.  If round-to-even is used (as  
///  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent
///  properties as well.  (That is, if `e` has one of these properties, so   
///  will `h`.)
#[inline]
pub fn grow_expansion_zeroelim(e: &[f64], b: f64, h: &mut [f64]) -> usize {
    let mut hindex = 0;
    let mut q = b;
    let mut eindex = 0;
    while eindex < e.len() {
        let [hh, q_new] = two_sum(q, e[eindex]);
        q = q_new;
        if hh != 0.0f64 {
            let fresh0 = hindex;
            hindex = hindex + 1;
            h[fresh0] = hh;
        }
        eindex += 1
    }
    if q != 0.0f64 || hindex == 0 {
        let fresh1 = hindex;
        hindex = hindex + 1;
        h[fresh1] = q;
    }
    hindex
}

///  Sums two expansions.
///                      
///  Sets `h = e + f`. See [the paper](http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf) for details.
///                                                                 
///  Maintains the nonoverlapping property.  If round-to-even is used (as
///  with IEEE 754), maintains the nonadjacent property as well.  (That is,
///  if `e` has one of these properties, so will `h`.)  Does NOT maintain the
///  strongly nonoverlapping property.                                  
#[inline]
pub fn expansion_sum(e: &[f64], f: &[f64], h: &mut [f64]) -> usize {
    let mut q = f[0];
    let mut hindex = 0;
    while hindex < e.len() {
        let [hh, qnew] = two_sum(q, e[hindex]);
        h[hindex] = hh;
        q = qnew;
        hindex += 1
    }
    h[hindex] = q;
    let mut hlast = hindex;
    let mut findex = 1;
    while findex < f.len() {
        q = f[findex];
        hindex = findex;
        while hindex <= hlast {
            let [hh, qnew] = two_sum(q, h[hindex]);
            h[hindex] = hh;
            q = qnew;
            hindex += 1
        }
        hlast += 1;
        h[hlast] = q;
        findex += 1
    }
    hlast + 1
}

///  Sums two expansions, eliminating zero components from the output expansion.
///                                                                          
///  Sets `h = e + f`. See [the
///  paper](http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf) for details.
///                                                                        
///  Maintains the nonoverlapping property.  If round-to-even is used (as  
///  with IEEE 754), maintains the nonadjacent property as well.  (That is,
///  if `e` has one of these properties, so will `h`.)  Does NOT maintain the  
///  strongly nonoverlapping property.                                     
#[inline]
pub fn expansion_sum_zeroelim1(e: &[f64], f: &[f64], h: &mut [f64]) -> usize {
    let mut q = f[0];
    let mut hindex = 0;
    while hindex < e.len() {
        let [hh, qnew] = two_sum(q, e[hindex]);
        h[hindex] = hh;
        q = qnew;
        hindex += 1
    }
    h[hindex] = q;
    let mut hlast = hindex;
    let mut findex = 1;
    while findex < f.len() {
        q = f[findex];
        hindex = findex;
        while hindex <= hlast {
            let [hh, qnew] = two_sum(q, h[hindex]);
            h[hindex] = hh;
            q = qnew;
            hindex += 1
        }
        hlast += 1;
        h[hlast] = q;
        findex += 1
    }
    let mut hindex: isize = -1;
    let mut index = 0;
    while index <= hlast {
        let hnow = h[index];
        if hnow != 0.0 {
            hindex += 1;
            h[hindex as usize] = hnow
        }
        index += 1
    }
    if hindex == -1 {
        1
    } else {
        hindex as usize + 1
    }
}
///  Sums two expansions, eliminating zero components from the output expansion.
///                                                                        
///  Sets `h = e + f`. See [the
///  paper](http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf) for details.
///                                                                        
///  Maintains the nonoverlapping property.  If round-to-even is used (as  
///  with IEEE 754), maintains the nonadjacent property as well.  (That is,
///  if `e` has one of these properties, so will `h`.)  Does NOT maintain the  
///  strongly nonoverlapping property.                                     
#[inline]
pub fn expansion_sum_zeroelim2(e: &[f64], f: &[f64], h: &mut [f64]) -> usize {
    let mut hindex = 0;
    let mut q = f[0];
    let mut eindex = 0;
    while eindex < e.len() {
        let [hh, qnew] = two_sum(q, e[eindex]);
        q = qnew;
        if hh != 0.0 {
            h[hindex] = hh;
            hindex += 1;
        }
        eindex += 1
    }
    h[hindex] = q;
    let mut hlast = hindex;
    let mut findex = 1;
    while findex < f.len() {
        hindex = 0;
        q = f[findex];
        eindex = 0;
        while eindex <= hlast {
            let [hh, qnew] = two_sum(q, h[eindex]);
            q = qnew;
            if hh != 0.0 {
                h[hindex] = hh;
                hindex += 1;
            }
            eindex += 1
        }
        h[hindex] = q;
        hlast = hindex;
        findex += 1
    }
    hlast + 1
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
/*
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
    enow = e[0];
    fnow = f[0];
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
            let [y, x] = Fast_Two_Sum(enow, Q);
            Qnew = x;
            h[0] = y;
            eindex += 1;
            enow = *e.offset(eindex as isize)
        } else {
            let [y, x] = Fast_Two_Sum(fnow, Q);
            Qnew = x;
            h[0] = y;
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
*/

///  Sums two expansions, eliminating zero components from the output expansion.
///                                                                        
///  Sets `h = e + f`.  See [the
///  paper](http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf) for details.
///                                                                        
///  If round-to-even is used (as with IEEE 754), maintains the strongly   
///  nonoverlapping property.  (That is, if `e` is strongly nonoverlapping, `h`
///  will be also.)  Does NOT maintain the nonoverlapping or nonadjacent   
///  properties.                                                           
#[inline]
pub fn fast_expansion_sum_zeroelim(e: &[f64], f: &[f64], h: &mut [f64]) -> usize {
    let mut q;
    let mut findex = 0;
    let mut eindex = findex;

    let enow = e[0];
    let fnow = f[0];
    if (fnow > enow) == (fnow > -enow) {
        q = enow;
        eindex += 1;
    } else {
        q = fnow;
        findex += 1;
    }

    let mut hindex = 0;
    if eindex < e.len() && findex < f.len() {
        let enow = e[eindex];
        let fnow = f[findex];
        let [hh, q_new] = if (fnow > enow) == (fnow > -enow) {
            eindex += 1;
            fast_two_sum(enow, q)
        } else {
            findex += 1;
            fast_two_sum(fnow, q)
        };
        q = q_new;
        if hh != 0.0f64 {
            h[hindex] = hh;
            hindex += 1;
        }
        while eindex < e.len() && findex < f.len() {
            let enow = e[eindex];
            let fnow = f[findex];
            let [hh, q_new] = if (fnow > enow) == (fnow > -enow) {
                eindex += 1;
                two_sum(q, enow)
            } else {
                findex += 1;
                two_sum(q, fnow)
            };
            q = q_new;
            if hh != 0.0f64 {
                h[hindex] = hh;
                hindex += 1;
            }
        }
    }
    while eindex < e.len() {
        let [hh, q_new] = two_sum(q, e[eindex]);
        eindex += 1;
        q = q_new;
        if hh != 0.0f64 {
            h[hindex] = hh;
            hindex += 1;
        }
    }
    while findex < f.len() {
        let [hh, q_new] = two_sum(q, f[findex]);
        findex += 1;
        q = q_new;
        if hh != 0.0f64 {
            h[hindex] = hh;
            hindex += 1;
        }
    }
    if q != 0.0f64 || hindex == 0 {
        h[hindex] = q;
        hindex += 1;
    }
    hindex
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
/*
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
    enow = e[0];
    fnow = f[0];
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
        let [y, x] = Fast_Two_Sum(enow, g0);
        Qnew = x;
        q = y;
        eindex += 1;
        enow = *e.offset(eindex as isize)
    } else {
        let [y, x] = Fast_Two_Sum(fnow, g0);
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
            let [y, x] = Fast_Two_Sum(enow, q);
            R = x;
            *h.offset(hindex as isize) = y;
            eindex += 1;
            enow = *e.offset(eindex as isize)
        } else {
            let [y, x] = Fast_Two_Sum(fnow, q);
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
*/
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

/*
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
    enow = e[0];
    fnow = f[0];
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
        let [y, x] = Fast_Two_Sum(enow, g0);
        Qnew = x;
        q = y;
        eindex += 1;
        enow = *e.offset(eindex as isize)
    } else {
        let [y, x] = Fast_Two_Sum(fnow, g0);
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
            let [y, x] = Fast_Two_Sum(enow, q);
            R = x;
            hh = y;
            eindex += 1;
            enow = *e.offset(eindex as isize)
        } else {
            let [y, x] = Fast_Two_Sum(fnow, q);
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
*/
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
/*
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
    Two_Product_Presplit(e[0], b, bhi, blo,
                         &mut Q, &mut h[0]);
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
*/

///  Multiply an expansion by a scalar, eliminating zero components from the
///  output expansion.
///                                                                       
///  Sets `h = be`. See either [\[1\]] or [\[2\]] for details.
///                                                                       
///  Maintains the nonoverlapping property.  If round-to-even is used (as
///  with IEEE 754), maintains the strongly nonoverlapping and nonadjacent
///  properties as well.  (That is, if `e` has one of these properties, so  
///  will `h`.)                                                             
///
/// [\[1\]]: http://www.cs.berkeley.edu/~jrs/papers/robustr.pdf
/// [\[2\]]: http://www.cs.berkeley.edu/~jrs/papers/robust-predicates.pdf
pub fn scale_expansion_zeroelim(e: &[f64], b: f64, h: &mut [f64]) -> usize {
    let [blo, bhi] = split(b);
    let [hh, mut q] = two_product_presplit(e[0], b, bhi, blo);

    let mut hindex = 0;
    if hh != 0.0f64 {
        h[hindex] = hh;
        hindex += 1;
    }
    for &enow in e.iter().skip(1) {
        let [product0, product1] = two_product_presplit(enow, b, bhi, blo);
        let [hh, sum] = two_sum(q, product0);
        if hh != 0.0f64 {
            h[hindex] = hh;
            hindex += 1;
        }
        let [hh, q_new] = fast_two_sum(product1, sum);
        q = q_new;
        if hh != 0.0f64 {
            h[hindex] = hh;
            hindex += 1;
        }
    }
    if q != 0.0f64 || hindex == 0 {
        h[hindex] = q;
        hindex += 1;
    }
    hindex
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
/*
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
        let [y, x] = Fast_Two_Sum(Q, enow);
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
        let [y, x] = Fast_Two_Sum(hnow, Q);
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
*/

/* ****************************************************************************/
/*                                                                           */
/*  orient2d_fast()   Approximate 2D orientation test.  Nonrobust.            */
/*  orient2d_exact()   Exact 2D orientation test.  Robust.                    */
/*  orient2d_slow()   Another exact 2D orientation test.  Robust.             */
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

/// Approximate 2D orientation test. Non-robust version of [`orient2d`].
#[inline]
pub fn orient2d_fast(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> f64 {
    let acx = pa[0] - pc[0];
    let bcx = pb[0] - pc[0];
    let acy = pa[1] - pc[1];
    let bcy = pb[1] - pc[1];
    acx * bcy - acy * bcx
}

#[inline]
pub fn orient2d_exact(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> f64 {
    let [axby0, axby1] = two_product(pa[0], pb[1]);
    let [axcy0, axcy1] = two_product(pa[0], pc[1]);
    let aterms = two_two_diff(axby1, axby0, axcy1, axcy0);
    let [bxcy0, bxcy1] = two_product(pb[0], pc[1]);
    let [bxay0, bxay1] = two_product(pb[0], pa[1]);
    let bterms = two_two_diff(bxcy1, bxcy0, bxay1, bxay0);
    let [cxay0, cxay1] = two_product(pc[0], pa[1]);
    let [cxby0, cxby1] = two_product(pc[0], pb[1]);
    let cterms = two_two_diff(cxay1, cxay0, cxby1, cxby0);
    let mut v = [0.; 8];
    let vlength = fast_expansion_sum_zeroelim(&aterms, &bterms, &mut v);
    let mut w = [0.; 12];
    let wlength = fast_expansion_sum_zeroelim(&v[..vlength], &cterms, &mut w);
    w[wlength - 1]
}

#[inline]
pub fn orient2d_slow(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2]) -> f64 {
    let [acxtail, acx] = two_diff(pa[0], pc[0]);
    let [acytail, acy] = two_diff(pa[1], pc[1]);
    let [bcxtail, bcx] = two_diff(pb[0], pc[0]);
    let [bcytail, bcy] = two_diff(pb[1], pc[1]);
    let axby = two_two_product(acx, acxtail, bcy, bcytail);
    let negate = -acy;
    let negatetail = -acytail;
    let bxay = two_two_product(bcx, bcxtail, negate, negatetail);
    let mut deter = [0.; 16];
    let deterlen = fast_expansion_sum_zeroelim(&axby, &bxay, &mut deter);
    deter[deterlen - 1]
}

#[inline]
pub fn orient2dadapt(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], detsum: f64) -> f64 {
    let acx = pa[0] - pc[0];
    let bcx = pb[0] - pc[0];
    let acy = pa[1] - pc[1];
    let bcy = pb[1] - pc[1];
    let [detlefttail, detleft] = two_product(acx, bcy);
    let [detrighttail, detright] = two_product(acy, bcx);
    let b = two_two_diff(detleft, detlefttail, detright, detrighttail);
    let mut det: f64 = b.iter().sum();
    let errbound = PARAMS.ccwerrbound_b * detsum;
    if det >= errbound || -det >= errbound {
        return det;
    }
    let acxtail = two_diff_tail(pa[0], pc[0], acx);
    let bcxtail = two_diff_tail(pb[0], pc[0], bcx);
    let acytail = two_diff_tail(pa[1], pc[1], acy);
    let bcytail = two_diff_tail(pb[1], pc[1], bcy);
    if acxtail == 0.0 && acytail == 0.0 && bcxtail == 0.0 && bcytail == 0.0 {
        return det;
    }
    let errbound = PARAMS.ccwerrbound_c * detsum + PARAMS.resulterrbound * abs(det);
    det += acx * bcytail + bcy * acxtail - (acy * bcxtail + bcx * acytail);
    if det >= errbound || -det >= errbound {
        return det;
    }
    let [s0, s1] = two_product(acxtail, bcy);
    let [t0, t1] = two_product(acytail, bcx);
    let u = two_two_diff(s1, s0, t1, t0);
    let mut c1: [f64; 8] = [0.; 8];
    let c1length = fast_expansion_sum_zeroelim(&b, &u, &mut c1);
    let [s0, s1] = two_product(acx, bcytail);
    let [t0, t1] = two_product(acy, bcxtail);
    let u = two_two_diff(s1, s0, t1, t0);
    let mut c2: [f64; 12] = [0.; 12];
    let c2length = fast_expansion_sum_zeroelim(&c1[..c1length], &u, &mut c2);
    let [s0, s1] = two_product(acxtail, bcytail);
    let [t0, t1] = two_product(acytail, bcxtail);
    let u = two_two_diff(s1, s0, t1, t0);
    let mut d: [f64; 16] = [0.; 16];
    let dlength = fast_expansion_sum_zeroelim(&c2[..c2length], &u, &mut d);
    d[dlength - 1]
}

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
    let detleft = (pa[0] - pc[0]) * (pb[1] - pc[1]);
    let detright = (pa[1] - pc[1]) * (pb[0] - pc[0]);
    let det = detleft - detright;
    let detsum = if detleft > 0.0 {
        if detright <= 0.0 {
            return det;
        } else {
            detleft + detright
        }
    } else if detleft < 0.0 {
        if detright >= 0.0 {
            return det;
        } else {
            -detleft - detright
        }
    } else {
        return det;
    };
    let errbound = PARAMS.ccwerrbound_a * detsum;
    if det >= errbound || -det >= errbound {
        return det;
    }
    orient2dadapt(pa, pb, pc, detsum)
}

/* ****************************************************************************/
/*                                                                           */
/*  orient3d_fast()   Approximate 3D orientation test.  Nonrobust.            */
/*  orient3d_exact()   Exact 3D orientation test.  Robust.                    */
/*  orient3d_slow()   Another exact 3D orientation test.  Robust.             */
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

/// Approximate 3D orientation test. Non-robust version of [`orient3d`].
#[inline]
pub fn orient3d_fast(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> f64 {
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];
    let adz = pa[2] - pd[2];
    let bdz = pb[2] - pd[2];
    let cdz = pc[2] - pd[2];
    adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) + cdx * (ady * bdz - adz * bdy)
}

#[inline]
pub fn orient3d_exact(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> f64 {
    let [axby0, axby1] = two_product(pa[0], pb[1]);
    let [bxay0, bxay1] = two_product(pb[0], pa[1]);
    let ab = two_two_diff(axby1, axby0, bxay1, bxay0);
    let [bxcy0, bxcy1] = two_product(pb[0], pc[1]);
    let [cxby0, cxby1] = two_product(pc[0], pb[1]);
    let bc = two_two_diff(bxcy1, bxcy0, cxby1, cxby0);
    let [cxdy0, cxdy1] = two_product(pc[0], pd[1]);
    let [dxcy0, dxcy1] = two_product(pd[0], pc[1]);
    let cd = two_two_diff(cxdy1, cxdy0, dxcy1, dxcy0);
    let [dxay0, dxay1] = two_product(pd[0], pa[1]);
    let [axdy0, axdy1] = two_product(pa[0], pd[1]);
    let da = two_two_diff(dxay1, dxay0, axdy1, axdy0);
    let [axcy0, axcy1] = two_product(pa[0], pc[1]);
    let [cxay0, cxay1] = two_product(pc[0], pa[1]);
    let mut ac = two_two_diff(axcy1, axcy0, cxay1, cxay0);
    let [bxdy0, bxdy1] = two_product(pb[0], pd[1]);
    let [dxby0, dxby1] = two_product(pd[0], pb[1]);
    let mut bd = two_two_diff(bxdy1, bxdy0, dxby1, dxby0);

    let mut temp8 = [0.; 8];
    let mut abc = [0.; 12];
    let mut bcd = [0.; 12];
    let mut cda = [0.; 12];
    let mut dab = [0.; 12];
    let mut adet = [0.; 24];
    let mut bdet = [0.; 24];
    let mut cdet = [0.; 24];
    let mut ddet = [0.; 24];
    let mut abdet = [0.; 48];
    let mut cddet = [0.; 48];
    let mut deter = [0.; 96];

    let templen = fast_expansion_sum_zeroelim(&cd, &da, &mut temp8);
    let cdalen = fast_expansion_sum_zeroelim(&temp8[..templen], &ac, &mut cda);
    let templen = fast_expansion_sum_zeroelim(&da, &ab, &mut temp8);
    let dablen = fast_expansion_sum_zeroelim(&temp8[..templen], &bd, &mut dab);
    bd.iter_mut().for_each(|x| *x = -*x);
    ac.iter_mut().for_each(|x| *x = -*x);
    let templen = fast_expansion_sum_zeroelim(&ab, &bc, &mut temp8);
    let abclen = fast_expansion_sum_zeroelim(&temp8[..templen], &ac, &mut abc);
    let templen = fast_expansion_sum_zeroelim(&bc, &cd, &mut temp8);
    let bcdlen = fast_expansion_sum_zeroelim(&temp8[..templen], &bd, &mut bcd);
    let alen = scale_expansion_zeroelim(&bcd[..bcdlen], pa[2], &mut adet);
    let blen = scale_expansion_zeroelim(&cda[..cdalen], -pb[2], &mut bdet);
    let clen = scale_expansion_zeroelim(&dab[..dablen], pc[2], &mut cdet);
    let dlen = scale_expansion_zeroelim(&abc[..abclen], -pd[2], &mut ddet);
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let cdlen = fast_expansion_sum_zeroelim(&cdet[..clen], &ddet[..dlen], &mut cddet);
    let deterlen = fast_expansion_sum_zeroelim(&abdet[..ablen], &cddet[..cdlen], &mut deter);
    deter[deterlen - 1]
}

#[inline]
pub fn orient3d_slow(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3]) -> f64 {
    let mut temp16: [f64; 16] = [0.; 16];
    let mut temp32: [f64; 32] = [0.; 32];
    let mut temp32t: [f64; 32] = [0.; 32];
    let mut adet: [f64; 64] = [0.; 64];
    let mut bdet: [f64; 64] = [0.; 64];
    let mut cdet: [f64; 64] = [0.; 64];
    let mut abdet: [f64; 128] = [0.; 128];
    let mut deter: [f64; 192] = [0.; 192];
    let [adxtail, adx] = two_diff(pa[0], pd[0]);
    let [adytail, ady] = two_diff(pa[1], pd[1]);
    let [adztail, adz] = two_diff(pa[2], pd[2]);
    let [bdxtail, bdx] = two_diff(pb[0], pd[0]);
    let [bdytail, bdy] = two_diff(pb[1], pd[1]);
    let [bdztail, bdz] = two_diff(pb[2], pd[2]);
    let [cdxtail, cdx] = two_diff(pc[0], pd[0]);
    let [cdytail, cdy] = two_diff(pc[1], pd[1]);
    let [cdztail, cdz] = two_diff(pc[2], pd[2]);
    let axby = two_two_product(adx, adxtail, bdy, bdytail);
    let negate = -ady;
    let negatetail = -adytail;
    let bxay = two_two_product(bdx, bdxtail, negate, negatetail);
    let bxcy = two_two_product(bdx, bdxtail, cdy, cdytail);
    let negate = -bdy;
    let negatetail = -bdytail;
    let cxby = two_two_product(cdx, cdxtail, negate, negatetail);
    let cxay = two_two_product(cdx, cdxtail, ady, adytail);
    let negate = -cdy;
    let negatetail = -cdytail;
    let axcy = two_two_product(adx, adxtail, negate, negatetail);
    let temp16len = fast_expansion_sum_zeroelim(&bxcy, &cxby, &mut temp16);
    let temp32len = scale_expansion_zeroelim(&temp16[..temp16len], adz, &mut temp32);
    let temp32tlen = scale_expansion_zeroelim(&temp16[..temp16len], adztail, &mut temp32t);
    let alen = fast_expansion_sum_zeroelim(&temp32[..temp32len], &temp32t[..temp32tlen], &mut adet);
    let temp16len = fast_expansion_sum_zeroelim(&cxay, &axcy, &mut temp16);
    let temp32len = scale_expansion_zeroelim(&temp16[..temp16len], bdz, &mut temp32);
    let temp32tlen = scale_expansion_zeroelim(&temp16[..temp16len], bdztail, &mut temp32t);
    let blen = fast_expansion_sum_zeroelim(&temp32[..temp32len], &temp32t[..temp32tlen], &mut bdet);
    let temp16len = fast_expansion_sum_zeroelim(&axby, &bxay, &mut temp16);
    let temp32len = scale_expansion_zeroelim(&temp16[..temp16len], cdz, &mut temp32);
    let temp32tlen = scale_expansion_zeroelim(&temp16[..temp16len], cdztail, &mut temp32t);
    let clen = fast_expansion_sum_zeroelim(&temp32[..temp32len], &temp32t[..temp32tlen], &mut cdet);
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let deterlen = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdet[..clen], &mut deter);
    deter[deterlen - 1]
}

#[inline]
pub fn orient3dadapt(
    pa: [f64; 3],
    pb: [f64; 3],
    pc: [f64; 3],
    pd: [f64; 3],
    permanent: f64,
) -> f64 {
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];
    let adz = pa[2] - pd[2];
    let bdz = pb[2] - pd[2];
    let cdz = pc[2] - pd[2];

    let [bdxcdy0, bdxcdy1] = two_product(bdx, cdy);
    let [cdxbdy0, cdxbdy1] = two_product(cdx, bdy);
    let bc = two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0);
    let mut adet = [0.; 8];
    let alen = scale_expansion_zeroelim(&bc, adz, &mut adet);

    let [cdxady0, cdxady1] = two_product(cdx, ady);
    let [adxcdy0, adxcdy1] = two_product(adx, cdy);
    let ca = two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0);
    let mut bdet = [0.; 8];
    let blen = scale_expansion_zeroelim(&ca, bdz, &mut bdet);

    let [adxbdy0, adxbdy1] = two_product(adx, bdy);
    let [bdxady0, bdxady1] = two_product(bdx, ady);
    let ab = two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0);
    let mut cdet = [0.; 8];
    let clen = scale_expansion_zeroelim(&ab, cdz, &mut cdet);

    let mut abdet = [0.; 16];
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let mut fin1 = [0.; 192];
    let mut finlength = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdet[..clen], &mut fin1);

    let mut det: f64 = fin1[..finlength].iter().sum();
    let errbound = PARAMS.o3derrbound_b * permanent;
    if det >= errbound || -det >= errbound {
        return det;
    }

    let adxtail = two_diff_tail(pa[0], pd[0], adx);
    let bdxtail = two_diff_tail(pb[0], pd[0], bdx);
    let cdxtail = two_diff_tail(pc[0], pd[0], cdx);
    let adytail = two_diff_tail(pa[1], pd[1], ady);
    let bdytail = two_diff_tail(pb[1], pd[1], bdy);
    let cdytail = two_diff_tail(pc[1], pd[1], cdy);
    let adztail = two_diff_tail(pa[2], pd[2], adz);
    let bdztail = two_diff_tail(pb[2], pd[2], bdz);
    let cdztail = two_diff_tail(pc[2], pd[2], cdz);
    if adxtail == 0.0
        && bdxtail == 0.0
        && cdxtail == 0.0
        && adytail == 0.0
        && bdytail == 0.0
        && cdytail == 0.0
        && adztail == 0.0
        && bdztail == 0.0
        && cdztail == 0.0
    {
        return det;
    }
    let errbound = PARAMS.o3derrbound_c * permanent + PARAMS.resulterrbound * abs(det);
    det += adz * (bdx * cdytail + cdy * bdxtail - (bdy * cdxtail + cdx * bdytail))
        + adztail * (bdx * cdy - bdy * cdx)
        + (bdz * (cdx * adytail + ady * cdxtail - (cdy * adxtail + adx * cdytail))
            + bdztail * (cdx * ady - cdy * adx))
        + (cdz * (adx * bdytail + bdy * adxtail - (ady * bdxtail + bdx * adytail))
            + cdztail * (adx * bdy - ady * bdx));

    if det >= errbound || -det >= errbound {
        return det;
    }

    let at_blen;
    let at_clen;
    let at_b;
    let at_c;
    if adxtail == 0.0 {
        if adytail == 0.0 {
            at_b = [0.; 4];
            at_blen = 1;
            at_c = [0.; 4];
            at_clen = 1;
        } else {
            let negate = -adytail;
            let [at_b0, at_blarge] = two_product(negate, bdx);
            at_b = [at_b0, at_blarge, 0., 0.];
            at_blen = 2;
            let [at_c0, at_clarge] = two_product(adytail, cdx);
            at_c = [at_c0, at_clarge, 0., 0.];
            at_clen = 2;
        }
    } else if adytail == 0.0 {
        let [at_b0, at_blarge] = two_product(adxtail, bdy);
        at_b = [at_b0, at_blarge, 0., 0.];
        at_blen = 2;
        let negate = -adxtail;
        let [at_c0, at_clarge] = two_product(negate, cdy);
        at_c = [at_c0, at_clarge, 0., 0.];
        at_clen = 2;
    } else {
        let [adxt_bdy0, adxt_bdy1] = two_product(adxtail, bdy);
        let [adyt_bdx0, adyt_bdx1] = two_product(adytail, bdx);
        at_b = two_two_diff(adxt_bdy1, adxt_bdy0, adyt_bdx1, adyt_bdx0);
        at_blen = 4;
        let [adyt_cdx0, adyt_cdx1] = two_product(adytail, cdx);
        let [adxt_cdy0, adxt_cdy1] = two_product(adxtail, cdy);
        at_c = two_two_diff(adyt_cdx1, adyt_cdx0, adxt_cdy1, adxt_cdy0);
        at_clen = 4;
    }
    let bt_clen;
    let bt_alen;
    let bt_c;
    let bt_a;
    if bdxtail == 0.0 {
        if bdytail == 0.0 {
            bt_c = [0.0; 4];
            bt_clen = 1;
            bt_a = [0.0; 4];
            bt_alen = 1;
        } else {
            let negate = -bdytail;
            let [bt_c0, bt_clarge] = two_product(negate, cdx);
            bt_c = [bt_c0, bt_clarge, 0., 0.];
            bt_clen = 2;
            let [bt_a0, bt_alarge] = two_product(bdytail, adx);
            bt_a = [bt_a0, bt_alarge, 0., 0.];
            bt_alen = 2;
        }
    } else if bdytail == 0.0 {
        let [bt_c0, bt_clarge] = two_product(bdxtail, cdy);
        bt_c = [bt_c0, bt_clarge, 0., 0.];
        bt_clen = 2;
        let negate = -bdxtail;
        let [bt_a0, bt_alarge] = two_product(negate, ady);
        bt_a = [bt_a0, bt_alarge, 0., 0.];
        bt_alen = 2
    } else {
        let [bdxt_cdy0, bdxt_cdy1] = two_product(bdxtail, cdy);
        let [bdyt_cdx0, bdyt_cdx1] = two_product(bdytail, cdx);
        bt_c = two_two_diff(bdxt_cdy1, bdxt_cdy0, bdyt_cdx1, bdyt_cdx0);
        bt_clen = 4;
        let [bdyt_adx0, bdyt_adx1] = two_product(bdytail, adx);
        let [bdxt_ady0, bdxt_ady1] = two_product(bdxtail, ady);
        bt_a = two_two_diff(bdyt_adx1, bdyt_adx0, bdxt_ady1, bdxt_ady0);
        bt_alen = 4;
    }
    let ct_alen;
    let ct_blen;
    let ct_a;
    let ct_b;
    if cdxtail == 0.0 {
        if cdytail == 0.0 {
            ct_a = [0.; 4];
            ct_alen = 1;
            ct_b = [0.; 4];
            ct_blen = 1;
        } else {
            let negate = -cdytail;
            let [ct_a0, ct_alarge] = two_product(negate, adx);
            ct_a = [ct_a0, ct_alarge, 0., 0.];
            ct_alen = 2;
            let [ct_b0, ct_blarge] = two_product(cdytail, bdx);
            ct_b = [ct_b0, ct_blarge, 0., 0.];
            ct_blen = 2;
        }
    } else if cdytail == 0.0 {
        let [ct_a0, ct_alarge] = two_product(cdxtail, ady);
        ct_a = [ct_a0, ct_alarge, 0., 0.];
        ct_alen = 2;
        let negate = -cdxtail;
        let [ct_b0, ct_blarge] = two_product(negate, bdy);
        ct_b = [ct_b0, ct_blarge, 0., 0.];
        ct_blen = 2;
    } else {
        let [cdxt_ady0, cdxt_ady1] = two_product(cdxtail, ady);
        let [cdyt_adx0, cdyt_adx1] = two_product(cdytail, adx);
        ct_a = two_two_diff(cdxt_ady1, cdxt_ady0, cdyt_adx1, cdyt_adx0);
        ct_alen = 4;
        let [cdyt_bdx0, cdyt_bdx1] = two_product(cdytail, bdx);
        let [cdxt_bdy0, cdxt_bdy1] = two_product(cdxtail, bdy);
        ct_b = two_two_diff(cdyt_bdx1, cdyt_bdx0, cdxt_bdy1, cdxt_bdy0);
        ct_blen = 4;
    }

    let mut fin2 = [0.; 192];

    let mut w = [0.; 16];

    let mut bct = [0.; 8];
    let bctlen = fast_expansion_sum_zeroelim(&bt_c[..bt_clen], &ct_b[..ct_blen], &mut bct);
    let wlength = scale_expansion_zeroelim(&bct[..bctlen], adz, &mut w);
    finlength = fast_expansion_sum_zeroelim(&fin1[..finlength], &w[..wlength], &mut fin2);

    let mut cat = [0.; 8];
    let catlen = fast_expansion_sum_zeroelim(&ct_a[..ct_alen], &at_c[..at_clen], &mut cat);
    let wlength = scale_expansion_zeroelim(&cat[..catlen], bdz, &mut w);
    finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &w[..wlength], &mut fin1);

    let mut abt = [0.; 8];
    let abtlen = fast_expansion_sum_zeroelim(&at_b[..at_blen], &bt_a[..bt_alen], &mut abt);
    let wlength = scale_expansion_zeroelim(&abt[..abtlen], cdz, &mut w);
    finlength = fast_expansion_sum_zeroelim(&fin1[..finlength], &w[..wlength], &mut fin2);

    let mut v = [0.; 12];

    // TODO: replace these swaps with destructuring assignment when it is stable;
    // https://github.com/rust-lang/rfcs/pull/2909
    let (mut fin1, fin2) = if adztail != 0.0 {
        let vlength = scale_expansion_zeroelim(&bc, adztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &v[..vlength], &mut fin1);
        (fin2, fin1)
    } else {
        (fin1, fin2)
    };
    let (mut fin1, fin2) = if bdztail != 0.0 {
        let vlength = scale_expansion_zeroelim(&ca, bdztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &v[..vlength], &mut fin1);
        (fin2, fin1)
    } else {
        (fin1, fin2)
    };
    let (mut fin1, fin2) = if cdztail != 0.0 {
        let vlength = scale_expansion_zeroelim(&ab, cdztail, &mut v);
        finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &v[..vlength], &mut fin1);
        (fin2, fin1)
    } else {
        (fin1, fin2)
    };
    let (mut fin1, fin2) = if adxtail != 0.0 {
        let (mut fin1, fin2) = if bdytail != 0.0 {
            let [adxt_bdyt0, adxt_bdyt1] = two_product(adxtail, bdytail);
            let u = two_one_product(adxt_bdyt1, adxt_bdyt0, cdz);
            finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
            let (mut fin1, fin2) = (fin2, fin1);
            if cdztail != 0.0 {
                let u = two_one_product(adxt_bdyt1, adxt_bdyt0, cdztail);
                finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
                (fin2, fin1)
            } else {
                (fin1, fin2)
            }
        } else {
            (fin1, fin2)
        };
        if cdytail != 0.0 {
            let negate = -adxtail;
            let [adxt_cdyt0, adxt_cdyt1] = two_product(negate, cdytail);
            let u = two_one_product(adxt_cdyt1, adxt_cdyt0, bdz);
            finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
            let (mut fin1, fin2) = (fin2, fin1);
            if bdztail != 0.0 {
                let u = two_one_product(adxt_cdyt1, adxt_cdyt0, bdztail);
                finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
                (fin2, fin1)
            } else {
                (fin1, fin2)
            }
        } else {
            (fin1, fin2)
        }
    } else {
        (fin1, fin2)
    };
    let (mut fin1, fin2) = if bdxtail != 0.0 {
        let (mut fin1, fin2) = if cdytail != 0.0 {
            let [bdxt_cdyt0, bdxt_cdyt1] = two_product(bdxtail, cdytail);
            let u = two_one_product(bdxt_cdyt1, bdxt_cdyt0, adz);
            finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
            let (mut fin1, fin2) = (fin2, fin1);
            if adztail != 0.0 {
                let u = two_one_product(bdxt_cdyt1, bdxt_cdyt0, adztail);
                finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
                (fin2, fin1)
            } else {
                (fin1, fin2)
            }
        } else {
            (fin1, fin2)
        };
        if adytail != 0.0 {
            let negate = -bdxtail;
            let [bdxt_adyt0, bdxt_adyt1] = two_product(negate, adytail);
            let u = two_one_product(bdxt_adyt1, bdxt_adyt0, cdz);
            finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
            let (mut fin1, fin2) = (fin2, fin1);
            if cdztail != 0.0 {
                let u = two_one_product(bdxt_adyt1, bdxt_adyt0, cdztail);
                finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
                (fin2, fin1)
            } else {
                (fin1, fin2)
            }
        } else {
            (fin1, fin2)
        }
    } else {
        (fin1, fin2)
    };
    let (mut fin1, fin2) = if cdxtail != 0.0 {
        let (mut fin1, fin2) = if adytail != 0.0 {
            let [cdxt_adyt0, cdxt_adyt1] = two_product(cdxtail, adytail);
            let u = two_one_product(cdxt_adyt1, cdxt_adyt0, bdz);
            finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
            let (mut fin1, fin2) = (fin2, fin1);
            if bdztail != 0.0 {
                let u = two_one_product(cdxt_adyt1, cdxt_adyt0, bdztail);
                finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
                (fin2, fin1)
            } else {
                (fin1, fin2)
            }
        } else {
            (fin1, fin2)
        };
        if bdytail != 0.0 {
            let negate = -cdxtail;
            let [cdxt_bdyt0, cdxt_bdyt1] = two_product(negate, bdytail);
            let u = two_one_product(cdxt_bdyt1, cdxt_bdyt0, adz);
            finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
            let (mut fin1, fin2) = (fin2, fin1);
            if adztail != 0.0 {
                let u = two_one_product(cdxt_bdyt1, cdxt_bdyt0, adztail);
                finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &u, &mut fin1);
                (fin2, fin1)
            } else {
                (fin1, fin2)
            }
        } else {
            (fin1, fin2)
        }
    } else {
        (fin1, fin2)
    };
    let (mut fin1, fin2) = if adztail != 0.0 {
        let wlength = scale_expansion_zeroelim(&bct[..bctlen], adztail, &mut w);
        finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &w[..wlength], &mut fin1);
        (fin2, fin1)
    } else {
        (fin1, fin2)
    };
    let (mut fin1, fin2) = if bdztail != 0.0 {
        let wlength = scale_expansion_zeroelim(&cat[..catlen], bdztail, &mut w);
        finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &w[..wlength], &mut fin1);
        (fin2, fin1)
    } else {
        (fin1, fin2)
    };
    let fin2 = if cdztail != 0.0 {
        let wlength = scale_expansion_zeroelim(&abt[..abtlen], cdztail, &mut w);
        finlength = fast_expansion_sum_zeroelim(&fin2[..finlength], &w[..wlength], &mut fin1);
        fin1
    } else {
        fin2
    };
    fin2[finlength - 1]
}

/**
 * Adaptive exact 3D orientation test. Robust.
 *
 * Returns a positive value if the point `pd` lies below the
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
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];
    let adz = pa[2] - pd[2];
    let bdz = pb[2] - pd[2];
    let cdz = pc[2] - pd[2];
    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
    let permanent = (abs(bdxcdy) + abs(cdxbdy)) * abs(adz)
        + (abs(cdxady) + abs(adxcdy)) * abs(bdz)
        + (abs(adxbdy) + abs(bdxady)) * abs(cdz);
    let errbound = PARAMS.o3derrbound_a * permanent;
    if det > errbound || -det > errbound {
        return det;
    }
    orient3dadapt(pa, pb, pc, pd, permanent)
}

/* ****************************************************************************/
/*                                                                           */
/*  incircle_fast()   Approximate 2D incircle test.  Nonrobust.               */
/*  incircle_exact()   Exact 2D incircle test.  Robust.                       */
/*  incircle_slow()   Another exact 2D incircle test.  Robust.                */
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

/// Approximate 2D incircle test. Non-robust version of [`incircle`].
#[inline]
pub fn incircle_fast(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> f64 {
    let adx = pa[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdx = pb[0] - pd[0];
    let bdy = pb[1] - pd[1];
    let cdx = pc[0] - pd[0];
    let cdy = pc[1] - pd[1];
    let abdet = adx * bdy - bdx * ady;
    let bcdet = bdx * cdy - cdx * bdy;
    let cadet = cdx * ady - adx * cdy;
    let alift = adx * adx + ady * ady;
    let blift = bdx * bdx + bdy * bdy;
    let clift = cdx * cdx + cdy * cdy;
    alift * bcdet + blift * cadet + clift * abdet
}

#[inline]
pub fn incircle_exact(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> f64 {
    let [axby0, axby1] = two_product(pa[0], pb[1]);
    let [bxay0, bxay1] = two_product(pb[0], pa[1]);
    let ab = two_two_diff(axby1, axby0, bxay1, bxay0);
    let [bxcy0, bxcy1] = two_product(pb[0], pc[1]);
    let [cxby0, cxby1] = two_product(pc[0], pb[1]);
    let bc = two_two_diff(bxcy1, bxcy0, cxby1, cxby0);
    let [cxdy0, cxdy1] = two_product(pc[0], pd[1]);
    let [dxcy0, dxcy1] = two_product(pd[0], pc[1]);
    let cd = two_two_diff(cxdy1, cxdy0, dxcy1, dxcy0);
    let [dxay0, dxay1] = two_product(pd[0], pa[1]);
    let [axdy0, axdy1] = two_product(pa[0], pd[1]);
    let da = two_two_diff(dxay1, dxay0, axdy1, axdy0);
    let [axcy0, axcy1] = two_product(pa[0], pc[1]);
    let [cxay0, cxay1] = two_product(pc[0], pa[1]);
    let mut ac = two_two_diff(axcy1, axcy0, cxay1, cxay0);
    let [bxdy0, bxdy1] = two_product(pb[0], pd[1]);
    let [dxby0, dxby1] = two_product(pd[0], pb[1]);
    let mut bd = two_two_diff(bxdy1, bxdy0, dxby1, dxby0);
    let mut temp8 = [0.; 8];
    let templen = fast_expansion_sum_zeroelim(&cd, &da, &mut temp8);
    let mut cda = [0.; 12];
    let cdalen = fast_expansion_sum_zeroelim(&temp8[..templen], &ac, &mut cda);
    let templen = fast_expansion_sum_zeroelim(&da, &ab, &mut temp8);
    let mut dab = [0.; 12];
    let dablen = fast_expansion_sum_zeroelim(&temp8[..templen], &bd, &mut dab);
    bd.iter_mut().for_each(|x| *x = -*x);
    ac.iter_mut().for_each(|x| *x = -*x);
    let templen = fast_expansion_sum_zeroelim(&ab, &bc, &mut temp8);
    let mut abc = [0.; 12];
    let abclen = fast_expansion_sum_zeroelim(&temp8[..templen], &ac, &mut abc);
    let templen = fast_expansion_sum_zeroelim(&bc, &cd, &mut temp8);
    let mut bcd = [0.; 12];
    let bcdlen = fast_expansion_sum_zeroelim(&temp8[..templen], &bd, &mut bcd);
    let mut det24x = [0.; 24];
    let xlen = scale_expansion_zeroelim(&bcd[..bcdlen], pa[0], &mut det24x);
    let mut det48x = [0.; 48];
    let xlen = scale_expansion_zeroelim(&det24x[..xlen], pa[0], &mut det48x);
    let mut det24y = [0.; 24];
    let ylen = scale_expansion_zeroelim(&bcd[..bcdlen], pa[1], &mut det24y);
    let mut det48y = [0.; 48];
    let ylen = scale_expansion_zeroelim(&det24y[..ylen], pa[1], &mut det48y);
    let mut adet = [0.; 96];
    let alen = fast_expansion_sum_zeroelim(&det48x[..xlen], &det48y[..ylen], &mut adet);
    let xlen = scale_expansion_zeroelim(&cda[..cdalen], pb[0], &mut det24x);
    let xlen = scale_expansion_zeroelim(&det24x[..xlen], -pb[0], &mut det48x);
    let ylen = scale_expansion_zeroelim(&cda[..cdalen], pb[1], &mut det24y);
    let ylen = scale_expansion_zeroelim(&det24y[..ylen], -pb[1], &mut det48y);
    let mut bdet = [0.; 96];
    let blen = fast_expansion_sum_zeroelim(&det48x[..xlen], &det48y[..ylen], &mut bdet);
    let xlen = scale_expansion_zeroelim(&dab[..dablen], pc[0], &mut det24x);
    let xlen = scale_expansion_zeroelim(&det24x[..xlen], pc[0], &mut det48x);
    let ylen = scale_expansion_zeroelim(&dab[..dablen], pc[1], &mut det24y);
    let ylen = scale_expansion_zeroelim(&det24y[..ylen], pc[1], &mut det48y);
    let mut cdet = [0.; 96];
    let clen = fast_expansion_sum_zeroelim(&det48x[..xlen], &det48y[..ylen], &mut cdet);
    let xlen = scale_expansion_zeroelim(&abc[..abclen], pd[0], &mut det24x);
    let xlen = scale_expansion_zeroelim(&det24x[..xlen], -pd[0], &mut det48x);
    let ylen = scale_expansion_zeroelim(&abc[..abclen], pd[1], &mut det24y);
    let ylen = scale_expansion_zeroelim(&det24y[..ylen], -pd[1], &mut det48y);
    let mut ddet = [0.; 96];
    let dlen = fast_expansion_sum_zeroelim(&det48x[..xlen], &det48y[..ylen], &mut ddet);
    let mut abdet = [0.; 192];
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let mut cddet = [0.; 192];
    let cdlen = fast_expansion_sum_zeroelim(&cdet[..clen], &ddet[..dlen], &mut cddet);
    let mut deter = [0.; 384];
    let deterlen = fast_expansion_sum_zeroelim(&abdet[..ablen], &cddet[..cdlen], &mut deter);
    deter[deterlen - 1]
}

#[inline]
pub fn incircle_slow(pa: [f64; 2], pb: [f64; 2], pc: [f64; 2], pd: [f64; 2]) -> f64 {
    let [adxtail, adx] = two_diff(pa[0], pd[0]);
    let [adytail, ady] = two_diff(pa[1], pd[1]);
    let [bdxtail, bdx] = two_diff(pb[0], pd[0]);
    let [bdytail, bdy] = two_diff(pb[1], pd[1]);
    let [cdxtail, cdx] = two_diff(pc[0], pd[0]);
    let [cdytail, cdy] = two_diff(pc[1], pd[1]);
    let axby = two_two_product(adx, adxtail, bdy, bdytail);
    let negate = -ady;
    let negatetail = -adytail;
    let bxay = two_two_product(bdx, bdxtail, negate, negatetail);
    let bxcy = two_two_product(bdx, bdxtail, cdy, cdytail);
    let negate = -bdy;
    let negatetail = -bdytail;
    let cxby = two_two_product(cdx, cdxtail, negate, negatetail);
    let cxay = two_two_product(cdx, cdxtail, ady, adytail);
    let negate = -cdy;
    let negatetail = -cdytail;
    let axcy = two_two_product(adx, adxtail, negate, negatetail);
    let mut temp16 = [0.; 16];
    let temp16len = fast_expansion_sum_zeroelim(&bxcy, &cxby, &mut temp16);
    let mut detx = [0.; 32];
    let xlen = scale_expansion_zeroelim(&temp16[..temp16len], adx, &mut detx);
    let mut detxx = [0.; 64];
    let xxlen = scale_expansion_zeroelim(&detx[..xlen], adx, &mut detxx);
    let mut detxt = [0.; 32];
    let xtlen = scale_expansion_zeroelim(&temp16[..temp16len], adxtail, &mut detxt);
    let mut detxxt = [0.; 64];
    let xxtlen = scale_expansion_zeroelim(&detxt[..xtlen], adx, &mut detxxt);
    detxxt[..xxtlen].iter_mut().for_each(|x| *x *= 2.0);
    let mut detxtxt = [0.; 64];
    let xtxtlen = scale_expansion_zeroelim(&detxt[..xtlen], adxtail, &mut detxtxt);
    let mut x1 = [0.; 128];
    let x1len = fast_expansion_sum_zeroelim(&detxx[..xxlen], &detxxt[..xxtlen], &mut x1);
    let mut x2 = [0.; 192];
    let x2len = fast_expansion_sum_zeroelim(&x1[..x1len], &detxtxt[..xtxtlen], &mut x2);
    let mut dety = [0.; 32];
    let ylen = scale_expansion_zeroelim(&temp16[..temp16len], ady, &mut dety);
    let mut detyy = [0.; 64];
    let yylen = scale_expansion_zeroelim(&dety[..ylen], ady, &mut detyy);
    let mut detyt = [0.; 32];
    let ytlen = scale_expansion_zeroelim(&temp16[..temp16len], adytail, &mut detyt);
    let mut detyyt = [0.; 64];
    let yytlen = scale_expansion_zeroelim(&detyt[..ytlen], ady, &mut detyyt);
    detyyt[..yytlen].iter_mut().for_each(|x| *x *= 2.0);
    let mut detytyt = [0.; 64];
    let ytytlen = scale_expansion_zeroelim(&detyt[..ytlen], adytail, &mut detytyt);
    let mut y1 = [0.; 128];
    let y1len = fast_expansion_sum_zeroelim(&detyy[..yylen], &detyyt[..yytlen], &mut y1);
    let mut y2 = [0.; 192];
    let y2len = fast_expansion_sum_zeroelim(&y1[..y1len], &detytyt[..ytytlen], &mut y2);
    let mut adet = [0.; 384];
    let alen = fast_expansion_sum_zeroelim(&x2[..x2len], &y2[..y2len], &mut adet);
    let temp16len = fast_expansion_sum_zeroelim(&cxay, &axcy, &mut temp16);
    let xlen = scale_expansion_zeroelim(&temp16[..temp16len], bdx, &mut detx);
    let xxlen = scale_expansion_zeroelim(&detx[..xlen], bdx, &mut detxx);
    let xtlen = scale_expansion_zeroelim(&temp16[..temp16len], bdxtail, &mut detxt);
    let xxtlen = scale_expansion_zeroelim(&detxt[..xtlen], bdx, &mut detxxt);
    detxxt[..xxtlen].iter_mut().for_each(|x| *x *= 2.0);
    let xtxtlen = scale_expansion_zeroelim(&detxt[..xtlen], bdxtail, &mut detxtxt);
    let x1len = fast_expansion_sum_zeroelim(&detxx[..xxlen], &detxxt[..xxtlen], &mut x1);
    let x2len = fast_expansion_sum_zeroelim(&x1[..x1len], &detxtxt[..xtxtlen], &mut x2);
    let ylen = scale_expansion_zeroelim(&temp16[..temp16len], bdy, &mut dety);
    let yylen = scale_expansion_zeroelim(&dety[..ylen], bdy, &mut detyy);
    let ytlen = scale_expansion_zeroelim(&temp16[..temp16len], bdytail, &mut detyt);
    let yytlen = scale_expansion_zeroelim(&detyt[..ytlen], bdy, &mut detyyt);
    detyyt[..yytlen].iter_mut().for_each(|x| *x *= 2.0);
    let ytytlen = scale_expansion_zeroelim(&detyt[..ytlen], bdytail, &mut detytyt);
    let y1len = fast_expansion_sum_zeroelim(&detyy[..yylen], &detyyt[..yytlen], &mut y1);
    let y2len = fast_expansion_sum_zeroelim(&y1[..y1len], &detytyt[..ytytlen], &mut y2);
    let mut bdet = [0.; 384];
    let blen = fast_expansion_sum_zeroelim(&x2[..x2len], &y2[..y2len], &mut bdet);
    let temp16len = fast_expansion_sum_zeroelim(&axby, &bxay, &mut temp16);
    let xlen = scale_expansion_zeroelim(&temp16[..temp16len], cdx, &mut detx);
    let xxlen = scale_expansion_zeroelim(&detx[..xlen], cdx, &mut detxx);
    let xtlen = scale_expansion_zeroelim(&temp16[..temp16len], cdxtail, &mut detxt);
    let xxtlen = scale_expansion_zeroelim(&detxt[..xtlen], cdx, &mut detxxt);
    detxxt[..xxtlen].iter_mut().for_each(|x| *x *= 2.0);
    let xtxtlen = scale_expansion_zeroelim(&detxt[..xtlen], cdxtail, &mut detxtxt);
    let x1len = fast_expansion_sum_zeroelim(&detxx[..xxlen], &detxxt[..xxtlen], &mut x1);
    let x2len = fast_expansion_sum_zeroelim(&x1[..x1len], &detxtxt[..xtxtlen], &mut x2);
    let ylen = scale_expansion_zeroelim(&temp16[..temp16len], cdy, &mut dety);
    let yylen = scale_expansion_zeroelim(&dety[..ylen], cdy, &mut detyy);
    let ytlen = scale_expansion_zeroelim(&temp16[..temp16len], cdytail, &mut detyt);
    let yytlen = scale_expansion_zeroelim(&detyt[..ytlen], cdy, &mut detyyt);
    detyyt[..yytlen].iter_mut().for_each(|x| *x *= 2.0);
    let ytytlen = scale_expansion_zeroelim(&detyt[..ytlen], cdytail, &mut detytyt);
    let y1len = fast_expansion_sum_zeroelim(&detyy[..yylen], &detyyt[..yytlen], &mut y1);
    let y2len = fast_expansion_sum_zeroelim(&y1[..y1len], &detytyt[..ytytlen], &mut y2);
    let mut cdet = [0.; 384];
    let clen = fast_expansion_sum_zeroelim(&x2[..x2len], &y2[..y2len], &mut cdet);
    let mut abdet = [0.; 768];
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let mut deter = [0.; 1152];
    let deterlen = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdet[..clen], &mut deter);
    deter[deterlen - 1]
}

#[inline]
pub fn incircleadapt(
    pa: [f64; 2],
    pb: [f64; 2],
    pc: [f64; 2],
    pd: [f64; 2],
    permanent: f64,
) -> f64 {
    let mut axbc = [0.; 8];
    let mut axxbc = [0.; 16];
    let mut aybc = [0.; 8];
    let mut ayybc = [0.; 16];
    let mut adet = [0.; 32];
    let mut bxca = [0.; 8];
    let mut bxxca = [0.; 16];
    let mut byca = [0.; 8];
    let mut byyca = [0.; 16];
    let mut bdet = [0.; 32];
    let mut cxab = [0.; 8];
    let mut cxxab = [0.; 16];
    let mut cyab = [0.; 8];
    let mut cyyab = [0.; 16];
    let mut cdet = [0.; 32];
    let mut abdet = [0.; 64];
    let mut fin1 = [0.; 1152];
    let mut fin2 = [0.; 1152];
    let mut finnow: &mut [f64];
    let mut finother: &mut [f64];
    let mut temp8 = [0.; 8];
    let mut temp16a = [0.; 16];
    let mut temp16b = [0.; 16];
    let mut temp16c = [0.; 16];
    let mut temp32a = [0.; 32];
    let mut temp32b = [0.; 32];
    let mut temp48 = [0.; 48];
    let mut temp64 = [0.; 64];
    let mut axtbb = [0.; 8];
    let mut axtcc = [0.; 8];
    let mut aytbb = [0.; 8];
    let mut aytcc = [0.; 8];
    let mut bxtaa = [0.; 8];
    let mut bxtcc = [0.; 8];
    let mut bytaa = [0.; 8];
    let mut bytcc = [0.; 8];
    let mut cxtaa = [0.; 8];
    let mut cxtbb = [0.; 8];
    let mut cytaa = [0.; 8];
    let mut cytbb = [0.; 8];
    let mut axtbct = [0.; 16];
    let mut aytbct = [0.; 16];
    let mut bxtcat = [0.; 16];
    let mut bytcat = [0.; 16];
    let mut cxtabt = [0.; 16];
    let mut cytabt = [0.; 16];
    let mut axtbctt = [0.; 8];
    let mut aytbctt = [0.; 8];
    let mut bxtcatt = [0.; 8];
    let mut bytcatt = [0.; 8];
    let mut cxtabtt = [0.; 8];
    let mut cytabtt = [0.; 8];
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];
    let [bdxcdy0, bdxcdy1] = two_product(bdx, cdy);
    let [cdxbdy0, cdxbdy1] = two_product(cdx, bdy);
    let bc = two_two_diff(bdxcdy1, bdxcdy0, cdxbdy1, cdxbdy0);
    let axbclen = scale_expansion_zeroelim(&bc, adx, &mut axbc);
    let axxbclen = scale_expansion_zeroelim(&axbc[..axbclen], adx, &mut axxbc);
    let aybclen = scale_expansion_zeroelim(&bc, ady, &mut aybc);
    let ayybclen = scale_expansion_zeroelim(&aybc[..aybclen], ady, &mut ayybc);
    let alen = fast_expansion_sum_zeroelim(&axxbc[..axxbclen], &ayybc[..ayybclen], &mut adet);
    let [cdxady0, cdxady1] = two_product(cdx, ady);
    let [adxcdy0, adxcdy1] = two_product(adx, cdy);
    let ca = two_two_diff(cdxady1, cdxady0, adxcdy1, adxcdy0);
    let bxcalen = scale_expansion_zeroelim(&ca, bdx, &mut bxca);
    let bxxcalen = scale_expansion_zeroelim(&bxca[..bxcalen], bdx, &mut bxxca);
    let bycalen = scale_expansion_zeroelim(&ca, bdy, &mut byca);
    let byycalen = scale_expansion_zeroelim(&byca[..bycalen], bdy, &mut byyca);
    let blen = fast_expansion_sum_zeroelim(&bxxca[..bxxcalen], &byyca[..byycalen], &mut bdet);
    let [adxbdy0, adxbdy1] = two_product(adx, bdy);
    let [bdxady0, bdxady1] = two_product(bdx, ady);
    let ab = two_two_diff(adxbdy1, adxbdy0, bdxady1, bdxady0);
    let cxablen = scale_expansion_zeroelim(&ab, cdx, &mut cxab);
    let cxxablen = scale_expansion_zeroelim(&cxab[..cxablen], cdx, &mut cxxab);
    let cyablen = scale_expansion_zeroelim(&ab, cdy, &mut cyab);
    let cyyablen = scale_expansion_zeroelim(&cyab[..cyablen], cdy, &mut cyyab);
    let clen = fast_expansion_sum_zeroelim(&cxxab[..cxxablen], &cyyab[..cyyablen], &mut cdet);
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let mut finlength = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdet[..clen], &mut fin1);
    let mut det: f64 = fin1[..finlength].iter().sum();
    let errbound = PARAMS.iccerrbound_b * permanent;
    if det >= errbound || -det >= errbound {
        return det;
    }
    let adxtail = two_diff_tail(pa[0], pd[0], adx);
    let adytail = two_diff_tail(pa[1], pd[1], ady);
    let bdxtail = two_diff_tail(pb[0], pd[0], bdx);
    let bdytail = two_diff_tail(pb[1], pd[1], bdy);
    let cdxtail = two_diff_tail(pc[0], pd[0], cdx);
    let cdytail = two_diff_tail(pc[1], pd[1], cdy);
    if adxtail == 0.0
        && bdxtail == 0.0
        && cdxtail == 0.0
        && adytail == 0.0
        && bdytail == 0.0
        && cdytail == 0.0
    {
        return det;
    }
    let errbound = PARAMS.iccerrbound_c * permanent + PARAMS.resulterrbound * abs(det);
    det += (adx * adx + ady * ady)
        * (bdx * cdytail + cdy * bdxtail - (bdy * cdxtail + cdx * bdytail))
        + 2.0 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx)
        + ((bdx * bdx + bdy * bdy)
            * (cdx * adytail + ady * cdxtail - (cdy * adxtail + adx * cdytail))
            + 2.0 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx))
        + ((cdx * cdx + cdy * cdy)
            * (adx * bdytail + bdy * adxtail - (ady * bdxtail + bdx * adytail))
            + 2.0 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx));
    if det >= errbound || -det >= errbound {
        return det;
    }
    finnow = &mut fin1;
    finother = &mut fin2;
    let aa = if bdxtail != 0.0 || bdytail != 0.0 || cdxtail != 0.0 || cdytail != 0.0 {
        let [adxadx0, adxadx1] = square(adx);
        let [adyady0, adyady1] = square(ady);
        two_two_sum(adxadx1, adxadx0, adyady1, adyady0)
    } else {
        [0.; 4]
    };
    let bb = if cdxtail != 0.0 || cdytail != 0.0 || adxtail != 0.0 || adytail != 0.0 {
        let [bdxbdx0, bdxbdx1] = square(bdx);
        let [bdybdy0, bdybdy1] = square(bdy);
        two_two_sum(bdxbdx1, bdxbdx0, bdybdy1, bdybdy0)
    } else {
        [0.; 4]
    };
    let cc = if adxtail != 0.0 || adytail != 0.0 || bdxtail != 0.0 || bdytail != 0.0 {
        let [cdxcdx0, cdxcdx1] = square(cdx);
        let [cdycdy0, cdycdy1] = square(cdy);
        two_two_sum(cdxcdx1, cdxcdx0, cdycdy1, cdycdy0)
    } else {
        [0.; 4]
    };
    let mut axtbclen = 8;
    let mut axtbc = [0.; 8];
    if adxtail != 0.0 {
        axtbclen = scale_expansion_zeroelim(&bc, adxtail, &mut axtbc);
        let temp16alen = scale_expansion_zeroelim(&axtbc[..axtbclen], 2.0 * adx, &mut temp16a);
        let axtcclen = scale_expansion_zeroelim(&cc, adxtail, &mut axtcc);
        let temp16blen = scale_expansion_zeroelim(&axtcc[..axtcclen], bdy, &mut temp16b);
        let axtbblen = scale_expansion_zeroelim(&bb, adxtail, &mut axtbb);
        let temp16clen = scale_expansion_zeroelim(&axtbb[..axtbblen], -cdy, &mut temp16c);
        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
        core::mem::swap(&mut finnow, &mut finother);
    }
    let mut aytbclen = 8;
    let mut aytbc = [0.; 8];
    if adytail != 0.0 {
        aytbclen = scale_expansion_zeroelim(&bc, adytail, &mut aytbc);
        let temp16alen = scale_expansion_zeroelim(&aytbc[..aytbclen], 2.0 * ady, &mut temp16a);
        let aytbblen = scale_expansion_zeroelim(&bb, adytail, &mut aytbb);
        let temp16blen = scale_expansion_zeroelim(&aytbb[..aytbblen], cdx, &mut temp16b);
        let aytcclen = scale_expansion_zeroelim(&cc, adytail, &mut aytcc);
        let temp16clen = scale_expansion_zeroelim(&aytcc[..aytcclen], -bdx, &mut temp16c);
        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
        core::mem::swap(&mut finnow, &mut finother);
    }
    let mut bxtcalen = 8;
    let mut bxtca = [0.; 8];
    if bdxtail != 0.0 {
        bxtcalen = scale_expansion_zeroelim(&ca, bdxtail, &mut bxtca);
        let temp16alen = scale_expansion_zeroelim(&bxtca[..bxtcalen], 2.0 * bdx, &mut temp16a);
        let bxtaalen = scale_expansion_zeroelim(&aa, bdxtail, &mut bxtaa);
        let temp16blen = scale_expansion_zeroelim(&bxtaa[..bxtaalen], cdy, &mut temp16b);
        let bxtcclen = scale_expansion_zeroelim(&cc, bdxtail, &mut bxtcc);
        let temp16clen = scale_expansion_zeroelim(&bxtcc[..bxtcclen], -ady, &mut temp16c);
        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
        core::mem::swap(&mut finnow, &mut finother);
    }
    let mut bytcalen = 8;
    let mut bytca = [0.; 8];
    if bdytail != 0.0 {
        bytcalen = scale_expansion_zeroelim(&ca, bdytail, &mut bytca);
        let temp16alen = scale_expansion_zeroelim(&bytca[..bytcalen], 2.0 * bdy, &mut temp16a);
        let bytcclen = scale_expansion_zeroelim(&cc, bdytail, &mut bytcc);
        let temp16blen = scale_expansion_zeroelim(&bytcc[..bytcclen], adx, &mut temp16b);
        let bytaalen = scale_expansion_zeroelim(&aa, bdytail, &mut bytaa);
        let temp16clen = scale_expansion_zeroelim(&bytaa[..bytaalen], -cdx, &mut temp16c);
        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
        core::mem::swap(&mut finnow, &mut finother);
    }
    let cxtablen = 8;
    let mut cxtab = [0.; 8];
    if cdxtail != 0.0 {
        let cxtablen = scale_expansion_zeroelim(&ab, cdxtail, &mut cxtab);
        let temp16alen = scale_expansion_zeroelim(&cxtab[..cxtablen], 2.0 * cdx, &mut temp16a);
        let cxtbblen = scale_expansion_zeroelim(&bb, cdxtail, &mut cxtbb);
        let temp16blen = scale_expansion_zeroelim(&cxtbb[..cxtbblen], ady, &mut temp16b);
        let cxtaalen = scale_expansion_zeroelim(&aa, cdxtail, &mut cxtaa);
        let temp16clen = scale_expansion_zeroelim(&cxtaa[..cxtaalen], -bdy, &mut temp16c);
        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
        core::mem::swap(&mut finnow, &mut finother);
    }
    let mut cytablen = 8;
    let mut cytab = [0.; 8];
    if cdytail != 0.0 {
        cytablen = scale_expansion_zeroelim(&ab, cdytail, &mut cytab);
        let temp16alen = scale_expansion_zeroelim(&cytab[..cytablen], 2.0 * cdy, &mut temp16a);
        let cytaalen = scale_expansion_zeroelim(&aa, cdytail, &mut cytaa);
        let temp16blen = scale_expansion_zeroelim(&cytaa[..cytaalen], bdx, &mut temp16b);
        let cytbblen = scale_expansion_zeroelim(&bb, cdytail, &mut cytbb);
        let temp16clen = scale_expansion_zeroelim(&cytbb[..cytbblen], -adx, &mut temp16c);
        let temp32alen = fast_expansion_sum_zeroelim(
            &temp16a[..temp16alen],
            &temp16b[..temp16blen],
            &mut temp32a,
        );
        let temp48len = fast_expansion_sum_zeroelim(
            &temp16c[..temp16clen],
            &temp32a[..temp32alen],
            &mut temp48,
        );
        finlength =
            fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
        core::mem::swap(&mut finnow, &mut finother);
    }
    if adxtail != 0.0 || adytail != 0.0 {
        let mut bctlen = 1;
        let mut bcttlen = 1;
        let mut bct: [f64; 8] = [0.; 8];
        let mut bctt: [f64; 4] = [0.; 4];
        if bdxtail != 0.0 || bdytail != 0.0 || cdxtail != 0.0 || cdytail != 0.0 {
            let [ti0, ti1] = two_product(bdxtail, cdy);
            let [tj0, tj1] = two_product(bdx, cdytail);
            let u = two_two_sum(ti1, ti0, tj1, tj0);
            let negate = -bdy;
            let [ti0, ti1] = two_product(cdxtail, negate);
            let negate = -bdytail;
            let [tj0, tj1] = two_product(cdx, negate);
            let v = two_two_sum(ti1, ti0, tj1, tj0);
            bctlen = fast_expansion_sum_zeroelim(&u, &v, &mut bct);
            let [ti0, ti1] = two_product(bdxtail, cdytail);
            let [tj0, tj1] = two_product(cdxtail, bdytail);
            bctt = two_two_diff(ti1, ti0, tj1, tj0);
            bcttlen = 4;
        }
        if adxtail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&axtbc[..axtbclen], adxtail, &mut temp16a);
            let axtbctlen = scale_expansion_zeroelim(&bct[..bctlen], adxtail, &mut axtbct);
            let temp32alen =
                scale_expansion_zeroelim(&axtbct[..axtbctlen], 2.0 * adx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
            core::mem::swap(&mut finnow, &mut finother);
            if bdytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&cc, adxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], bdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &finnow[..finlength],
                    &temp16a[..temp16alen],
                    finother,
                );
                core::mem::swap(&mut finnow, &mut finother);
            }
            if cdytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&bb, -adxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], cdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &finnow[..finlength],
                    &temp16a[..temp16alen],
                    finother,
                );
                core::mem::swap(&mut finnow, &mut finother);
            }
            let temp32alen = scale_expansion_zeroelim(&axtbct[..axtbctlen], adxtail, &mut temp32a);
            let axtbcttlen = scale_expansion_zeroelim(&bctt[..bcttlen], adxtail, &mut axtbctt);
            let temp16alen =
                scale_expansion_zeroelim(&axtbctt[..axtbcttlen], 2.0 * adx, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&axtbctt[..axtbcttlen], adxtail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp64[..temp64len], finother);
            core::mem::swap(&mut finnow, &mut finother);
        }
        if adytail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&aytbc[..aytbclen], adytail, &mut temp16a);
            let aytbctlen = scale_expansion_zeroelim(&bct[..bctlen], adytail, &mut aytbct);
            let temp32alen =
                scale_expansion_zeroelim(&aytbct[..aytbctlen], 2.0 * ady, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
            core::mem::swap(&mut finnow, &mut finother);
            let temp32alen = scale_expansion_zeroelim(&aytbct[..aytbctlen], adytail, &mut temp32a);
            let aytbcttlen = scale_expansion_zeroelim(&bctt[..bcttlen], adytail, &mut aytbctt);
            let temp16alen =
                scale_expansion_zeroelim(&aytbctt[..aytbcttlen], 2.0 * ady, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&aytbctt[..aytbcttlen], adytail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp64[..temp64len], finother);
            core::mem::swap(&mut finnow, &mut finother);
        }
    }
    if bdxtail != 0.0 || bdytail != 0.0 {
        let mut catlen = 1;
        let mut cattlen = 1;
        let mut cat: [f64; 8] = [0.; 8];
        let mut catt: [f64; 4] = [0.; 4];
        if cdxtail != 0.0 || cdytail != 0.0 || adxtail != 0.0 || adytail != 0.0 {
            let [ti0, ti1] = two_product(cdxtail, ady);
            let [tj0, tj1] = two_product(cdx, adytail);
            let u = two_two_sum(ti1, ti0, tj1, tj0);
            let negate = -cdy;
            let [ti0, ti1] = two_product(adxtail, negate);
            let negate = -cdytail;
            let [tj0, tj1] = two_product(adx, negate);
            let v = two_two_sum(ti1, ti0, tj1, tj0);
            catlen = fast_expansion_sum_zeroelim(&u, &v, &mut cat);
            let [ti0, ti1] = two_product(cdxtail, adytail);
            let [tj0, tj1] = two_product(adxtail, cdytail);
            catt = two_two_diff(ti1, ti0, tj1, tj0);
            cattlen = 4;
        }
        if bdxtail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&bxtca[..bxtcalen], bdxtail, &mut temp16a);
            let bxtcatlen = scale_expansion_zeroelim(&cat[..catlen], bdxtail, &mut bxtcat);
            let temp32alen =
                scale_expansion_zeroelim(&bxtcat[..bxtcatlen], 2.0 * bdx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
            core::mem::swap(&mut finnow, &mut finother);
            if cdytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&aa, bdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], cdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &finnow[..finlength],
                    &temp16a[..temp16alen],
                    finother,
                );
                core::mem::swap(&mut finnow, &mut finother);
            }
            if adytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&cc, -bdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], adytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &finnow[..finlength],
                    &temp16a[..temp16alen],
                    finother,
                );
                core::mem::swap(&mut finnow, &mut finother);
            }
            let temp32alen = scale_expansion_zeroelim(&bxtcat[..bxtcatlen], bdxtail, &mut temp32a);
            let bxtcattlen = scale_expansion_zeroelim(&catt[..cattlen], bdxtail, &mut bxtcatt);
            let temp16alen =
                scale_expansion_zeroelim(&bxtcatt[..bxtcattlen], 2.0 * bdx, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&bxtcatt[..bxtcattlen], bdxtail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp64[..temp64len], finother);
            core::mem::swap(&mut finnow, &mut finother);
        }
        if bdytail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&bytca[..bytcalen], bdytail, &mut temp16a);
            let bytcatlen = scale_expansion_zeroelim(&cat[..catlen], bdytail, &mut bytcat);
            let temp32alen =
                scale_expansion_zeroelim(&bytcat[..bytcatlen], 2.0 * bdy, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
            core::mem::swap(&mut finnow, &mut finother);
            let temp32alen = scale_expansion_zeroelim(&bytcat[..bytcatlen], bdytail, &mut temp32a);
            let bytcattlen = scale_expansion_zeroelim(&catt[..cattlen], bdytail, &mut bytcatt);
            let temp16alen =
                scale_expansion_zeroelim(&bytcatt[..bytcattlen], 2.0 * bdy, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&bytcatt[..bytcattlen], bdytail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp64[..temp64len], finother);
            core::mem::swap(&mut finnow, &mut finother);
        }
    }
    if cdxtail != 0.0 || cdytail != 0.0 {
        let mut abtlen = 1;
        let mut abttlen = 1;
        let mut abt: [f64; 8] = [0.; 8];
        let mut abtt: [f64; 4] = [0.; 4];
        if adxtail != 0.0 || adytail != 0.0 || bdxtail != 0.0 || bdytail != 0.0 {
            let [ti0, ti1] = two_product(adxtail, bdy);
            let [tj0, tj1] = two_product(adx, bdytail);
            let u = two_two_sum(ti1, ti0, tj1, tj0);
            let negate = -ady;
            let [ti0, ti1] = two_product(bdxtail, negate);
            let negate = -adytail;
            let [tj0, tj1] = two_product(bdx, negate);
            let v = two_two_sum(ti1, ti0, tj1, tj0);
            abtlen = fast_expansion_sum_zeroelim(&u, &v, &mut abt);
            let [ti0, ti1] = two_product(adxtail, bdytail);
            let [tj0, tj1] = two_product(bdxtail, adytail);
            abtt = two_two_diff(ti1, ti0, tj1, tj0);
            abttlen = 4;
        }
        if cdxtail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&cxtab[..cxtablen], cdxtail, &mut temp16a);
            let cxtabtlen = scale_expansion_zeroelim(&abt[..abtlen], cdxtail, &mut cxtabt);
            let temp32alen =
                scale_expansion_zeroelim(&cxtabt[..cxtabtlen], 2.0 * cdx, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
            core::mem::swap(&mut finnow, &mut finother);
            if adytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&bb, cdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], adytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &finnow[..finlength],
                    &temp16a[..temp16alen],
                    finother,
                );
                core::mem::swap(&mut finnow, &mut finother);
            }
            if bdytail != 0.0 {
                let temp8len = scale_expansion_zeroelim(&aa, -cdxtail, &mut temp8);
                let temp16alen =
                    scale_expansion_zeroelim(&temp8[..temp8len], bdytail, &mut temp16a);
                finlength = fast_expansion_sum_zeroelim(
                    &finnow[..finlength],
                    &temp16a[..temp16alen],
                    finother,
                );
                core::mem::swap(&mut finnow, &mut finother);
            }
            let temp32alen = scale_expansion_zeroelim(&cxtabt[..cxtabtlen], cdxtail, &mut temp32a);
            let cxtabttlen = scale_expansion_zeroelim(&abtt[..abttlen], cdxtail, &mut cxtabtt);
            let temp16alen =
                scale_expansion_zeroelim(&cxtabtt[..cxtabttlen], 2.0 * cdx, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&cxtabtt[..cxtabttlen], cdxtail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp64[..temp64len], finother);
            core::mem::swap(&mut finnow, &mut finother);
        }
        if cdytail != 0.0 {
            let temp16alen = scale_expansion_zeroelim(&cytab[..cytablen], cdytail, &mut temp16a);
            let cytabtlen = scale_expansion_zeroelim(&abt[..abtlen], cdytail, &mut cytabt);
            let temp32alen =
                scale_expansion_zeroelim(&cytabt[..cytabtlen], 2.0 * cdy, &mut temp32a);
            let temp48len = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp32a[..temp32alen],
                &mut temp48,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp48[..temp48len], finother);
            core::mem::swap(&mut finnow, &mut finother);
            let temp32alen = scale_expansion_zeroelim(&cytabt[..cytabtlen], cdytail, &mut temp32a);
            let cytabttlen = scale_expansion_zeroelim(&abtt[..abttlen], cdytail, &mut cytabtt);
            let temp16alen =
                scale_expansion_zeroelim(&cytabtt[..cytabttlen], 2.0 * cdy, &mut temp16a);
            let temp16blen =
                scale_expansion_zeroelim(&cytabtt[..cytabttlen], cdytail, &mut temp16b);
            let temp32blen = fast_expansion_sum_zeroelim(
                &temp16a[..temp16alen],
                &temp16b[..temp16blen],
                &mut temp32b,
            );
            let temp64len = fast_expansion_sum_zeroelim(
                &temp32a[..temp32alen],
                &temp32b[..temp32blen],
                &mut temp64,
            );
            finlength =
                fast_expansion_sum_zeroelim(&finnow[..finlength], &temp64[..temp64len], finother);
            core::mem::swap(&mut finnow, &mut finother);
        }
    }
    finnow[finlength - 1]
}

/**
 * Adaptive exact 2D incircle test. Robust.
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
    let adx = pa[0] - pd[0];
    let bdx = pb[0] - pd[0];
    let cdx = pc[0] - pd[0];
    let ady = pa[1] - pd[1];
    let bdy = pb[1] - pd[1];
    let cdy = pc[1] - pd[1];
    let bdxcdy = bdx * cdy;
    let cdxbdy = cdx * bdy;
    let alift = adx * adx + ady * ady;
    let cdxady = cdx * ady;
    let adxcdy = adx * cdy;
    let blift = bdx * bdx + bdy * bdy;
    let adxbdy = adx * bdy;
    let bdxady = bdx * ady;
    let clift = cdx * cdx + cdy * cdy;
    let det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
    let permanent = (abs(bdxcdy) + abs(cdxbdy)) * alift
        + (abs(cdxady) + abs(adxcdy)) * blift
        + (abs(adxbdy) + abs(bdxady)) * clift;
    let errbound = PARAMS.iccerrbound_a * permanent;
    if det > errbound || -det > errbound {
        return det;
    }
    incircleadapt(pa, pb, pc, pd, permanent)
}

/* ****************************************************************************/
/*                                                                           */
/*  insphere_fast()   Approximate 3D insphere test.  Nonrobust.               */
/*  insphere_exact()   Exact 3D insphere test.  Robust.                       */
/*  insphere_slow()   Another exact 3D insphere test.  Robust.                */
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

/// Approximate 3D insphere test. Non-robust version of [`insphere`].
#[inline]
pub fn insphere_fast(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3], pe: [f64; 3]) -> f64 {
    let aex = pa[0] - pe[0];
    let bex = pb[0] - pe[0];
    let cex = pc[0] - pe[0];
    let dex = pd[0] - pe[0];
    let aey = pa[1] - pe[1];
    let bey = pb[1] - pe[1];
    let cey = pc[1] - pe[1];
    let dey = pd[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bez = pb[2] - pe[2];
    let cez = pc[2] - pe[2];
    let dez = pd[2] - pe[2];
    let ab = aex * bey - bex * aey;
    let bc = bex * cey - cex * bey;
    let cd = cex * dey - dex * cey;
    let da = dex * aey - aex * dey;
    let ac = aex * cey - cex * aey;
    let bd = bex * dey - dex * bey;
    let abc = aez * bc - bez * ac + cez * ab;
    let bcd = bez * cd - cez * bd + dez * bc;
    let cda = cez * da + dez * ac + aez * cd;
    let dab = dez * ab + aez * bd + bez * da;
    let alift = aex * aex + aey * aey + aez * aez;
    let blift = bex * bex + bey * bey + bez * bez;
    let clift = cex * cex + cey * cey + cez * cez;
    let dlift = dex * dex + dey * dey + dez * dez;
    dlift * abc - clift * dab + (blift * cda - alift * bcd)
}

#[inline]
pub fn insphere_exact(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3], pe: [f64; 3]) -> f64 {
    let mut temp8a = [0.; 8];
    let mut temp8b = [0.; 8];
    let mut temp16 = [0.; 16];
    let mut abc = [0.; 24];
    let mut bcd = [0.; 24];
    let mut cde = [0.; 24];
    let mut dea = [0.; 24];
    let mut eab = [0.; 24];
    let mut abd = [0.; 24];
    let mut bce = [0.; 24];
    let mut cda = [0.; 24];
    let mut deb = [0.; 24];
    let mut eac = [0.; 24];
    let mut temp48a = [0.; 48];
    let mut temp48b = [0.; 48];
    let mut abcd = [0.; 96];
    let mut bcde = [0.; 96];
    let mut cdea = [0.; 96];
    let mut deab = [0.; 96];
    let mut eabc = [0.; 96];
    let mut temp192 = [0.; 192];
    let mut det384x = [0.; 384];
    let mut det384y = [0.; 384];
    let mut det384z = [0.; 384];
    let mut detxy = [0.; 768];
    let mut adet = [0.; 1152];
    let mut bdet = [0.; 1152];
    let mut cdet = [0.; 1152];
    let mut ddet = [0.; 1152];
    let mut edet = [0.; 1152];
    let mut abdet = [0.; 2304];
    let mut cddet = [0.; 2304];
    let mut cdedet = [0.; 3456];
    let [axby0, axby1] = two_product(pa[0], pb[1]);
    let [bxay0, bxay1] = two_product(pb[0], pa[1]);
    let ab = two_two_diff(axby1, axby0, bxay1, bxay0);
    let [bxcy0, bxcy1] = two_product(pb[0], pc[1]);
    let [cxby0, cxby1] = two_product(pc[0], pb[1]);
    let bc = two_two_diff(bxcy1, bxcy0, cxby1, cxby0);
    let [cxdy0, cxdy1] = two_product(pc[0], pd[1]);
    let [dxcy0, dxcy1] = two_product(pd[0], pc[1]);
    let cd = two_two_diff(cxdy1, cxdy0, dxcy1, dxcy0);
    let [dxey0, dxey1] = two_product(pd[0], pe[1]);
    let [exdy0, exdy1] = two_product(pe[0], pd[1]);
    let de = two_two_diff(dxey1, dxey0, exdy1, exdy0);
    let [exay0, exay1] = two_product(pe[0], pa[1]);
    let [axey0, axey1] = two_product(pa[0], pe[1]);
    let ea = two_two_diff(exay1, exay0, axey1, axey0);
    let [axcy0, axcy1] = two_product(pa[0], pc[1]);
    let [cxay0, cxay1] = two_product(pc[0], pa[1]);
    let ac = two_two_diff(axcy1, axcy0, cxay1, cxay0);
    let [bxdy0, bxdy1] = two_product(pb[0], pd[1]);
    let [dxby0, dxby1] = two_product(pd[0], pb[1]);
    let bd = two_two_diff(bxdy1, bxdy0, dxby1, dxby0);
    let [cxey0, cxey1] = two_product(pc[0], pe[1]);
    let [excy0, excy1] = two_product(pe[0], pc[1]);
    let ce = two_two_diff(cxey1, cxey0, excy1, excy0);
    let [dxay0, dxay1] = two_product(pd[0], pa[1]);
    let [axdy0, axdy1] = two_product(pa[0], pd[1]);
    let da = two_two_diff(dxay1, dxay0, axdy1, axdy0);
    let [exby0, exby1] = two_product(pe[0], pb[1]);
    let [bxey0, bxey1] = two_product(pb[0], pe[1]);
    let eb = two_two_diff(exby1, exby0, bxey1, bxey0);
    let temp8alen = scale_expansion_zeroelim(&bc, pa[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&ac, -pb[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&ab, pc[2], &mut temp8a);
    let abclen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut abc);
    let temp8alen = scale_expansion_zeroelim(&cd, pb[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&bd, -pc[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&bc, pd[2], &mut temp8a);
    let bcdlen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut bcd);
    let temp8alen = scale_expansion_zeroelim(&de, pc[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&ce, -pd[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&cd, pe[2], &mut temp8a);
    let cdelen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut cde);
    let temp8alen = scale_expansion_zeroelim(&ea, pd[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&da, -pe[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&de, pa[2], &mut temp8a);
    let dealen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut dea);
    let temp8alen = scale_expansion_zeroelim(&ab, pe[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&eb, -pa[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&ea, pb[2], &mut temp8a);
    let eablen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut eab);
    let temp8alen = scale_expansion_zeroelim(&bd, pa[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&da, pb[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&ab, pd[2], &mut temp8a);
    let abdlen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut abd);
    let temp8alen = scale_expansion_zeroelim(&ce, pb[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&eb, pc[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&bc, pe[2], &mut temp8a);
    let bcelen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut bce);
    let temp8alen = scale_expansion_zeroelim(&da, pc[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&ac, pd[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&cd, pa[2], &mut temp8a);
    let cdalen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut cda);
    let temp8alen = scale_expansion_zeroelim(&eb, pd[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&bd, pe[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&de, pb[2], &mut temp8a);
    let deblen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut deb);
    let temp8alen = scale_expansion_zeroelim(&ac, pe[2], &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&ce, pa[2], &mut temp8b);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp8alen = scale_expansion_zeroelim(&ea, pc[2], &mut temp8a);
    let eaclen = fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp16[..temp16len], &mut eac);
    let temp48alen = fast_expansion_sum_zeroelim(&cde[..cdelen], &bce[..bcelen], &mut temp48a);
    let temp48blen = fast_expansion_sum_zeroelim(&deb[..deblen], &bcd[..bcdlen], &mut temp48b);
    temp48b[..temp48blen].iter_mut().for_each(|x| *x = -*x);
    let bcdelen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut bcde);
    let xlen = scale_expansion_zeroelim(&bcde[..bcdelen], pa[0], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..xlen], pa[0], &mut det384x);
    let ylen = scale_expansion_zeroelim(&bcde[..bcdelen], pa[1], &mut temp192);
    let ylen = scale_expansion_zeroelim(&temp192[..ylen], pa[1], &mut det384y);
    let zlen = scale_expansion_zeroelim(&bcde[..bcdelen], pa[2], &mut temp192);
    let zlen = scale_expansion_zeroelim(&temp192[..zlen], pa[2], &mut det384z);
    let xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let alen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut adet);
    let temp48alen = fast_expansion_sum_zeroelim(&dea[..dealen], &cda[..cdalen], &mut temp48a);
    let temp48blen = fast_expansion_sum_zeroelim(&eac[..eaclen], &cde[..cdelen], &mut temp48b);
    temp48b[..temp48blen].iter_mut().for_each(|x| *x = -*x);
    let cdealen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut cdea);
    let xlen = scale_expansion_zeroelim(&cdea[..cdealen], pb[0], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..xlen], pb[0], &mut det384x);
    let ylen = scale_expansion_zeroelim(&cdea[..cdealen], pb[1], &mut temp192);
    let ylen = scale_expansion_zeroelim(&temp192[..ylen], pb[1], &mut det384y);
    let zlen = scale_expansion_zeroelim(&cdea[..cdealen], pb[2], &mut temp192);
    let zlen = scale_expansion_zeroelim(&temp192[..zlen], pb[2], &mut det384z);
    let xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let blen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut bdet);
    let temp48alen = fast_expansion_sum_zeroelim(&eab[..eablen], &deb[..deblen], &mut temp48a);
    let temp48blen = fast_expansion_sum_zeroelim(&abd[..abdlen], &dea[..dealen], &mut temp48b);
    temp48b[..temp48blen].iter_mut().for_each(|x| *x = -*x);
    let deablen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut deab);
    let xlen = scale_expansion_zeroelim(&deab[..deablen], pc[0], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..xlen], pc[0], &mut det384x);
    let ylen = scale_expansion_zeroelim(&deab[..deablen], pc[1], &mut temp192);
    let ylen = scale_expansion_zeroelim(&temp192[..ylen], pc[1], &mut det384y);
    let zlen = scale_expansion_zeroelim(&deab[..deablen], pc[2], &mut temp192);
    let zlen = scale_expansion_zeroelim(&temp192[..zlen], pc[2], &mut det384z);
    let xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let clen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut cdet);
    let temp48alen = fast_expansion_sum_zeroelim(&abc[..abclen], &eac[..eaclen], &mut temp48a);
    let temp48blen = fast_expansion_sum_zeroelim(&bce[..bcelen], &eab[..eablen], &mut temp48b);
    temp48b[..temp48blen].iter_mut().for_each(|x| *x = -*x);
    let eabclen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut eabc);
    let xlen = scale_expansion_zeroelim(&eabc[..eabclen], pd[0], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..xlen], pd[0], &mut det384x);
    let ylen = scale_expansion_zeroelim(&eabc[..eabclen], pd[1], &mut temp192);
    let ylen = scale_expansion_zeroelim(&temp192[..ylen], pd[1], &mut det384y);
    let zlen = scale_expansion_zeroelim(&eabc[..eabclen], pd[2], &mut temp192);
    let zlen = scale_expansion_zeroelim(&temp192[..zlen], pd[2], &mut det384z);
    let xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let dlen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut ddet);
    let temp48alen = fast_expansion_sum_zeroelim(&bcd[..bcdlen], &abd[..abdlen], &mut temp48a);
    let temp48blen = fast_expansion_sum_zeroelim(&cda[..cdalen], &abc[..abclen], &mut temp48b);
    temp48b[..temp48blen].iter_mut().for_each(|x| *x = -*x);
    let abcdlen =
        fast_expansion_sum_zeroelim(&temp48a[..temp48alen], &temp48b[..temp48blen], &mut abcd);
    let xlen = scale_expansion_zeroelim(&abcd[..abcdlen], pe[0], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..xlen], pe[0], &mut det384x);
    let ylen = scale_expansion_zeroelim(&abcd[..abcdlen], pe[1], &mut temp192);
    let ylen = scale_expansion_zeroelim(&temp192[..ylen], pe[1], &mut det384y);
    let zlen = scale_expansion_zeroelim(&abcd[..abcdlen], pe[2], &mut temp192);
    let zlen = scale_expansion_zeroelim(&temp192[..zlen], pe[2], &mut det384z);
    let xylen = fast_expansion_sum_zeroelim(&det384x[..xlen], &det384y[..ylen], &mut detxy);
    let elen = fast_expansion_sum_zeroelim(&detxy[..xylen], &det384z[..zlen], &mut edet);
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let cdlen = fast_expansion_sum_zeroelim(&cdet[..clen], &ddet[..dlen], &mut cddet);
    let cdelen = fast_expansion_sum_zeroelim(&cddet[..cdlen], &edet[..elen], &mut cdedet);
    let mut deter = [0.; 5760];
    let deterlen = fast_expansion_sum_zeroelim(&abdet[..ablen], &cdedet[..cdelen], &mut deter);
    deter[deterlen - 1]
}

#[inline]
pub fn insphere_slow(pa: [f64; 3], pb: [f64; 3], pc: [f64; 3], pd: [f64; 3], pe: [f64; 3]) -> f64 {
    let mut ab = [0.; 16];
    let mut bc = [0.; 16];
    let mut cd = [0.; 16];
    let mut da = [0.; 16];
    let mut ac = [0.; 16];
    let mut bd = [0.; 16];
    let mut temp32a = [0.; 32];
    let mut temp32b = [0.; 32];
    let mut temp64a = [0.; 64];
    let mut temp64b = [0.; 64];
    let mut temp64c = [0.; 64];
    let mut temp128 = [0.; 128];
    let mut temp192 = [0.; 192];
    let mut detx = [0.; 384];
    let mut detxx = [0.; 768];
    let mut detxt = [0.; 384];
    let mut detxxt = [0.; 768];
    let mut detxtxt = [0.; 768];
    let mut x1 = [0.; 1536];
    let mut x2 = [0.; 2304];
    let mut dety = [0.; 384];
    let mut detyy = [0.; 768];
    let mut detyt = [0.; 384];
    let mut detyyt = [0.; 768];
    let mut detytyt = [0.; 768];
    let mut y1 = [0.; 1536];
    let mut y2 = [0.; 2304];
    let mut detz = [0.; 384];
    let mut detzz = [0.; 768];
    let mut detzt = [0.; 384];
    let mut detzzt = [0.; 768];
    let mut detztzt = [0.; 768];
    let mut z1 = [0.; 1536];
    let mut z2 = [0.; 2304];
    let mut detxy = [0.; 4608];
    let mut adet = [0.; 6912];
    let mut bdet = [0.; 6912];
    let mut cdet = [0.; 6912];
    let mut ddet = [0.; 6912];
    let mut abdet = [0.; 13824];
    let mut cddet = [0.; 13824];
    let mut deter = [0.; 27648];
    let [aextail, aex] = two_diff(pa[0], pe[0]);
    let [aeytail, aey] = two_diff(pa[1], pe[1]);
    let [aeztail, aez] = two_diff(pa[2], pe[2]);
    let [bextail, bex] = two_diff(pb[0], pe[0]);
    let [beytail, bey] = two_diff(pb[1], pe[1]);
    let [beztail, bez] = two_diff(pb[2], pe[2]);
    let [cextail, cex] = two_diff(pc[0], pe[0]);
    let [ceytail, cey] = two_diff(pc[1], pe[1]);
    let [ceztail, cez] = two_diff(pc[2], pe[2]);
    let [dextail, dex] = two_diff(pd[0], pe[0]);
    let [deytail, dey] = two_diff(pd[1], pe[1]);
    let [deztail, dez] = two_diff(pd[2], pe[2]);
    let axby = two_two_product(aex, aextail, bey, beytail);
    let negate = -aey;
    let negatetail = -aeytail;
    let bxay = two_two_product(bex, bextail, negate, negatetail);
    let ablen = fast_expansion_sum_zeroelim(&axby, &bxay, &mut ab);
    let bxcy = two_two_product(bex, bextail, cey, ceytail);
    let negate = -bey;
    let negatetail = -beytail;
    let cxby = two_two_product(cex, cextail, negate, negatetail);
    let bclen = fast_expansion_sum_zeroelim(&bxcy, &cxby, &mut bc);
    let cxdy = two_two_product(cex, cextail, dey, deytail);
    let negate = -cey;
    let negatetail = -ceytail;
    let dxcy = two_two_product(dex, dextail, negate, negatetail);
    let cdlen = fast_expansion_sum_zeroelim(&cxdy, &dxcy, &mut cd);
    let dxay = two_two_product(dex, dextail, aey, aeytail);
    let negate = -dey;
    let negatetail = -deytail;
    let axdy = two_two_product(aex, aextail, negate, negatetail);
    let dalen = fast_expansion_sum_zeroelim(&dxay, &axdy, &mut da);
    let axcy = two_two_product(aex, aextail, cey, ceytail);
    let negate = -aey;
    let negatetail = -aeytail;
    let cxay = two_two_product(cex, cextail, negate, negatetail);
    let aclen = fast_expansion_sum_zeroelim(&axcy, &cxay, &mut ac);
    let bxdy = two_two_product(bex, bextail, dey, deytail);
    let negate = -bey;
    let negatetail = -beytail;
    let dxby = two_two_product(dex, dextail, negate, negatetail);
    let bdlen = fast_expansion_sum_zeroelim(&bxdy, &dxby, &mut bd);
    let temp32alen = scale_expansion_zeroelim(&cd[..cdlen], -bez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&cd[..cdlen], -beztail, &mut temp32b);
    let temp64alen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64a);
    let temp32alen = scale_expansion_zeroelim(&bd[..bdlen], cez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&bd[..bdlen], ceztail, &mut temp32b);
    let temp64blen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64b);
    let temp32alen = scale_expansion_zeroelim(&bc[..bclen], -dez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&bc[..bclen], -deztail, &mut temp32b);
    let temp64clen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64c);
    let temp128len =
        fast_expansion_sum_zeroelim(&temp64a[..temp64alen], &temp64b[..temp64blen], &mut temp128);
    let temp192len =
        fast_expansion_sum_zeroelim(&temp64c[..temp64clen], &temp128[..temp128len], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..temp192len], aex, &mut detx);
    let xxlen = scale_expansion_zeroelim(&detx[..xlen], aex, &mut detxx);
    let xtlen = scale_expansion_zeroelim(&temp192[..temp192len], aextail, &mut detxt);
    let xxtlen = scale_expansion_zeroelim(&detxt[..xtlen], aex, &mut detxxt);
    detxxt[..xxtlen].iter_mut().for_each(|x| *x *= 2.0);
    let xtxtlen = scale_expansion_zeroelim(&detxt[..xtlen], aextail, &mut detxtxt);
    let x1len = fast_expansion_sum_zeroelim(&detxx[..xxlen], &detxxt[..xxtlen], &mut x1);
    let x2len = fast_expansion_sum_zeroelim(&x1[..x1len], &detxtxt[..xtxtlen], &mut x2);
    let ylen = scale_expansion_zeroelim(&temp192[..temp192len], aey, &mut dety);
    let yylen = scale_expansion_zeroelim(&dety[..ylen], aey, &mut detyy);
    let ytlen = scale_expansion_zeroelim(&temp192[..temp192len], aeytail, &mut detyt);
    let yytlen = scale_expansion_zeroelim(&detyt[..ytlen], aey, &mut detyyt);
    detyyt[..yytlen].iter_mut().for_each(|x| *x *= 2.0);
    let ytytlen = scale_expansion_zeroelim(&detyt[..ytlen], aeytail, &mut detytyt);
    let y1len = fast_expansion_sum_zeroelim(&detyy[..yylen], &detyyt[..yytlen], &mut y1);
    let y2len = fast_expansion_sum_zeroelim(&y1[..y1len], &detytyt[..ytytlen], &mut y2);
    let zlen = scale_expansion_zeroelim(&temp192[..temp192len], aez, &mut detz);
    let zzlen = scale_expansion_zeroelim(&detz[..zlen], aez, &mut detzz);
    let ztlen = scale_expansion_zeroelim(&temp192[..temp192len], aeztail, &mut detzt);
    let zztlen = scale_expansion_zeroelim(&detzt[..ztlen], aez, &mut detzzt);
    detzzt[..zztlen].iter_mut().for_each(|x| *x *= 2.0);
    let ztztlen = scale_expansion_zeroelim(&detzt[..ztlen], aeztail, &mut detztzt);
    let z1len = fast_expansion_sum_zeroelim(&detzz[..zzlen], &detzzt[..zztlen], &mut z1);
    let z2len = fast_expansion_sum_zeroelim(&z1[..z1len], &detztzt[..ztztlen], &mut z2);
    let xylen = fast_expansion_sum_zeroelim(&x2[..x2len], &y2[..y2len], &mut detxy);
    let alen = fast_expansion_sum_zeroelim(&z2[..z2len], &detxy[..xylen], &mut adet);
    let temp32alen = scale_expansion_zeroelim(&da[..dalen], cez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&da[..dalen], ceztail, &mut temp32b);
    let temp64alen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64a);
    let temp32alen = scale_expansion_zeroelim(&ac[..aclen], dez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&ac[..aclen], deztail, &mut temp32b);
    let temp64blen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64b);
    let temp32alen = scale_expansion_zeroelim(&cd[..cdlen], aez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&cd[..cdlen], aeztail, &mut temp32b);
    let temp64clen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64c);
    let temp128len =
        fast_expansion_sum_zeroelim(&temp64a[..temp64alen], &temp64b[..temp64blen], &mut temp128);
    let temp192len =
        fast_expansion_sum_zeroelim(&temp64c[..temp64clen], &temp128[..temp128len], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..temp192len], bex, &mut detx);
    let xxlen = scale_expansion_zeroelim(&detx[..xlen], bex, &mut detxx);
    let xtlen = scale_expansion_zeroelim(&temp192[..temp192len], bextail, &mut detxt);
    let xxtlen = scale_expansion_zeroelim(&detxt[..xtlen], bex, &mut detxxt);
    detxxt[..xxtlen].iter_mut().for_each(|x| *x *= 2.0);
    let xtxtlen = scale_expansion_zeroelim(&detxt[..xtlen], bextail, &mut detxtxt);
    let x1len = fast_expansion_sum_zeroelim(&detxx[..xxlen], &detxxt[..xxtlen], &mut x1);
    let x2len = fast_expansion_sum_zeroelim(&x1[..x1len], &detxtxt[..xtxtlen], &mut x2);
    let ylen = scale_expansion_zeroelim(&temp192[..temp192len], bey, &mut dety);
    let yylen = scale_expansion_zeroelim(&dety[..ylen], bey, &mut detyy);
    let ytlen = scale_expansion_zeroelim(&temp192[..temp192len], beytail, &mut detyt);
    let yytlen = scale_expansion_zeroelim(&detyt[..ytlen], bey, &mut detyyt);
    detyyt[..yytlen].iter_mut().for_each(|x| *x *= 2.0);
    let ytytlen = scale_expansion_zeroelim(&detyt[..ytlen], beytail, &mut detytyt);
    let y1len = fast_expansion_sum_zeroelim(&detyy[..yylen], &detyyt[..yytlen], &mut y1);
    let y2len = fast_expansion_sum_zeroelim(&y1[..y1len], &detytyt[..ytytlen], &mut y2);
    let zlen = scale_expansion_zeroelim(&temp192[..temp192len], bez, &mut detz);
    let zzlen = scale_expansion_zeroelim(&detz[..zlen], bez, &mut detzz);
    let ztlen = scale_expansion_zeroelim(&temp192[..temp192len], beztail, &mut detzt);
    let zztlen = scale_expansion_zeroelim(&detzt[..ztlen], bez, &mut detzzt);
    detzzt[..zztlen].iter_mut().for_each(|x| *x *= 2.0);
    let ztztlen = scale_expansion_zeroelim(&detzt[..ztlen], beztail, &mut detztzt);
    let z1len = fast_expansion_sum_zeroelim(&detzz[..zzlen], &detzzt[..zztlen], &mut z1);
    let z2len = fast_expansion_sum_zeroelim(&z1[..z1len], &detztzt[..ztztlen], &mut z2);
    let xylen = fast_expansion_sum_zeroelim(&x2[..x2len], &y2[..y2len], &mut detxy);
    let blen = fast_expansion_sum_zeroelim(&z2[..z2len], &detxy[..xylen], &mut bdet);
    let temp32alen = scale_expansion_zeroelim(&ab[..ablen], -dez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&ab[..ablen], -deztail, &mut temp32b);
    let temp64alen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64a);
    let temp32alen = scale_expansion_zeroelim(&bd[..bdlen], -aez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&bd[..bdlen], -aeztail, &mut temp32b);
    let temp64blen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64b);
    let temp32alen = scale_expansion_zeroelim(&da[..dalen], -bez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&da[..dalen], -beztail, &mut temp32b);
    let temp64clen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64c);
    let temp128len =
        fast_expansion_sum_zeroelim(&temp64a[..temp64alen], &temp64b[..temp64blen], &mut temp128);
    let temp192len =
        fast_expansion_sum_zeroelim(&temp64c[..temp64clen], &temp128[..temp128len], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..temp192len], cex, &mut detx);
    let xxlen = scale_expansion_zeroelim(&detx[..xlen], cex, &mut detxx);
    let xtlen = scale_expansion_zeroelim(&temp192[..temp192len], cextail, &mut detxt);
    let xxtlen = scale_expansion_zeroelim(&detxt[..xtlen], cex, &mut detxxt);
    detxxt[..xxtlen].iter_mut().for_each(|x| *x *= 2.0);
    let xtxtlen = scale_expansion_zeroelim(&detxt[..xtlen], cextail, &mut detxtxt);
    let x1len = fast_expansion_sum_zeroelim(&detxx[..xxlen], &detxxt[..xxtlen], &mut x1);
    let x2len = fast_expansion_sum_zeroelim(&x1[..x1len], &detxtxt[..xtxtlen], &mut x2);
    let ylen = scale_expansion_zeroelim(&temp192[..temp192len], cey, &mut dety);
    let yylen = scale_expansion_zeroelim(&dety[..ylen], cey, &mut detyy);
    let ytlen = scale_expansion_zeroelim(&temp192[..temp192len], ceytail, &mut detyt);
    let yytlen = scale_expansion_zeroelim(&detyt[..ytlen], cey, &mut detyyt);
    detyyt[..yytlen].iter_mut().for_each(|x| *x *= 2.0);
    let ytytlen = scale_expansion_zeroelim(&detyt[..ytlen], ceytail, &mut detytyt);
    let y1len = fast_expansion_sum_zeroelim(&detyy[..yylen], &detyyt[..yytlen], &mut y1);
    let y2len = fast_expansion_sum_zeroelim(&y1[..y1len], &detytyt[..ytytlen], &mut y2);
    let zlen = scale_expansion_zeroelim(&temp192[..temp192len], cez, &mut detz);
    let zzlen = scale_expansion_zeroelim(&detz[..zlen], cez, &mut detzz);
    let ztlen = scale_expansion_zeroelim(&temp192[..temp192len], ceztail, &mut detzt);
    let zztlen = scale_expansion_zeroelim(&detzt[..ztlen], cez, &mut detzzt);
    detzzt[..zztlen].iter_mut().for_each(|x| *x *= 2.0);
    let ztztlen = scale_expansion_zeroelim(&detzt[..ztlen], ceztail, &mut detztzt);
    let z1len = fast_expansion_sum_zeroelim(&detzz[..zzlen], &detzzt[..zztlen], &mut z1);
    let z2len = fast_expansion_sum_zeroelim(&z1[..z1len], &detztzt[..ztztlen], &mut z2);
    let xylen = fast_expansion_sum_zeroelim(&x2[..x2len], &y2[..y2len], &mut detxy);
    let clen = fast_expansion_sum_zeroelim(&z2[..z2len], &detxy[..xylen], &mut cdet);
    let temp32alen = scale_expansion_zeroelim(&bc[..bclen], aez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&bc[..bclen], aeztail, &mut temp32b);
    let temp64alen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64a);
    let temp32alen = scale_expansion_zeroelim(&ac[..aclen], -bez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&ac[..aclen], -beztail, &mut temp32b);
    let temp64blen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64b);
    let temp32alen = scale_expansion_zeroelim(&ab[..ablen], cez, &mut temp32a);
    let temp32blen = scale_expansion_zeroelim(&ab[..ablen], ceztail, &mut temp32b);
    let temp64clen =
        fast_expansion_sum_zeroelim(&temp32a[..temp32alen], &temp32b[..temp32blen], &mut temp64c);
    let temp128len =
        fast_expansion_sum_zeroelim(&temp64a[..temp64alen], &temp64b[..temp64blen], &mut temp128);
    let temp192len =
        fast_expansion_sum_zeroelim(&temp64c[..temp64clen], &temp128[..temp128len], &mut temp192);
    let xlen = scale_expansion_zeroelim(&temp192[..temp192len], dex, &mut detx);
    let xxlen = scale_expansion_zeroelim(&detx[..xlen], dex, &mut detxx);
    let xtlen = scale_expansion_zeroelim(&temp192[..temp192len], dextail, &mut detxt);
    let xxtlen = scale_expansion_zeroelim(&detxt[..xtlen], dex, &mut detxxt);
    detxxt[..xxtlen].iter_mut().for_each(|x| *x *= 2.0);
    let xtxtlen = scale_expansion_zeroelim(&detxt[..xtlen], dextail, &mut detxtxt);
    let x1len = fast_expansion_sum_zeroelim(&detxx[..xxlen], &detxxt[..xxtlen], &mut x1);
    let x2len = fast_expansion_sum_zeroelim(&x1[..x1len], &detxtxt[..xtxtlen], &mut x2);
    let ylen = scale_expansion_zeroelim(&temp192[..temp192len], dey, &mut dety);
    let yylen = scale_expansion_zeroelim(&dety[..ylen], dey, &mut detyy);
    let ytlen = scale_expansion_zeroelim(&temp192[..temp192len], deytail, &mut detyt);
    let yytlen = scale_expansion_zeroelim(&detyt[..ytlen], dey, &mut detyyt);
    detyyt[..yytlen].iter_mut().for_each(|x| *x *= 2.0);
    let ytytlen = scale_expansion_zeroelim(&detyt[..ytlen], deytail, &mut detytyt);
    let y1len = fast_expansion_sum_zeroelim(&detyy[..yylen], &detyyt[..yytlen], &mut y1);
    let y2len = fast_expansion_sum_zeroelim(&y1[..y1len], &detytyt[..ytytlen], &mut y2);
    let zlen = scale_expansion_zeroelim(&temp192[..temp192len], dez, &mut detz);
    let zzlen = scale_expansion_zeroelim(&detz[..zlen], dez, &mut detzz);
    let ztlen = scale_expansion_zeroelim(&temp192[..temp192len], deztail, &mut detzt);
    let zztlen = scale_expansion_zeroelim(&detzt[..ztlen], dez, &mut detzzt);
    detzzt[..zztlen].iter_mut().for_each(|x| *x *= 2.0);
    let ztztlen = scale_expansion_zeroelim(&detzt[..ztlen], deztail, &mut detztzt);
    let z1len = fast_expansion_sum_zeroelim(&detzz[..zzlen], &detzzt[..zztlen], &mut z1);
    let z2len = fast_expansion_sum_zeroelim(&z1[..z1len], &detztzt[..ztztlen], &mut z2);
    let xylen = fast_expansion_sum_zeroelim(&x2[..x2len], &y2[..y2len], &mut detxy);
    let dlen = fast_expansion_sum_zeroelim(&z2[..z2len], &detxy[..xylen], &mut ddet);
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let cdlen = fast_expansion_sum_zeroelim(&cdet[..clen], &ddet[..dlen], &mut cddet);
    let deterlen = fast_expansion_sum_zeroelim(&abdet[..ablen], &cddet[..cdlen], &mut deter);
    deter[deterlen - 1]
}

#[inline]
pub fn insphereadapt(
    pa: [f64; 3],
    pb: [f64; 3],
    pc: [f64; 3],
    pd: [f64; 3],
    pe: [f64; 3],
    permanent: f64,
) -> f64 {
    let mut temp8a: [f64; 8] = [0.; 8];
    let mut temp8b: [f64; 8] = [0.; 8];
    let mut temp8c: [f64; 8] = [0.; 8];
    let mut temp16: [f64; 16] = [0.; 16];
    let mut temp24: [f64; 24] = [0.; 24];
    let mut temp48: [f64; 48] = [0.; 48];
    let mut xdet: [f64; 96] = [0.; 96];
    let mut ydet: [f64; 96] = [0.; 96];
    let mut zdet: [f64; 96] = [0.; 96];
    let mut xydet: [f64; 192] = [0.; 192];
    let mut adet: [f64; 288] = [0.; 288];
    let mut bdet: [f64; 288] = [0.; 288];
    let mut cdet: [f64; 288] = [0.; 288];
    let mut ddet: [f64; 288] = [0.; 288];
    let mut abdet: [f64; 576] = [0.; 576];
    let mut cddet: [f64; 576] = [0.; 576];
    let mut fin1: [f64; 1152] = [0.; 1152];
    let aex = pa[0] - pe[0];
    let bex = pb[0] - pe[0];
    let cex = pc[0] - pe[0];
    let dex = pd[0] - pe[0];
    let aey = pa[1] - pe[1];
    let bey = pb[1] - pe[1];
    let cey = pc[1] - pe[1];
    let dey = pd[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bez = pb[2] - pe[2];
    let cez = pc[2] - pe[2];
    let dez = pd[2] - pe[2];
    let [aexbey0, aexbey1] = two_product(aex, bey);
    let [bexaey0, bexaey1] = two_product(bex, aey);
    let ab = two_two_diff(aexbey1, aexbey0, bexaey1, bexaey0);
    let [bexcey0, bexcey1] = two_product(bex, cey);
    let [cexbey0, cexbey1] = two_product(cex, bey);
    let bc = two_two_diff(bexcey1, bexcey0, cexbey1, cexbey0);
    let [cexdey0, cexdey1] = two_product(cex, dey);
    let [dexcey0, dexcey1] = two_product(dex, cey);
    let cd = two_two_diff(cexdey1, cexdey0, dexcey1, dexcey0);
    let [dexaey0, dexaey1] = two_product(dex, aey);
    let [aexdey0, aexdey1] = two_product(aex, dey);
    let da = two_two_diff(dexaey1, dexaey0, aexdey1, aexdey0);
    let [aexcey0, aexcey1] = two_product(aex, cey);
    let [cexaey0, cexaey1] = two_product(cex, aey);
    let ac = two_two_diff(aexcey1, aexcey0, cexaey1, cexaey0);
    let [bexdey0, bexdey1] = two_product(bex, dey);
    let [dexbey0, dexbey1] = two_product(dex, bey);
    let bd = two_two_diff(bexdey1, bexdey0, dexbey1, dexbey0);
    let temp8alen = scale_expansion_zeroelim(&cd, bez, &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&bd, -cez, &mut temp8b);
    let temp8clen = scale_expansion_zeroelim(&bc, dez, &mut temp8c);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp24len =
        fast_expansion_sum_zeroelim(&temp8c[..temp8clen], &temp16[..temp16len], &mut temp24);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], aex, &mut temp48);
    let xlen = scale_expansion_zeroelim(&temp48[..temp48len], -aex, &mut xdet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], aey, &mut temp48);
    let ylen = scale_expansion_zeroelim(&temp48[..temp48len], -aey, &mut ydet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], aez, &mut temp48);
    let zlen = scale_expansion_zeroelim(&temp48[..temp48len], -aez, &mut zdet);
    let xylen = fast_expansion_sum_zeroelim(&xdet[..xlen], &ydet[..ylen], &mut xydet);
    let alen = fast_expansion_sum_zeroelim(&xydet[..xylen], &zdet[..zlen], &mut adet);
    let temp8alen = scale_expansion_zeroelim(&da, cez, &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&ac, dez, &mut temp8b);
    let temp8clen = scale_expansion_zeroelim(&cd, aez, &mut temp8c);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp24len =
        fast_expansion_sum_zeroelim(&temp8c[..temp8clen], &temp16[..temp16len], &mut temp24);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], bex, &mut temp48);
    let xlen = scale_expansion_zeroelim(&temp48[..temp48len], bex, &mut xdet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], bey, &mut temp48);
    let ylen = scale_expansion_zeroelim(&temp48[..temp48len], bey, &mut ydet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], bez, &mut temp48);
    let zlen = scale_expansion_zeroelim(&temp48[..temp48len], bez, &mut zdet);
    let xylen = fast_expansion_sum_zeroelim(&xdet[..xlen], &ydet[..ylen], &mut xydet);
    let blen = fast_expansion_sum_zeroelim(&xydet[..xylen], &zdet[..zlen], &mut bdet);
    let temp8alen = scale_expansion_zeroelim(&ab, dez, &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&bd, aez, &mut temp8b);
    let temp8clen = scale_expansion_zeroelim(&da, bez, &mut temp8c);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp24len =
        fast_expansion_sum_zeroelim(&temp8c[..temp8clen], &temp16[..temp16len], &mut temp24);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], cex, &mut temp48);
    let xlen = scale_expansion_zeroelim(&temp48[..temp48len], -cex, &mut xdet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], cey, &mut temp48);
    let ylen = scale_expansion_zeroelim(&temp48[..temp48len], -cey, &mut ydet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], cez, &mut temp48);
    let zlen = scale_expansion_zeroelim(&temp48[..temp48len], -cez, &mut zdet);
    let xylen = fast_expansion_sum_zeroelim(&xdet[..xlen], &ydet[..ylen], &mut xydet);
    let clen = fast_expansion_sum_zeroelim(&xydet[..xylen], &zdet[..zlen], &mut cdet);
    let temp8alen = scale_expansion_zeroelim(&bc, aez, &mut temp8a);
    let temp8blen = scale_expansion_zeroelim(&ac, -bez, &mut temp8b);
    let temp8clen = scale_expansion_zeroelim(&ab, cez, &mut temp8c);
    let temp16len =
        fast_expansion_sum_zeroelim(&temp8a[..temp8alen], &temp8b[..temp8blen], &mut temp16);
    let temp24len =
        fast_expansion_sum_zeroelim(&temp8c[..temp8clen], &temp16[..temp16len], &mut temp24);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], dex, &mut temp48);
    let xlen = scale_expansion_zeroelim(&temp48[..temp48len], dex, &mut xdet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], dey, &mut temp48);
    let ylen = scale_expansion_zeroelim(&temp48[..temp48len], dey, &mut ydet);
    let temp48len = scale_expansion_zeroelim(&temp24[..temp24len], dez, &mut temp48);
    let zlen = scale_expansion_zeroelim(&temp48[..temp48len], dez, &mut zdet);
    let xylen = fast_expansion_sum_zeroelim(&xdet[..xlen], &ydet[..ylen], &mut xydet);
    let dlen = fast_expansion_sum_zeroelim(&xydet[..xylen], &zdet[..zlen], &mut ddet);
    let ablen = fast_expansion_sum_zeroelim(&adet[..alen], &bdet[..blen], &mut abdet);
    let cdlen = fast_expansion_sum_zeroelim(&cdet[..clen], &ddet[..dlen], &mut cddet);
    let finlength = fast_expansion_sum_zeroelim(&abdet[..ablen], &cddet[..cdlen], &mut fin1);
    let mut det: f64 = fin1[..finlength].iter().sum();
    let errbound = PARAMS.isperrbound_b * permanent;
    if det >= errbound || -det >= errbound {
        return det;
    }
    let aextail = two_diff_tail(pa[0], pe[0], aex);
    let aeytail = two_diff_tail(pa[1], pe[1], aey);
    let aeztail = two_diff_tail(pa[2], pe[2], aez);
    let bextail = two_diff_tail(pb[0], pe[0], bex);
    let beytail = two_diff_tail(pb[1], pe[1], bey);
    let beztail = two_diff_tail(pb[2], pe[2], bez);
    let cextail = two_diff_tail(pc[0], pe[0], cex);
    let ceytail = two_diff_tail(pc[1], pe[1], cey);
    let ceztail = two_diff_tail(pc[2], pe[2], cez);
    let dextail = two_diff_tail(pd[0], pe[0], dex);
    let deytail = two_diff_tail(pd[1], pe[1], dey);
    let deztail = two_diff_tail(pd[2], pe[2], dez);
    if aextail == 0.0
        && aeytail == 0.0
        && aeztail == 0.0
        && bextail == 0.0
        && beytail == 0.0
        && beztail == 0.0
        && cextail == 0.0
        && ceytail == 0.0
        && ceztail == 0.0
        && dextail == 0.0
        && deytail == 0.0
        && deztail == 0.0
    {
        return det;
    }
    let errbound = PARAMS.isperrbound_c * permanent + PARAMS.resulterrbound * abs(det);
    let abeps = aex * beytail + bey * aextail - (aey * bextail + bex * aeytail);
    let bceps = bex * ceytail + cey * bextail - (bey * cextail + cex * beytail);
    let cdeps = cex * deytail + dey * cextail - (cey * dextail + dex * ceytail);
    let daeps = dex * aeytail + aey * dextail - (dey * aextail + aex * deytail);
    let aceps = aex * ceytail + cey * aextail - (aey * cextail + cex * aeytail);
    let bdeps = bex * deytail + dey * bextail - (bey * dextail + dex * beytail);
    det += (bex * bex + bey * bey + bez * bez)
        * (cez * daeps
            + dez * aceps
            + aez * cdeps
            + (ceztail * da[3] + deztail * ac[3] + aeztail * cd[3]))
        + (dex * dex + dey * dey + dez * dez)
            * (aez * bceps - bez * aceps
                + cez * abeps
                + (aeztail * bc[3] - beztail * ac[3] + ceztail * ab[3]))
        - ((aex * aex + aey * aey + aez * aez)
            * (bez * cdeps - cez * bdeps
                + dez * bceps
                + (beztail * cd[3] - ceztail * bd[3] + deztail * bc[3]))
            + (cex * cex + cey * cey + cez * cez)
                * (dez * abeps
                    + aez * bdeps
                    + bez * daeps
                    + (deztail * ab[3] + aeztail * bd[3] + beztail * da[3])))
        + 2.0
            * ((bex * bextail + bey * beytail + bez * beztail)
                * (cez * da[3] + dez * ac[3] + aez * cd[3])
                + (dex * dextail + dey * deytail + dez * deztail)
                    * (aez * bc[3] - bez * ac[3] + cez * ab[3])
                - ((aex * aextail + aey * aeytail + aez * aeztail)
                    * (bez * cd[3] - cez * bd[3] + dez * bc[3])
                    + (cex * cextail + cey * ceytail + cez * ceztail)
                        * (dez * ab[3] + aez * bd[3] + bez * da[3])));
    if det >= errbound || -det >= errbound {
        return det;
    }
    insphere_exact(pa, pb, pc, pd, pe)
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
    let aex = pa[0] - pe[0];
    let bex = pb[0] - pe[0];
    let cex = pc[0] - pe[0];
    let dex = pd[0] - pe[0];
    let aey = pa[1] - pe[1];
    let bey = pb[1] - pe[1];
    let cey = pc[1] - pe[1];
    let dey = pd[1] - pe[1];
    let aez = pa[2] - pe[2];
    let bez = pb[2] - pe[2];
    let cez = pc[2] - pe[2];
    let dez = pd[2] - pe[2];
    let aexbey = aex * bey;
    let bexaey = bex * aey;
    let ab = aexbey - bexaey;
    let bexcey = bex * cey;
    let cexbey = cex * bey;
    let bc = bexcey - cexbey;
    let cexdey = cex * dey;
    let dexcey = dex * cey;
    let cd = cexdey - dexcey;
    let dexaey = dex * aey;
    let aexdey = aex * dey;
    let da = dexaey - aexdey;
    let aexcey = aex * cey;
    let cexaey = cex * aey;
    let ac = aexcey - cexaey;
    let bexdey = bex * dey;
    let dexbey = dex * bey;
    let bd = bexdey - dexbey;
    let abc = aez * bc - bez * ac + cez * ab;
    let bcd = bez * cd - cez * bd + dez * bc;
    let cda = cez * da + dez * ac + aez * cd;
    let dab = dez * ab + aez * bd + bez * da;
    let alift = aex * aex + aey * aey + aez * aez;
    let blift = bex * bex + bey * bey + bez * bez;
    let clift = cex * cex + cey * cey + cez * cez;
    let dlift = dex * dex + dey * dey + dez * dez;
    let det = dlift * abc - clift * dab + (blift * cda - alift * bcd);
    let aezplus = abs(aez);
    let bezplus = abs(bez);
    let cezplus = abs(cez);
    let dezplus = abs(dez);
    let aexbeyplus = abs(aexbey);
    let bexaeyplus = abs(bexaey);
    let bexceyplus = abs(bexcey);
    let cexbeyplus = abs(cexbey);
    let cexdeyplus = abs(cexdey);
    let dexceyplus = abs(dexcey);
    let dexaeyplus = abs(dexaey);
    let aexdeyplus = abs(aexdey);
    let aexceyplus = abs(aexcey);
    let cexaeyplus = abs(cexaey);
    let bexdeyplus = abs(bexdey);
    let dexbeyplus = abs(dexbey);
    let permanent = ((cexdeyplus + dexceyplus) * bezplus
        + (dexbeyplus + bexdeyplus) * cezplus
        + (bexceyplus + cexbeyplus) * dezplus)
        * alift
        + ((dexaeyplus + aexdeyplus) * cezplus
            + (aexceyplus + cexaeyplus) * dezplus
            + (cexdeyplus + dexceyplus) * aezplus)
            * blift
        + ((aexbeyplus + bexaeyplus) * dezplus
            + (bexdeyplus + dexbeyplus) * aezplus
            + (dexaeyplus + aexdeyplus) * bezplus)
            * clift
        + ((bexceyplus + cexbeyplus) * aezplus
            + (cexaeyplus + aexceyplus) * bezplus
            + (aexbeyplus + bexaeyplus) * cezplus)
            * dlift;
    let errbound = PARAMS.isperrbound_a * permanent;
    if det > errbound || -det > errbound {
        return det;
    }
    insphereadapt(pa, pb, pc, pd, pe, permanent)
}
