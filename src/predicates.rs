#![allow(non_snake_case)] // disable variable naming warnings
#![allow(unused_mut)] // disable unused mut warnings
#![allow(unused_assignments)] // disable unused assignment warnings
#![allow(dead_code)] // allow unused functions for now

use std::mem::uninitialized; // allow uninitialized variables

#[inline(always)]
fn Absolute(a : f64) -> f64 {
    f64::abs(a)
}

pub struct RawGeometryPredicates {
    splitter : f64,
    resulterrbound : f64,
    ccwerrboundA : f64,
    ccwerrboundB : f64,
    ccwerrboundC : f64,
    o3derrboundA : f64,
    o3derrboundB : f64,
    o3derrboundC : f64,
    iccerrboundA : f64,
    iccerrboundB : f64,
    iccerrboundC : f64,
    isperrboundA : f64,
    isperrboundB : f64,
    isperrboundC : f64,
}


impl RawGeometryPredicates {
    pub fn exactinit() -> RawGeometryPredicates {

        let mut half : f64;
        let mut check : f64;
        let mut lastcheck : f64;
        let mut every_other : i32;
        every_other = 1i32;
        half = 0.5f64;
        let mut epsilon = 1.0f64;
        let mut splitter = 1.0f64;
        check = 1.0f64;
        'loop1: loop {
            lastcheck = check;
            epsilon = epsilon * half;
            if every_other != 0 {
                splitter = splitter * 2.0f64;
            }
            every_other = (every_other == 0) as (i32);
            check = 1.0f64 + epsilon;
            if !(check != 1.0f64 && (check != lastcheck)) {
                break;
            }
        }
        let splitter = splitter + 1.0f64;
        let resulterrbound = (3.0f64 + 8.0f64 * epsilon) * epsilon;
        let ccwerrboundA = (3.0f64 + 16.0f64 * epsilon) * epsilon;
        let ccwerrboundB = (2.0f64 + 12.0f64 * epsilon) * epsilon;
        let ccwerrboundC = (9.0f64 + 64.0f64 * epsilon) * epsilon * epsilon;
        let o3derrboundA = (7.0f64 + 56.0f64 * epsilon) * epsilon;
        let o3derrboundB = (3.0f64 + 28.0f64 * epsilon) * epsilon;
        let o3derrboundC = (26.0f64 + 288.0f64 * epsilon) * epsilon * epsilon;
        let iccerrboundA = (10.0f64 + 96.0f64 * epsilon) * epsilon;
        let iccerrboundB = (4.0f64 + 48.0f64 * epsilon) * epsilon;
        let iccerrboundC = (44.0f64 + 576.0f64 * epsilon) * epsilon * epsilon;
        let isperrboundA = (16.0f64 + 224.0f64 * epsilon) * epsilon;
        let isperrboundB = (5.0f64 + 72.0f64 * epsilon) * epsilon;
        let isperrboundC = (71.0f64 + 1408.0f64 * epsilon) * epsilon * epsilon;
        RawGeometryPredicates {
            splitter,
            resulterrbound,
            ccwerrboundA,
            ccwerrboundB,
            ccwerrboundC,
            o3derrboundA,
            o3derrboundB,
            o3derrboundC,
            iccerrboundA,
            iccerrboundB,
            iccerrboundC,
            isperrboundA,
            isperrboundB,
            isperrboundC,
        }
    }

    pub unsafe fn scale_expansion(&self,
        mut elen : i32, mut e : *const f64, mut b : f64, mut h : *mut f64
    ) -> i32 {
        let mut Q : f64;
        let mut sum : f64;
        let mut product1 : f64;
        let mut product0 : f64;
        let mut eindex : i32;
        let mut hindex : i32;
        let mut enow : f64;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        c = self.splitter * b;
        abig = c - b;
        bhi = c - abig;
        blo = b - bhi;
        Q = *e.offset(0isize) * b;
        c = self.splitter * *e.offset(0isize);
        abig = c - *e.offset(0isize);
        ahi = c - abig;
        alo = *e.offset(0isize) - ahi;
        err1 = Q - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        *h.offset(0isize) = alo * blo - err3;
        hindex = 1i32;
        eindex = 1i32;
        'loop1: loop {
            if !(eindex < elen) {
                break;
            }
            enow = *e.offset(eindex as (isize));
            product1 = enow * b;
            c = self.splitter * enow;
            abig = c - enow;
            ahi = c - abig;
            alo = enow - ahi;
            err1 = product1 - ahi * bhi;
            err2 = err1 - alo * bhi;
            err3 = err2 - ahi * blo;
            product0 = alo * blo - err3;
            sum = Q + product0;
            bvirt = sum - Q;
            avirt = sum - bvirt;
            bround = product0 - bvirt;
            around = Q - avirt;
            *h.offset(hindex as (isize)) = around + bround;
            hindex = hindex + 1;
            Q = product1 + sum;
            bvirt = Q - product1;
            avirt = Q - bvirt;
            bround = sum - bvirt;
            around = product1 - avirt;
            *h.offset(hindex as (isize)) = around + bround;
            hindex = hindex + 1;
            eindex = eindex + 1;
        }
        *h.offset(hindex as (isize)) = Q;
        elen + elen
    }

    
    pub unsafe fn scale_expansion_zeroelim(&self,
        mut elen : i32, mut e : *const f64, mut b : f64, mut h : *mut f64
    ) -> i32 {
        let mut Q : f64;
        let mut sum : f64;
        let mut hh : f64;
        let mut product1 : f64;
        let mut product0 : f64;
        let mut eindex : i32;
        let mut hindex : i32;
        let mut enow : f64;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        c = self.splitter * b;
        abig = c - b;
        bhi = c - abig;
        blo = b - bhi;
        Q = *e.offset(0isize) * b;
        c = self.splitter * *e.offset(0isize);
        abig = c - *e.offset(0isize);
        ahi = c - abig;
        alo = *e.offset(0isize) - ahi;
        err1 = Q - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        hh = alo * blo - err3;
        hindex = 0i32;
        if hh != 0i32 as (f64) {
            *h.offset(
                 {
                     let _old = hindex;
                     hindex = hindex + 1;
                     _old
                 } as (isize)
             ) = hh;
        }
        eindex = 1i32;
        'loop3: loop {
            if !(eindex < elen) {
                break;
            }
            enow = *e.offset(eindex as (isize));
            product1 = enow * b;
            c = self.splitter * enow;
            abig = c - enow;
            ahi = c - abig;
            alo = enow - ahi;
            err1 = product1 - ahi * bhi;
            err2 = err1 - alo * bhi;
            err3 = err2 - ahi * blo;
            product0 = alo * blo - err3;
            sum = Q + product0;
            bvirt = sum - Q;
            avirt = sum - bvirt;
            bround = product0 - bvirt;
            around = Q - avirt;
            hh = around + bround;
            if hh != 0i32 as (f64) {
                *h.offset(
                     {
                         let _old = hindex;
                         hindex = hindex + 1;
                         _old
                     } as (isize)
                 ) = hh;
            }
            Q = product1 + sum;
            bvirt = Q - product1;
            hh = sum - bvirt;
            if hh != 0i32 as (f64) {
                *h.offset(
                     {
                         let _old = hindex;
                         hindex = hindex + 1;
                         _old
                     } as (isize)
                 ) = hh;
            }
            eindex = eindex + 1;
        }
        if Q != 0.0f64 || hindex == 0i32 {
            *h.offset(
                 {
                     let _old = hindex;
                     hindex = hindex + 1;
                     _old
                 } as (isize)
             ) = Q;
        }
        hindex
    }

    
    pub unsafe fn orient2dexact(&self,
        mut pa : *const f64, mut pb : *const f64, mut pc : *const f64
    ) -> f64 {
        let mut axby1 : f64;
        let mut axcy1 : f64;
        let mut bxcy1 : f64;
        let mut bxay1 : f64;
        let mut cxay1 : f64;
        let mut cxby1 : f64;
        let mut axby0 : f64;
        let mut axcy0 : f64;
        let mut bxcy0 : f64;
        let mut bxay0 : f64;
        let mut cxay0 : f64;
        let mut cxby0 : f64;
        let mut aterms : [f64; 4];
        let mut bterms : [f64; 4];
        let mut cterms : [f64; 4];
        let mut aterms3 : f64;
        let mut bterms3 : f64;
        let mut cterms3 : f64;
        let mut v : [f64; 8];
        let mut w : [f64; 12];
        let mut vlength : i32;
        let mut wlength : i32;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _0 : f64;

        axby1 = uninitialized();
        axcy1 = uninitialized();
        bxcy1 = uninitialized();
        bxay1 = uninitialized();
        cxay1 = uninitialized();
        cxby1 = uninitialized();
        axby0 = uninitialized();
        axcy0 = uninitialized();
        bxcy0 = uninitialized();
        bxay0 = uninitialized();
        cxay0 = uninitialized();
        cxby0 = uninitialized();
        aterms = uninitialized();
        bterms = uninitialized();
        cterms = uninitialized();
        aterms3 = uninitialized();
        bterms3 = uninitialized();
        cterms3 = uninitialized();
        v = uninitialized();
        w = uninitialized();
        vlength = uninitialized();
        wlength = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        ahi = uninitialized();
        alo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _0 = uninitialized();
        axby1 = *pa.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = axby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axby0 = alo * blo - err3;
        axcy1 = *pa.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = axcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axcy0 = alo * blo - err3;
        _i = axby0 - axcy0;
        bvirt = axby0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axcy0;
        around = axby0 - avirt;
        aterms[0usize] = around + bround;
        _j = axby1 + _i;
        bvirt = _j - axby1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = axby1 - avirt;
        _0 = around + bround;
        _i = _0 - axcy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axcy1;
        around = _0 - avirt;
        aterms[1usize] = around + bround;
        aterms3 = _j + _i;
        bvirt = aterms3 - _j;
        avirt = aterms3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        aterms[2usize] = around + bround;
        aterms[3usize] = aterms3;
        bxcy1 = *pb.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = bxcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxcy0 = alo * blo - err3;
        bxay1 = *pb.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = bxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxay0 = alo * blo - err3;
        _i = bxcy0 - bxay0;
        bvirt = bxcy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxay0;
        around = bxcy0 - avirt;
        bterms[0usize] = around + bround;
        _j = bxcy1 + _i;
        bvirt = _j - bxcy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bxcy1 - avirt;
        _0 = around + bround;
        _i = _0 - bxay1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxay1;
        around = _0 - avirt;
        bterms[1usize] = around + bround;
        bterms3 = _j + _i;
        bvirt = bterms3 - _j;
        avirt = bterms3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bterms[2usize] = around + bround;
        bterms[3usize] = bterms3;
        cxay1 = *pc.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = cxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxay0 = alo * blo - err3;
        cxby1 = *pc.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = cxby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxby0 = alo * blo - err3;
        _i = cxay0 - cxby0;
        bvirt = cxay0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxby0;
        around = cxay0 - avirt;
        cterms[0usize] = around + bround;
        _j = cxay1 + _i;
        bvirt = _j - cxay1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = cxay1 - avirt;
        _0 = around + bround;
        _i = _0 - cxby1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxby1;
        around = _0 - avirt;
        cterms[1usize] = around + bround;
        cterms3 = _j + _i;
        bvirt = cterms3 - _j;
        avirt = cterms3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        cterms[2usize] = around + bround;
        cterms[3usize] = cterms3;
        vlength = fast_expansion_sum_zeroelim(
                      4i32,
                      aterms.as_ptr(),
                      4i32,
                      bterms.as_ptr(),
                      v.as_mut_ptr()
                  );
        wlength = fast_expansion_sum_zeroelim(
                      vlength,
                      v.as_ptr(),
                      4i32,
                      cterms.as_ptr(),
                      w.as_mut_ptr()
                  );
        w[(wlength - 1i32) as (usize)]
    }

    
    pub unsafe fn orient2dslow(&self,
        mut pa : *const f64, mut pb : *const f64, mut pc : *const f64
    ) -> f64 {
        let mut acx : f64;
        let mut acy : f64;
        let mut bcx : f64;
        let mut bcy : f64;
        let mut acxtail : f64;
        let mut acytail : f64;
        let mut bcxtail : f64;
        let mut bcytail : f64;
        let mut negate : f64;
        let mut negatetail : f64;
        let mut axby : [f64; 8];
        let mut bxay : [f64; 8];
        let mut axby7 : f64;
        let mut bxay7 : f64;
        let mut deter : [f64; 16];
        let mut deterlen : i32;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut a0hi : f64;
        let mut a0lo : f64;
        let mut a1hi : f64;
        let mut a1lo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _k : f64;
        let mut _l : f64;
        let mut _m : f64;
        let mut _n : f64;
        let mut _0 : f64;
        let mut _1 : f64;
        let mut _2 : f64;

        acx = uninitialized();
        acy = uninitialized();
        bcx = uninitialized();
        bcy = uninitialized();
        acxtail = uninitialized();
        acytail = uninitialized();
        bcxtail = uninitialized();
        bcytail = uninitialized();
        negate = uninitialized();
        negatetail = uninitialized();
        axby = uninitialized();
        bxay = uninitialized();
        axby7 = uninitialized();
        bxay7 = uninitialized();
        deter = uninitialized();
        deterlen = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        a0hi = uninitialized();
        a0lo = uninitialized();
        a1hi = uninitialized();
        a1lo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _k = uninitialized();
        _l = uninitialized();
        _m = uninitialized();
        _n = uninitialized();
        _0 = uninitialized();
        _1 = uninitialized();
        _2 = uninitialized();
        acx = *pa.offset(0isize) - *pc.offset(0isize);
        bvirt = *pa.offset(0isize) - acx;
        avirt = acx + bvirt;
        bround = bvirt - *pc.offset(0isize);
        around = *pa.offset(0isize) - avirt;
        acxtail = around + bround;
        acy = *pa.offset(1isize) - *pc.offset(1isize);
        bvirt = *pa.offset(1isize) - acy;
        avirt = acy + bvirt;
        bround = bvirt - *pc.offset(1isize);
        around = *pa.offset(1isize) - avirt;
        acytail = around + bround;
        bcx = *pb.offset(0isize) - *pc.offset(0isize);
        bvirt = *pb.offset(0isize) - bcx;
        avirt = bcx + bvirt;
        bround = bvirt - *pc.offset(0isize);
        around = *pb.offset(0isize) - avirt;
        bcxtail = around + bround;
        bcy = *pb.offset(1isize) - *pc.offset(1isize);
        bvirt = *pb.offset(1isize) - bcy;
        avirt = bcy + bvirt;
        bround = bvirt - *pc.offset(1isize);
        around = *pb.offset(1isize) - avirt;
        bcytail = around + bround;
        c = self.splitter * acxtail;
        abig = c - acxtail;
        a0hi = c - abig;
        a0lo = acxtail - a0hi;
        c = self.splitter * bcytail;
        abig = c - bcytail;
        bhi = c - abig;
        blo = bcytail - bhi;
        _i = acxtail * bcytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        axby[0usize] = a0lo * blo - err3;
        c = self.splitter * acx;
        abig = c - acx;
        a1hi = c - abig;
        a1lo = acx - a1hi;
        _j = acx * bcytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * bcy;
        abig = c - bcy;
        bhi = c - abig;
        blo = bcy - bhi;
        _i = acxtail * bcy;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = acx * bcy;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        axby[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        axby[5usize] = around + bround;
        axby7 = _m + _k;
        bvirt = axby7 - _m;
        avirt = axby7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        axby[6usize] = around + bround;
        axby[7usize] = axby7;
        negate = -acy;
        negatetail = -acytail;
        c = self.splitter * bcxtail;
        abig = c - bcxtail;
        a0hi = c - abig;
        a0lo = bcxtail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = bcxtail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        bxay[0usize] = a0lo * blo - err3;
        c = self.splitter * bcx;
        abig = c - bcx;
        a1hi = c - abig;
        a1lo = bcx - a1hi;
        _j = bcx * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = bcxtail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = bcx * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        bxay[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        bxay[5usize] = around + bround;
        bxay7 = _m + _k;
        bvirt = bxay7 - _m;
        avirt = bxay7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        bxay[6usize] = around + bround;
        bxay[7usize] = bxay7;
        deterlen = fast_expansion_sum_zeroelim(
                       8i32,
                       axby.as_ptr(),
                       8i32,
                       bxay.as_ptr(),
                       deter.as_mut_ptr()
                   );
        deter[(deterlen - 1i32) as (usize)]
    }

    
    pub unsafe fn orient2dadapt(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut detsum : f64
    ) -> f64 {
        let mut acx : f64;
        let mut acy : f64;
        let mut bcx : f64;
        let mut bcy : f64;
        let mut acxtail : f64;
        let mut acytail : f64;
        let mut bcxtail : f64;
        let mut bcytail : f64;
        let mut detleft : f64;
        let mut detright : f64;
        let mut detlefttail : f64;
        let mut detrighttail : f64;
        let mut det : f64;
        let mut errbound : f64;
        let mut B : [f64; 4];
        let mut C1 : [f64; 8];
        let mut C2 : [f64; 12];
        let mut D : [f64; 16];
        let mut B3 : f64;
        let mut C1length : i32;
        let mut C2length : i32;
        let mut Dlength : i32;
        let mut u : [f64; 4];
        let mut u3 : f64;
        let mut s1 : f64;
        let mut t1 : f64;
        let mut s0 : f64;
        let mut t0 : f64;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _0 : f64;

        acx = uninitialized();
        acy = uninitialized();
        bcx = uninitialized();
        bcy = uninitialized();
        acxtail = uninitialized();
        acytail = uninitialized();
        bcxtail = uninitialized();
        bcytail = uninitialized();
        detleft = uninitialized();
        detright = uninitialized();
        detlefttail = uninitialized();
        detrighttail = uninitialized();
        det = uninitialized();
        errbound = uninitialized();
        B = uninitialized();
        C1 = uninitialized();
        C2 = uninitialized();
        D = uninitialized();
        B3 = uninitialized();
        C1length = uninitialized();
        C2length = uninitialized();
        Dlength = uninitialized();
        u = uninitialized();
        u3 = uninitialized();
        s1 = uninitialized();
        t1 = uninitialized();
        s0 = uninitialized();
        t0 = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        ahi = uninitialized();
        alo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _0 = uninitialized();
        acx = *pa.offset(0isize) - *pc.offset(0isize);
        bcx = *pb.offset(0isize) - *pc.offset(0isize);
        acy = *pa.offset(1isize) - *pc.offset(1isize);
        bcy = *pb.offset(1isize) - *pc.offset(1isize);
        detleft = acx * bcy;
        c = self.splitter * acx;
        abig = c - acx;
        ahi = c - abig;
        alo = acx - ahi;
        c = self.splitter * bcy;
        abig = c - bcy;
        bhi = c - abig;
        blo = bcy - bhi;
        err1 = detleft - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        detlefttail = alo * blo - err3;
        detright = acy * bcx;
        c = self.splitter * acy;
        abig = c - acy;
        ahi = c - abig;
        alo = acy - ahi;
        c = self.splitter * bcx;
        abig = c - bcx;
        bhi = c - abig;
        blo = bcx - bhi;
        err1 = detright - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        detrighttail = alo * blo - err3;
        _i = detlefttail - detrighttail;
        bvirt = detlefttail - _i;
        avirt = _i + bvirt;
        bround = bvirt - detrighttail;
        around = detlefttail - avirt;
        B[0usize] = around + bround;
        _j = detleft + _i;
        bvirt = _j - detleft;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = detleft - avirt;
        _0 = around + bround;
        _i = _0 - detright;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - detright;
        around = _0 - avirt;
        B[1usize] = around + bround;
        B3 = _j + _i;
        bvirt = B3 - _j;
        avirt = B3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        B[2usize] = around + bround;
        B[3usize] = B3;
        det = estimate(4i32,B.as_ptr());
        errbound = self.ccwerrboundB * detsum;
        if det >= errbound || -det >= errbound {
            det
        } else {
            bvirt = *pa.offset(0isize) - acx;
            avirt = acx + bvirt;
            bround = bvirt - *pc.offset(0isize);
            around = *pa.offset(0isize) - avirt;
            acxtail = around + bround;
            bvirt = *pb.offset(0isize) - bcx;
            avirt = bcx + bvirt;
            bround = bvirt - *pc.offset(0isize);
            around = *pb.offset(0isize) - avirt;
            bcxtail = around + bround;
            bvirt = *pa.offset(1isize) - acy;
            avirt = acy + bvirt;
            bround = bvirt - *pc.offset(1isize);
            around = *pa.offset(1isize) - avirt;
            acytail = around + bround;
            bvirt = *pb.offset(1isize) - bcy;
            avirt = bcy + bvirt;
            bround = bvirt - *pc.offset(1isize);
            around = *pb.offset(1isize) - avirt;
            bcytail = around + bround;
            (if acxtail == 0.0f64 && (acytail == 0.0f64) && (bcxtail == 0.0f64) && (bcytail == 0.0f64) {
                 det
             } else {
                 errbound = self.ccwerrboundC * detsum + self.resulterrbound * Absolute(det);
                 det = det + (acx * bcytail + bcy * acxtail - (acy * bcxtail + bcx * acytail));
                 (if det >= errbound || -det >= errbound {
                      det
                  } else {
                      s1 = acxtail * bcy;
                      c = self.splitter * acxtail;
                      abig = c - acxtail;
                      ahi = c - abig;
                      alo = acxtail - ahi;
                      c = self.splitter * bcy;
                      abig = c - bcy;
                      bhi = c - abig;
                      blo = bcy - bhi;
                      err1 = s1 - ahi * bhi;
                      err2 = err1 - alo * bhi;
                      err3 = err2 - ahi * blo;
                      s0 = alo * blo - err3;
                      t1 = acytail * bcx;
                      c = self.splitter * acytail;
                      abig = c - acytail;
                      ahi = c - abig;
                      alo = acytail - ahi;
                      c = self.splitter * bcx;
                      abig = c - bcx;
                      bhi = c - abig;
                      blo = bcx - bhi;
                      err1 = t1 - ahi * bhi;
                      err2 = err1 - alo * bhi;
                      err3 = err2 - ahi * blo;
                      t0 = alo * blo - err3;
                      _i = s0 - t0;
                      bvirt = s0 - _i;
                      avirt = _i + bvirt;
                      bround = bvirt - t0;
                      around = s0 - avirt;
                      u[0usize] = around + bround;
                      _j = s1 + _i;
                      bvirt = _j - s1;
                      avirt = _j - bvirt;
                      bround = _i - bvirt;
                      around = s1 - avirt;
                      _0 = around + bround;
                      _i = _0 - t1;
                      bvirt = _0 - _i;
                      avirt = _i + bvirt;
                      bround = bvirt - t1;
                      around = _0 - avirt;
                      u[1usize] = around + bround;
                      u3 = _j + _i;
                      bvirt = u3 - _j;
                      avirt = u3 - bvirt;
                      bround = _i - bvirt;
                      around = _j - avirt;
                      u[2usize] = around + bround;
                      u[3usize] = u3;
                      C1length = fast_expansion_sum_zeroelim(
                                     4i32,
                                     B.as_ptr(),
                                     4i32,
                                     u.as_ptr(),
                                     C1.as_mut_ptr()
                                 );
                      s1 = acx * bcytail;
                      c = self.splitter * acx;
                      abig = c - acx;
                      ahi = c - abig;
                      alo = acx - ahi;
                      c = self.splitter * bcytail;
                      abig = c - bcytail;
                      bhi = c - abig;
                      blo = bcytail - bhi;
                      err1 = s1 - ahi * bhi;
                      err2 = err1 - alo * bhi;
                      err3 = err2 - ahi * blo;
                      s0 = alo * blo - err3;
                      t1 = acy * bcxtail;
                      c = self.splitter * acy;
                      abig = c - acy;
                      ahi = c - abig;
                      alo = acy - ahi;
                      c = self.splitter * bcxtail;
                      abig = c - bcxtail;
                      bhi = c - abig;
                      blo = bcxtail - bhi;
                      err1 = t1 - ahi * bhi;
                      err2 = err1 - alo * bhi;
                      err3 = err2 - ahi * blo;
                      t0 = alo * blo - err3;
                      _i = s0 - t0;
                      bvirt = s0 - _i;
                      avirt = _i + bvirt;
                      bround = bvirt - t0;
                      around = s0 - avirt;
                      u[0usize] = around + bround;
                      _j = s1 + _i;
                      bvirt = _j - s1;
                      avirt = _j - bvirt;
                      bround = _i - bvirt;
                      around = s1 - avirt;
                      _0 = around + bround;
                      _i = _0 - t1;
                      bvirt = _0 - _i;
                      avirt = _i + bvirt;
                      bround = bvirt - t1;
                      around = _0 - avirt;
                      u[1usize] = around + bround;
                      u3 = _j + _i;
                      bvirt = u3 - _j;
                      avirt = u3 - bvirt;
                      bround = _i - bvirt;
                      around = _j - avirt;
                      u[2usize] = around + bround;
                      u[3usize] = u3;
                      C2length = fast_expansion_sum_zeroelim(
                                     C1length,
                                     C1.as_ptr(),
                                     4i32,
                                     u.as_ptr(),
                                     C2.as_mut_ptr()
                                 );
                      s1 = acxtail * bcytail;
                      c = self.splitter * acxtail;
                      abig = c - acxtail;
                      ahi = c - abig;
                      alo = acxtail - ahi;
                      c = self.splitter * bcytail;
                      abig = c - bcytail;
                      bhi = c - abig;
                      blo = bcytail - bhi;
                      err1 = s1 - ahi * bhi;
                      err2 = err1 - alo * bhi;
                      err3 = err2 - ahi * blo;
                      s0 = alo * blo - err3;
                      t1 = acytail * bcxtail;
                      c = self.splitter * acytail;
                      abig = c - acytail;
                      ahi = c - abig;
                      alo = acytail - ahi;
                      c = self.splitter * bcxtail;
                      abig = c - bcxtail;
                      bhi = c - abig;
                      blo = bcxtail - bhi;
                      err1 = t1 - ahi * bhi;
                      err2 = err1 - alo * bhi;
                      err3 = err2 - ahi * blo;
                      t0 = alo * blo - err3;
                      _i = s0 - t0;
                      bvirt = s0 - _i;
                      avirt = _i + bvirt;
                      bround = bvirt - t0;
                      around = s0 - avirt;
                      u[0usize] = around + bround;
                      _j = s1 + _i;
                      bvirt = _j - s1;
                      avirt = _j - bvirt;
                      bround = _i - bvirt;
                      around = s1 - avirt;
                      _0 = around + bround;
                      _i = _0 - t1;
                      bvirt = _0 - _i;
                      avirt = _i + bvirt;
                      bround = bvirt - t1;
                      around = _0 - avirt;
                      u[1usize] = around + bround;
                      u3 = _j + _i;
                      bvirt = u3 - _j;
                      avirt = u3 - bvirt;
                      bround = _i - bvirt;
                      around = _j - avirt;
                      u[2usize] = around + bround;
                      u[3usize] = u3;
                      Dlength = fast_expansion_sum_zeroelim(
                                    C2length,
                                    C2.as_ptr(),
                                    4i32,
                                    u.as_ptr(),
                                    D.as_mut_ptr()
                                );
                      D[(Dlength - 1i32) as (usize)]
                  })
             })
        }
    }

    
    pub unsafe fn orient2d(&self,
        mut pa : *const f64, mut pb : *const f64, mut pc : *const f64
    ) -> f64 {
        let mut detleft : f64;
        let mut detright : f64;
        let mut det : f64;
        let mut detsum : f64;
        let mut errbound : f64;
        detleft = (*pa.offset(0isize) - *pc.offset(0isize)) * (*pb.offset(
                                                                    1isize
                                                                ) - *pc.offset(1isize));
        detright = (*pa.offset(1isize) - *pc.offset(1isize)) * (*pb.offset(
                                                                     0isize
                                                                 ) - *pc.offset(0isize));
        det = detleft - detright;
        if detleft > 0.0f64 {
            if detright <= 0.0f64 {
                return det;
            } else {
                detsum = detleft + detright;
            }
        } else if detleft < 0.0f64 {
            if detright >= 0.0f64 {
                return det;
            } else {
                detsum = -detleft - detright;
            }
        } else {
            return det;
        }
        errbound = self.ccwerrboundA * detsum;
        if det >= errbound || -det >= errbound {
            det
        } else {
            self.orient2dadapt(pa,pb,pc,detsum)
        }
    }

    
    pub unsafe fn orient3dexact(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64
    ) -> f64 {
        let mut axby1 : f64;
        let mut bxcy1 : f64;
        let mut cxdy1 : f64;
        let mut dxay1 : f64;
        let mut axcy1 : f64;
        let mut bxdy1 : f64;
        let mut bxay1 : f64;
        let mut cxby1 : f64;
        let mut dxcy1 : f64;
        let mut axdy1 : f64;
        let mut cxay1 : f64;
        let mut dxby1 : f64;
        let mut axby0 : f64;
        let mut bxcy0 : f64;
        let mut cxdy0 : f64;
        let mut dxay0 : f64;
        let mut axcy0 : f64;
        let mut bxdy0 : f64;
        let mut bxay0 : f64;
        let mut cxby0 : f64;
        let mut dxcy0 : f64;
        let mut axdy0 : f64;
        let mut cxay0 : f64;
        let mut dxby0 : f64;
        let mut ab : [f64; 4];
        let mut bc : [f64; 4];
        let mut cd : [f64; 4];
        let mut da : [f64; 4];
        let mut ac : [f64; 4];
        let mut bd : [f64; 4];
        let mut temp8 : [f64; 8];
        let mut templen : i32;
        let mut abc : [f64; 12];
        let mut bcd : [f64; 12];
        let mut cda : [f64; 12];
        let mut dab : [f64; 12];
        let mut abclen : i32;
        let mut bcdlen : i32;
        let mut cdalen : i32;
        let mut dablen : i32;
        let mut adet : [f64; 24];
        let mut bdet : [f64; 24];
        let mut cdet : [f64; 24];
        let mut ddet : [f64; 24];
        let mut alen : i32;
        let mut blen : i32;
        let mut clen : i32;
        let mut dlen : i32;
        let mut abdet : [f64; 48];
        let mut cddet : [f64; 48];
        let mut ablen : i32;
        let mut cdlen : i32;
        let mut deter : [f64; 96];
        let mut deterlen : i32;
        let mut i : i32;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _0 : f64;

        axby1 = uninitialized();
        bxcy1 = uninitialized();
        cxdy1 = uninitialized();
        dxay1 = uninitialized();
        axcy1 = uninitialized();
        bxdy1 = uninitialized();
        bxay1 = uninitialized();
        cxby1 = uninitialized();
        dxcy1 = uninitialized();
        axdy1 = uninitialized();
        cxay1 = uninitialized();
        dxby1 = uninitialized();
        axby0 = uninitialized();
        bxcy0 = uninitialized();
        cxdy0 = uninitialized();
        dxay0 = uninitialized();
        axcy0 = uninitialized();
        bxdy0 = uninitialized();
        bxay0 = uninitialized();
        cxby0 = uninitialized();
        dxcy0 = uninitialized();
        axdy0 = uninitialized();
        cxay0 = uninitialized();
        dxby0 = uninitialized();
        ab = uninitialized();
        bc = uninitialized();
        cd = uninitialized();
        da = uninitialized();
        ac = uninitialized();
        bd = uninitialized();
        temp8 = uninitialized();
        templen = uninitialized();
        abc = uninitialized();
        bcd = uninitialized();
        cda = uninitialized();
        dab = uninitialized();
        abclen = uninitialized();
        bcdlen = uninitialized();
        cdalen = uninitialized();
        dablen = uninitialized();
        adet = uninitialized();
        bdet = uninitialized();
        cdet = uninitialized();
        ddet = uninitialized();
        alen = uninitialized();
        blen = uninitialized();
        clen = uninitialized();
        dlen = uninitialized();
        abdet = uninitialized();
        cddet = uninitialized();
        ablen = uninitialized();
        cdlen = uninitialized();
        deter = uninitialized();
        deterlen = uninitialized();
        i = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        ahi = uninitialized();
        alo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _0 = uninitialized();
        axby1 = *pa.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = axby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axby0 = alo * blo - err3;
        bxay1 = *pb.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = bxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxay0 = alo * blo - err3;
        _i = axby0 - bxay0;
        bvirt = axby0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxay0;
        around = axby0 - avirt;
        ab[0usize] = around + bround;
        _j = axby1 + _i;
        bvirt = _j - axby1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = axby1 - avirt;
        _0 = around + bround;
        _i = _0 - bxay1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxay1;
        around = _0 - avirt;
        ab[1usize] = around + bround;
        ab[3usize] = _j + _i;
        bvirt = ab[3usize] - _j;
        avirt = ab[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ab[2usize] = around + bround;
        bxcy1 = *pb.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = bxcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxcy0 = alo * blo - err3;
        cxby1 = *pc.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = cxby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxby0 = alo * blo - err3;
        _i = bxcy0 - cxby0;
        bvirt = bxcy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxby0;
        around = bxcy0 - avirt;
        bc[0usize] = around + bround;
        _j = bxcy1 + _i;
        bvirt = _j - bxcy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bxcy1 - avirt;
        _0 = around + bround;
        _i = _0 - cxby1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxby1;
        around = _0 - avirt;
        bc[1usize] = around + bround;
        bc[3usize] = _j + _i;
        bvirt = bc[3usize] - _j;
        avirt = bc[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bc[2usize] = around + bround;
        cxdy1 = *pc.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = cxdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxdy0 = alo * blo - err3;
        dxcy1 = *pd.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = dxcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxcy0 = alo * blo - err3;
        _i = cxdy0 - dxcy0;
        bvirt = cxdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxcy0;
        around = cxdy0 - avirt;
        cd[0usize] = around + bround;
        _j = cxdy1 + _i;
        bvirt = _j - cxdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = cxdy1 - avirt;
        _0 = around + bround;
        _i = _0 - dxcy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxcy1;
        around = _0 - avirt;
        cd[1usize] = around + bround;
        cd[3usize] = _j + _i;
        bvirt = cd[3usize] - _j;
        avirt = cd[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        cd[2usize] = around + bround;
        dxay1 = *pd.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = dxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxay0 = alo * blo - err3;
        axdy1 = *pa.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = axdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axdy0 = alo * blo - err3;
        _i = dxay0 - axdy0;
        bvirt = dxay0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axdy0;
        around = dxay0 - avirt;
        da[0usize] = around + bround;
        _j = dxay1 + _i;
        bvirt = _j - dxay1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = dxay1 - avirt;
        _0 = around + bround;
        _i = _0 - axdy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axdy1;
        around = _0 - avirt;
        da[1usize] = around + bround;
        da[3usize] = _j + _i;
        bvirt = da[3usize] - _j;
        avirt = da[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        da[2usize] = around + bround;
        axcy1 = *pa.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = axcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axcy0 = alo * blo - err3;
        cxay1 = *pc.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = cxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxay0 = alo * blo - err3;
        _i = axcy0 - cxay0;
        bvirt = axcy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxay0;
        around = axcy0 - avirt;
        ac[0usize] = around + bround;
        _j = axcy1 + _i;
        bvirt = _j - axcy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = axcy1 - avirt;
        _0 = around + bround;
        _i = _0 - cxay1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxay1;
        around = _0 - avirt;
        ac[1usize] = around + bround;
        ac[3usize] = _j + _i;
        bvirt = ac[3usize] - _j;
        avirt = ac[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ac[2usize] = around + bround;
        bxdy1 = *pb.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = bxdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxdy0 = alo * blo - err3;
        dxby1 = *pd.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = dxby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxby0 = alo * blo - err3;
        _i = bxdy0 - dxby0;
        bvirt = bxdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxby0;
        around = bxdy0 - avirt;
        bd[0usize] = around + bround;
        _j = bxdy1 + _i;
        bvirt = _j - bxdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bxdy1 - avirt;
        _0 = around + bround;
        _i = _0 - dxby1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxby1;
        around = _0 - avirt;
        bd[1usize] = around + bround;
        bd[3usize] = _j + _i;
        bvirt = bd[3usize] - _j;
        avirt = bd[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bd[2usize] = around + bround;
        templen = fast_expansion_sum_zeroelim(
                      4i32,
                      cd.as_ptr(),
                      4i32,
                      da.as_ptr(),
                      temp8.as_mut_ptr()
                  );
        cdalen = fast_expansion_sum_zeroelim(
                     templen,
                     temp8.as_ptr(),
                     4i32,
                     ac.as_ptr(),
                     cda.as_mut_ptr()
                 );
        templen = fast_expansion_sum_zeroelim(
                      4i32,
                      da.as_ptr(),
                      4i32,
                      ab.as_ptr(),
                      temp8.as_mut_ptr()
                  );
        dablen = fast_expansion_sum_zeroelim(
                     templen,
                     temp8.as_ptr(),
                     4i32,
                     bd.as_ptr(),
                     dab.as_mut_ptr()
                 );
        i = 0i32;
        'loop1: loop {
            if !(i < 4i32) {
                break;
            }
            bd[i as (usize)] = -bd[i as (usize)];
            ac[i as (usize)] = -ac[i as (usize)];
            i = i + 1;
        }
        templen = fast_expansion_sum_zeroelim(
                      4i32,
                      ab.as_ptr(),
                      4i32,
                      bc.as_ptr(),
                      temp8.as_mut_ptr()
                  );
        abclen = fast_expansion_sum_zeroelim(
                     templen,
                     temp8.as_ptr(),
                     4i32,
                     ac.as_ptr(),
                     abc.as_mut_ptr()
                 );
        templen = fast_expansion_sum_zeroelim(
                      4i32,
                      bc.as_ptr(),
                      4i32,
                      cd.as_ptr(),
                      temp8.as_mut_ptr()
                  );
        bcdlen = fast_expansion_sum_zeroelim(
                     templen,
                     temp8.as_ptr(),
                     4i32,
                     bd.as_ptr(),
                     bcd.as_mut_ptr()
                 );
        alen = self.scale_expansion_zeroelim(
                   bcdlen,
                   bcd.as_ptr(),
                   *pa.offset(2isize),
                   adet.as_mut_ptr()
               );
        blen = self.scale_expansion_zeroelim(
                   cdalen,
                   cda.as_ptr(),
                   -*pb.offset(2isize),
                   bdet.as_mut_ptr()
               );
        clen = self.scale_expansion_zeroelim(
                   dablen,
                   dab.as_ptr(),
                   *pc.offset(2isize),
                   cdet.as_mut_ptr()
               );
        dlen = self.scale_expansion_zeroelim(
                   abclen,
                   abc.as_ptr(),
                   -*pd.offset(2isize),
                   ddet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        cdlen = fast_expansion_sum_zeroelim(
                    clen,
                    cdet.as_ptr(),
                    dlen,
                    ddet.as_ptr(),
                    cddet.as_mut_ptr()
                );
        deterlen = fast_expansion_sum_zeroelim(
                       ablen,
                       abdet.as_ptr(),
                       cdlen,
                       cddet.as_ptr(),
                       deter.as_mut_ptr()
                   );
        deter[(deterlen - 1i32) as (usize)]
    }

    
    pub unsafe fn orient3dslow(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64
    ) -> f64 {
        let mut adx : f64;
        let mut ady : f64;
        let mut adz : f64;
        let mut bdx : f64;
        let mut bdy : f64;
        let mut bdz : f64;
        let mut cdx : f64;
        let mut cdy : f64;
        let mut cdz : f64;
        let mut adxtail : f64;
        let mut adytail : f64;
        let mut adztail : f64;
        let mut bdxtail : f64;
        let mut bdytail : f64;
        let mut bdztail : f64;
        let mut cdxtail : f64;
        let mut cdytail : f64;
        let mut cdztail : f64;
        let mut negate : f64;
        let mut negatetail : f64;
        let mut axby7 : f64;
        let mut bxcy7 : f64;
        let mut axcy7 : f64;
        let mut bxay7 : f64;
        let mut cxby7 : f64;
        let mut cxay7 : f64;
        let mut axby : [f64; 8];
        let mut bxcy : [f64; 8];
        let mut axcy : [f64; 8];
        let mut bxay : [f64; 8];
        let mut cxby : [f64; 8];
        let mut cxay : [f64; 8];
        let mut temp16 : [f64; 16];
        let mut temp32 : [f64; 32];
        let mut temp32t : [f64; 32];
        let mut temp16len : i32;
        let mut temp32len : i32;
        let mut temp32tlen : i32;
        let mut adet : [f64; 64];
        let mut bdet : [f64; 64];
        let mut cdet : [f64; 64];
        let mut alen : i32;
        let mut blen : i32;
        let mut clen : i32;
        let mut abdet : [f64; 128];
        let mut ablen : i32;
        let mut deter : [f64; 192];
        let mut deterlen : i32;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut a0hi : f64;
        let mut a0lo : f64;
        let mut a1hi : f64;
        let mut a1lo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _k : f64;
        let mut _l : f64;
        let mut _m : f64;
        let mut _n : f64;
        let mut _0 : f64;
        let mut _1 : f64;
        let mut _2 : f64;

        adx = uninitialized();
        ady = uninitialized();
        adz = uninitialized();
        bdx = uninitialized();
        bdy = uninitialized();
        bdz = uninitialized();
        cdx = uninitialized();
        cdy = uninitialized();
        cdz = uninitialized();
        adxtail = uninitialized();
        adytail = uninitialized();
        adztail = uninitialized();
        bdxtail = uninitialized();
        bdytail = uninitialized();
        bdztail = uninitialized();
        cdxtail = uninitialized();
        cdytail = uninitialized();
        cdztail = uninitialized();
        negate = uninitialized();
        negatetail = uninitialized();
        axby7 = uninitialized();
        bxcy7 = uninitialized();
        axcy7 = uninitialized();
        bxay7 = uninitialized();
        cxby7 = uninitialized();
        cxay7 = uninitialized();
        axby = uninitialized();
        bxcy = uninitialized();
        axcy = uninitialized();
        bxay = uninitialized();
        cxby = uninitialized();
        cxay = uninitialized();
        temp16 = uninitialized();
        temp32 = uninitialized();
        temp32t = uninitialized();
        temp16len = uninitialized();
        temp32len = uninitialized();
        temp32tlen = uninitialized();
        adet = uninitialized();
        bdet = uninitialized();
        cdet = uninitialized();
        alen = uninitialized();
        blen = uninitialized();
        clen = uninitialized();
        abdet = uninitialized();
        ablen = uninitialized();
        deter = uninitialized();
        deterlen = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        a0hi = uninitialized();
        a0lo = uninitialized();
        a1hi = uninitialized();
        a1lo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _k = uninitialized();
        _l = uninitialized();
        _m = uninitialized();
        _n = uninitialized();
        _0 = uninitialized();
        _1 = uninitialized();
        _2 = uninitialized();
        adx = *pa.offset(0isize) - *pd.offset(0isize);
        bvirt = *pa.offset(0isize) - adx;
        avirt = adx + bvirt;
        bround = bvirt - *pd.offset(0isize);
        around = *pa.offset(0isize) - avirt;
        adxtail = around + bround;
        ady = *pa.offset(1isize) - *pd.offset(1isize);
        bvirt = *pa.offset(1isize) - ady;
        avirt = ady + bvirt;
        bround = bvirt - *pd.offset(1isize);
        around = *pa.offset(1isize) - avirt;
        adytail = around + bround;
        adz = *pa.offset(2isize) - *pd.offset(2isize);
        bvirt = *pa.offset(2isize) - adz;
        avirt = adz + bvirt;
        bround = bvirt - *pd.offset(2isize);
        around = *pa.offset(2isize) - avirt;
        adztail = around + bround;
        bdx = *pb.offset(0isize) - *pd.offset(0isize);
        bvirt = *pb.offset(0isize) - bdx;
        avirt = bdx + bvirt;
        bround = bvirt - *pd.offset(0isize);
        around = *pb.offset(0isize) - avirt;
        bdxtail = around + bround;
        bdy = *pb.offset(1isize) - *pd.offset(1isize);
        bvirt = *pb.offset(1isize) - bdy;
        avirt = bdy + bvirt;
        bround = bvirt - *pd.offset(1isize);
        around = *pb.offset(1isize) - avirt;
        bdytail = around + bround;
        bdz = *pb.offset(2isize) - *pd.offset(2isize);
        bvirt = *pb.offset(2isize) - bdz;
        avirt = bdz + bvirt;
        bround = bvirt - *pd.offset(2isize);
        around = *pb.offset(2isize) - avirt;
        bdztail = around + bround;
        cdx = *pc.offset(0isize) - *pd.offset(0isize);
        bvirt = *pc.offset(0isize) - cdx;
        avirt = cdx + bvirt;
        bround = bvirt - *pd.offset(0isize);
        around = *pc.offset(0isize) - avirt;
        cdxtail = around + bround;
        cdy = *pc.offset(1isize) - *pd.offset(1isize);
        bvirt = *pc.offset(1isize) - cdy;
        avirt = cdy + bvirt;
        bround = bvirt - *pd.offset(1isize);
        around = *pc.offset(1isize) - avirt;
        cdytail = around + bround;
        cdz = *pc.offset(2isize) - *pd.offset(2isize);
        bvirt = *pc.offset(2isize) - cdz;
        avirt = cdz + bvirt;
        bround = bvirt - *pd.offset(2isize);
        around = *pc.offset(2isize) - avirt;
        cdztail = around + bround;
        c = self.splitter * adxtail;
        abig = c - adxtail;
        a0hi = c - abig;
        a0lo = adxtail - a0hi;
        c = self.splitter * bdytail;
        abig = c - bdytail;
        bhi = c - abig;
        blo = bdytail - bhi;
        _i = adxtail * bdytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        axby[0usize] = a0lo * blo - err3;
        c = self.splitter * adx;
        abig = c - adx;
        a1hi = c - abig;
        a1lo = adx - a1hi;
        _j = adx * bdytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * bdy;
        abig = c - bdy;
        bhi = c - abig;
        blo = bdy - bhi;
        _i = adxtail * bdy;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = adx * bdy;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        axby[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        axby[5usize] = around + bround;
        axby7 = _m + _k;
        bvirt = axby7 - _m;
        avirt = axby7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        axby[6usize] = around + bround;
        axby[7usize] = axby7;
        negate = -ady;
        negatetail = -adytail;
        c = self.splitter * bdxtail;
        abig = c - bdxtail;
        a0hi = c - abig;
        a0lo = bdxtail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = bdxtail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        bxay[0usize] = a0lo * blo - err3;
        c = self.splitter * bdx;
        abig = c - bdx;
        a1hi = c - abig;
        a1lo = bdx - a1hi;
        _j = bdx * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = bdxtail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = bdx * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        bxay[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        bxay[5usize] = around + bround;
        bxay7 = _m + _k;
        bvirt = bxay7 - _m;
        avirt = bxay7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        bxay[6usize] = around + bround;
        bxay[7usize] = bxay7;
        c = self.splitter * bdxtail;
        abig = c - bdxtail;
        a0hi = c - abig;
        a0lo = bdxtail - a0hi;
        c = self.splitter * cdytail;
        abig = c - cdytail;
        bhi = c - abig;
        blo = cdytail - bhi;
        _i = bdxtail * cdytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        bxcy[0usize] = a0lo * blo - err3;
        c = self.splitter * bdx;
        abig = c - bdx;
        a1hi = c - abig;
        a1lo = bdx - a1hi;
        _j = bdx * cdytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * cdy;
        abig = c - cdy;
        bhi = c - abig;
        blo = cdy - bhi;
        _i = bdxtail * cdy;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = bdx * cdy;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        bxcy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        bxcy[5usize] = around + bround;
        bxcy7 = _m + _k;
        bvirt = bxcy7 - _m;
        avirt = bxcy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        bxcy[6usize] = around + bround;
        bxcy[7usize] = bxcy7;
        negate = -bdy;
        negatetail = -bdytail;
        c = self.splitter * cdxtail;
        abig = c - cdxtail;
        a0hi = c - abig;
        a0lo = cdxtail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = cdxtail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        cxby[0usize] = a0lo * blo - err3;
        c = self.splitter * cdx;
        abig = c - cdx;
        a1hi = c - abig;
        a1lo = cdx - a1hi;
        _j = cdx * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = cdxtail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = cdx * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        cxby[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        cxby[5usize] = around + bround;
        cxby7 = _m + _k;
        bvirt = cxby7 - _m;
        avirt = cxby7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        cxby[6usize] = around + bround;
        cxby[7usize] = cxby7;
        c = self.splitter * cdxtail;
        abig = c - cdxtail;
        a0hi = c - abig;
        a0lo = cdxtail - a0hi;
        c = self.splitter * adytail;
        abig = c - adytail;
        bhi = c - abig;
        blo = adytail - bhi;
        _i = cdxtail * adytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        cxay[0usize] = a0lo * blo - err3;
        c = self.splitter * cdx;
        abig = c - cdx;
        a1hi = c - abig;
        a1lo = cdx - a1hi;
        _j = cdx * adytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * ady;
        abig = c - ady;
        bhi = c - abig;
        blo = ady - bhi;
        _i = cdxtail * ady;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = cdx * ady;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        cxay[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        cxay[5usize] = around + bround;
        cxay7 = _m + _k;
        bvirt = cxay7 - _m;
        avirt = cxay7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        cxay[6usize] = around + bround;
        cxay[7usize] = cxay7;
        negate = -cdy;
        negatetail = -cdytail;
        c = self.splitter * adxtail;
        abig = c - adxtail;
        a0hi = c - abig;
        a0lo = adxtail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = adxtail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        axcy[0usize] = a0lo * blo - err3;
        c = self.splitter * adx;
        abig = c - adx;
        a1hi = c - abig;
        a1lo = adx - a1hi;
        _j = adx * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = adxtail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = adx * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        axcy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        axcy[5usize] = around + bround;
        axcy7 = _m + _k;
        bvirt = axcy7 - _m;
        avirt = axcy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        axcy[6usize] = around + bround;
        axcy[7usize] = axcy7;
        temp16len = fast_expansion_sum_zeroelim(
                        8i32,
                        bxcy.as_ptr(),
                        8i32,
                        cxby.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp32len = self.scale_expansion_zeroelim(
                        temp16len,
                        temp16.as_ptr(),
                        adz,
                        temp32.as_mut_ptr()
                    );
        temp32tlen = self.scale_expansion_zeroelim(
                         temp16len,
                         temp16.as_ptr(),
                         adztail,
                         temp32t.as_mut_ptr()
                     );
        alen = fast_expansion_sum_zeroelim(
                   temp32len,
                   temp32.as_ptr(),
                   temp32tlen,
                   temp32t.as_ptr(),
                   adet.as_mut_ptr()
               );
        temp16len = fast_expansion_sum_zeroelim(
                        8i32,
                        cxay.as_ptr(),
                        8i32,
                        axcy.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp32len = self.scale_expansion_zeroelim(
                        temp16len,
                        temp16.as_ptr(),
                        bdz,
                        temp32.as_mut_ptr()
                    );
        temp32tlen = self.scale_expansion_zeroelim(
                         temp16len,
                         temp16.as_ptr(),
                         bdztail,
                         temp32t.as_mut_ptr()
                     );
        blen = fast_expansion_sum_zeroelim(
                   temp32len,
                   temp32.as_ptr(),
                   temp32tlen,
                   temp32t.as_ptr(),
                   bdet.as_mut_ptr()
               );
        temp16len = fast_expansion_sum_zeroelim(
                        8i32,
                        axby.as_ptr(),
                        8i32,
                        bxay.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp32len = self.scale_expansion_zeroelim(
                        temp16len,
                        temp16.as_ptr(),
                        cdz,
                        temp32.as_mut_ptr()
                    );
        temp32tlen = self.scale_expansion_zeroelim(
                         temp16len,
                         temp16.as_ptr(),
                         cdztail,
                         temp32t.as_mut_ptr()
                     );
        clen = fast_expansion_sum_zeroelim(
                   temp32len,
                   temp32.as_ptr(),
                   temp32tlen,
                   temp32t.as_ptr(),
                   cdet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        deterlen = fast_expansion_sum_zeroelim(
                       ablen,
                       abdet.as_ptr(),
                       clen,
                       cdet.as_ptr(),
                       deter.as_mut_ptr()
                   );
        deter[(deterlen - 1i32) as (usize)]
    }

    
    pub unsafe fn orient3dadapt(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64,
        mut permanent : f64
    ) -> f64 {
        let mut adx : f64;
        let mut bdx : f64;
        let mut cdx : f64;
        let mut ady : f64;
        let mut bdy : f64;
        let mut cdy : f64;
        let mut adz : f64;
        let mut bdz : f64;
        let mut cdz : f64;
        let mut det : f64;
        let mut errbound : f64;
        let mut bdxcdy1 : f64;
        let mut cdxbdy1 : f64;
        let mut cdxady1 : f64;
        let mut adxcdy1 : f64;
        let mut adxbdy1 : f64;
        let mut bdxady1 : f64;
        let mut bdxcdy0 : f64;
        let mut cdxbdy0 : f64;
        let mut cdxady0 : f64;
        let mut adxcdy0 : f64;
        let mut adxbdy0 : f64;
        let mut bdxady0 : f64;
        let mut bc : [f64; 4];
        let mut ca : [f64; 4];
        let mut ab : [f64; 4];
        let mut bc3 : f64;
        let mut ca3 : f64;
        let mut ab3 : f64;
        let mut adet : [f64; 8];
        let mut bdet : [f64; 8];
        let mut cdet : [f64; 8];
        let mut alen : i32;
        let mut blen : i32;
        let mut clen : i32;
        let mut abdet : [f64; 16];
        let mut ablen : i32;
        let mut finnow : *mut f64;
        let mut finother : *mut f64;
        let mut finswap : *mut f64;
        let mut fin1 : [f64; 192];
        let mut fin2 : [f64; 192];
        let mut finlength : i32;
        let mut adxtail : f64;
        let mut bdxtail : f64;
        let mut cdxtail : f64;
        let mut adytail : f64;
        let mut bdytail : f64;
        let mut cdytail : f64;
        let mut adztail : f64;
        let mut bdztail : f64;
        let mut cdztail : f64;
        let mut at_blarge : f64;
        let mut at_clarge : f64;
        let mut bt_clarge : f64;
        let mut bt_alarge : f64;
        let mut ct_alarge : f64;
        let mut ct_blarge : f64;
        let mut at_b : [f64; 4];
        let mut at_c : [f64; 4];
        let mut bt_c : [f64; 4];
        let mut bt_a : [f64; 4];
        let mut ct_a : [f64; 4];
        let mut ct_b : [f64; 4];
        let mut at_blen : i32;
        let mut at_clen : i32;
        let mut bt_clen : i32;
        let mut bt_alen : i32;
        let mut ct_alen : i32;
        let mut ct_blen : i32;
        let mut bdxt_cdy1 : f64;
        let mut cdxt_bdy1 : f64;
        let mut cdxt_ady1 : f64;
        let mut adxt_cdy1 : f64;
        let mut adxt_bdy1 : f64;
        let mut bdxt_ady1 : f64;
        let mut bdxt_cdy0 : f64;
        let mut cdxt_bdy0 : f64;
        let mut cdxt_ady0 : f64;
        let mut adxt_cdy0 : f64;
        let mut adxt_bdy0 : f64;
        let mut bdxt_ady0 : f64;
        let mut bdyt_cdx1 : f64;
        let mut cdyt_bdx1 : f64;
        let mut cdyt_adx1 : f64;
        let mut adyt_cdx1 : f64;
        let mut adyt_bdx1 : f64;
        let mut bdyt_adx1 : f64;
        let mut bdyt_cdx0 : f64;
        let mut cdyt_bdx0 : f64;
        let mut cdyt_adx0 : f64;
        let mut adyt_cdx0 : f64;
        let mut adyt_bdx0 : f64;
        let mut bdyt_adx0 : f64;
        let mut bct : [f64; 8];
        let mut cat : [f64; 8];
        let mut abt : [f64; 8];
        let mut bctlen : i32;
        let mut catlen : i32;
        let mut abtlen : i32;
        let mut bdxt_cdyt1 : f64;
        let mut cdxt_bdyt1 : f64;
        let mut cdxt_adyt1 : f64;
        let mut adxt_cdyt1 : f64;
        let mut adxt_bdyt1 : f64;
        let mut bdxt_adyt1 : f64;
        let mut bdxt_cdyt0 : f64;
        let mut cdxt_bdyt0 : f64;
        let mut cdxt_adyt0 : f64;
        let mut adxt_cdyt0 : f64;
        let mut adxt_bdyt0 : f64;
        let mut bdxt_adyt0 : f64;
        let mut u : [f64; 4];
        let mut v : [f64; 12];
        let mut w : [f64; 16];
        let mut u3 : f64;
        let mut vlength : i32;
        let mut wlength : i32;
        let mut negate : f64;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _k : f64;
        let mut _0 : f64;

        adx = uninitialized();
        bdx = uninitialized();
        cdx = uninitialized();
        ady = uninitialized();
        bdy = uninitialized();
        cdy = uninitialized();
        adz = uninitialized();
        bdz = uninitialized();
        cdz = uninitialized();
        det = uninitialized();
        errbound = uninitialized();
        bdxcdy1 = uninitialized();
        cdxbdy1 = uninitialized();
        cdxady1 = uninitialized();
        adxcdy1 = uninitialized();
        adxbdy1 = uninitialized();
        bdxady1 = uninitialized();
        bdxcdy0 = uninitialized();
        cdxbdy0 = uninitialized();
        cdxady0 = uninitialized();
        adxcdy0 = uninitialized();
        adxbdy0 = uninitialized();
        bdxady0 = uninitialized();
        bc = uninitialized();
        ca = uninitialized();
        ab = uninitialized();
        bc3 = uninitialized();
        ca3 = uninitialized();
        ab3 = uninitialized();
        adet = uninitialized();
        bdet = uninitialized();
        cdet = uninitialized();
        alen = uninitialized();
        blen = uninitialized();
        clen = uninitialized();
        abdet = uninitialized();
        ablen = uninitialized();
        finnow = uninitialized();
        finother = uninitialized();
        finswap = uninitialized();
        fin1 = uninitialized();
        fin2 = uninitialized();
        finlength = uninitialized();
        adxtail = uninitialized();
        bdxtail = uninitialized();
        cdxtail = uninitialized();
        adytail = uninitialized();
        bdytail = uninitialized();
        cdytail = uninitialized();
        adztail = uninitialized();
        bdztail = uninitialized();
        cdztail = uninitialized();
        at_blarge = uninitialized();
        at_clarge = uninitialized();
        bt_clarge = uninitialized();
        bt_alarge = uninitialized();
        ct_alarge = uninitialized();
        ct_blarge = uninitialized();
        at_b = uninitialized();
        at_c = uninitialized();
        bt_c = uninitialized();
        bt_a = uninitialized();
        ct_a = uninitialized();
        ct_b = uninitialized();
        at_blen = uninitialized();
        at_clen = uninitialized();
        bt_clen = uninitialized();
        bt_alen = uninitialized();
        ct_alen = uninitialized();
        ct_blen = uninitialized();
        bdxt_cdy1 = uninitialized();
        cdxt_bdy1 = uninitialized();
        cdxt_ady1 = uninitialized();
        adxt_cdy1 = uninitialized();
        adxt_bdy1 = uninitialized();
        bdxt_ady1 = uninitialized();
        bdxt_cdy0 = uninitialized();
        cdxt_bdy0 = uninitialized();
        cdxt_ady0 = uninitialized();
        adxt_cdy0 = uninitialized();
        adxt_bdy0 = uninitialized();
        bdxt_ady0 = uninitialized();
        bdyt_cdx1 = uninitialized();
        cdyt_bdx1 = uninitialized();
        cdyt_adx1 = uninitialized();
        adyt_cdx1 = uninitialized();
        adyt_bdx1 = uninitialized();
        bdyt_adx1 = uninitialized();
        bdyt_cdx0 = uninitialized();
        cdyt_bdx0 = uninitialized();
        cdyt_adx0 = uninitialized();
        adyt_cdx0 = uninitialized();
        adyt_bdx0 = uninitialized();
        bdyt_adx0 = uninitialized();
        bct = uninitialized();
        cat = uninitialized();
        abt = uninitialized();
        bctlen = uninitialized();
        catlen = uninitialized();
        abtlen = uninitialized();
        bdxt_cdyt1 = uninitialized();
        cdxt_bdyt1 = uninitialized();
        cdxt_adyt1 = uninitialized();
        adxt_cdyt1 = uninitialized();
        adxt_bdyt1 = uninitialized();
        bdxt_adyt1 = uninitialized();
        bdxt_cdyt0 = uninitialized();
        cdxt_bdyt0 = uninitialized();
        cdxt_adyt0 = uninitialized();
        adxt_cdyt0 = uninitialized();
        adxt_bdyt0 = uninitialized();
        bdxt_adyt0 = uninitialized();
        u = uninitialized();
        v = uninitialized();
        w = uninitialized();
        u3 = uninitialized();
        vlength = uninitialized();
        wlength = uninitialized();
        negate = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        ahi = uninitialized();
        alo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _k = uninitialized();
        _0 = uninitialized();
        adx = *pa.offset(0isize) - *pd.offset(0isize);
        bdx = *pb.offset(0isize) - *pd.offset(0isize);
        cdx = *pc.offset(0isize) - *pd.offset(0isize);
        ady = *pa.offset(1isize) - *pd.offset(1isize);
        bdy = *pb.offset(1isize) - *pd.offset(1isize);
        cdy = *pc.offset(1isize) - *pd.offset(1isize);
        adz = *pa.offset(2isize) - *pd.offset(2isize);
        bdz = *pb.offset(2isize) - *pd.offset(2isize);
        cdz = *pc.offset(2isize) - *pd.offset(2isize);
        bdxcdy1 = bdx * cdy;
        c = self.splitter * bdx;
        abig = c - bdx;
        ahi = c - abig;
        alo = bdx - ahi;
        c = self.splitter * cdy;
        abig = c - cdy;
        bhi = c - abig;
        blo = cdy - bhi;
        err1 = bdxcdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bdxcdy0 = alo * blo - err3;
        cdxbdy1 = cdx * bdy;
        c = self.splitter * cdx;
        abig = c - cdx;
        ahi = c - abig;
        alo = cdx - ahi;
        c = self.splitter * bdy;
        abig = c - bdy;
        bhi = c - abig;
        blo = bdy - bhi;
        err1 = cdxbdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cdxbdy0 = alo * blo - err3;
        _i = bdxcdy0 - cdxbdy0;
        bvirt = bdxcdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cdxbdy0;
        around = bdxcdy0 - avirt;
        bc[0usize] = around + bround;
        _j = bdxcdy1 + _i;
        bvirt = _j - bdxcdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bdxcdy1 - avirt;
        _0 = around + bround;
        _i = _0 - cdxbdy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cdxbdy1;
        around = _0 - avirt;
        bc[1usize] = around + bround;
        bc3 = _j + _i;
        bvirt = bc3 - _j;
        avirt = bc3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bc[2usize] = around + bround;
        bc[3usize] = bc3;
        alen = self.scale_expansion_zeroelim(
                   4i32,
                   bc.as_ptr(),
                   adz,
                   adet.as_mut_ptr()
               );
        cdxady1 = cdx * ady;
        c = self.splitter * cdx;
        abig = c - cdx;
        ahi = c - abig;
        alo = cdx - ahi;
        c = self.splitter * ady;
        abig = c - ady;
        bhi = c - abig;
        blo = ady - bhi;
        err1 = cdxady1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cdxady0 = alo * blo - err3;
        adxcdy1 = adx * cdy;
        c = self.splitter * adx;
        abig = c - adx;
        ahi = c - abig;
        alo = adx - ahi;
        c = self.splitter * cdy;
        abig = c - cdy;
        bhi = c - abig;
        blo = cdy - bhi;
        err1 = adxcdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        adxcdy0 = alo * blo - err3;
        _i = cdxady0 - adxcdy0;
        bvirt = cdxady0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - adxcdy0;
        around = cdxady0 - avirt;
        ca[0usize] = around + bround;
        _j = cdxady1 + _i;
        bvirt = _j - cdxady1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = cdxady1 - avirt;
        _0 = around + bround;
        _i = _0 - adxcdy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - adxcdy1;
        around = _0 - avirt;
        ca[1usize] = around + bround;
        ca3 = _j + _i;
        bvirt = ca3 - _j;
        avirt = ca3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ca[2usize] = around + bround;
        ca[3usize] = ca3;
        blen = self.scale_expansion_zeroelim(
                   4i32,
                   ca.as_ptr(),
                   bdz,
                   bdet.as_mut_ptr()
               );
        adxbdy1 = adx * bdy;
        c = self.splitter * adx;
        abig = c - adx;
        ahi = c - abig;
        alo = adx - ahi;
        c = self.splitter * bdy;
        abig = c - bdy;
        bhi = c - abig;
        blo = bdy - bhi;
        err1 = adxbdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        adxbdy0 = alo * blo - err3;
        bdxady1 = bdx * ady;
        c = self.splitter * bdx;
        abig = c - bdx;
        ahi = c - abig;
        alo = bdx - ahi;
        c = self.splitter * ady;
        abig = c - ady;
        bhi = c - abig;
        blo = ady - bhi;
        err1 = bdxady1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bdxady0 = alo * blo - err3;
        _i = adxbdy0 - bdxady0;
        bvirt = adxbdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bdxady0;
        around = adxbdy0 - avirt;
        ab[0usize] = around + bround;
        _j = adxbdy1 + _i;
        bvirt = _j - adxbdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = adxbdy1 - avirt;
        _0 = around + bround;
        _i = _0 - bdxady1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bdxady1;
        around = _0 - avirt;
        ab[1usize] = around + bround;
        ab3 = _j + _i;
        bvirt = ab3 - _j;
        avirt = ab3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ab[2usize] = around + bround;
        ab[3usize] = ab3;
        clen = self.scale_expansion_zeroelim(
                   4i32,
                   ab.as_ptr(),
                   cdz,
                   cdet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        finlength = fast_expansion_sum_zeroelim(
                        ablen,
                        abdet.as_ptr(),
                        clen,
                        cdet.as_ptr(),
                        fin1.as_mut_ptr()
                    );
        det = estimate(finlength,fin1.as_ptr());
        errbound = self.o3derrboundB * permanent;
        if det >= errbound || -det >= errbound {
            det
        } else {
            bvirt = *pa.offset(0isize) - adx;
            avirt = adx + bvirt;
            bround = bvirt - *pd.offset(0isize);
            around = *pa.offset(0isize) - avirt;
            adxtail = around + bround;
            bvirt = *pb.offset(0isize) - bdx;
            avirt = bdx + bvirt;
            bround = bvirt - *pd.offset(0isize);
            around = *pb.offset(0isize) - avirt;
            bdxtail = around + bround;
            bvirt = *pc.offset(0isize) - cdx;
            avirt = cdx + bvirt;
            bround = bvirt - *pd.offset(0isize);
            around = *pc.offset(0isize) - avirt;
            cdxtail = around + bround;
            bvirt = *pa.offset(1isize) - ady;
            avirt = ady + bvirt;
            bround = bvirt - *pd.offset(1isize);
            around = *pa.offset(1isize) - avirt;
            adytail = around + bround;
            bvirt = *pb.offset(1isize) - bdy;
            avirt = bdy + bvirt;
            bround = bvirt - *pd.offset(1isize);
            around = *pb.offset(1isize) - avirt;
            bdytail = around + bround;
            bvirt = *pc.offset(1isize) - cdy;
            avirt = cdy + bvirt;
            bround = bvirt - *pd.offset(1isize);
            around = *pc.offset(1isize) - avirt;
            cdytail = around + bround;
            bvirt = *pa.offset(2isize) - adz;
            avirt = adz + bvirt;
            bround = bvirt - *pd.offset(2isize);
            around = *pa.offset(2isize) - avirt;
            adztail = around + bround;
            bvirt = *pb.offset(2isize) - bdz;
            avirt = bdz + bvirt;
            bround = bvirt - *pd.offset(2isize);
            around = *pb.offset(2isize) - avirt;
            bdztail = around + bround;
            bvirt = *pc.offset(2isize) - cdz;
            avirt = cdz + bvirt;
            bround = bvirt - *pd.offset(2isize);
            around = *pc.offset(2isize) - avirt;
            cdztail = around + bround;
            (if adxtail == 0.0f64 && (bdxtail == 0.0f64) && (cdxtail == 0.0f64) && (adytail == 0.0f64) && (bdytail == 0.0f64) && (cdytail == 0.0f64) && (adztail == 0.0f64) && (bdztail == 0.0f64) && (cdztail == 0.0f64) {
                 det
             } else {
                 errbound = self.o3derrboundC * permanent + self.resulterrbound * Absolute(
                                                                            det
                                                                        );
                 det = det + (adz * (bdx * cdytail + cdy * bdxtail - (bdy * cdxtail + cdx * bdytail)) + adztail * (bdx * cdy - bdy * cdx) + (bdz * (cdx * adytail + ady * cdxtail - (cdy * adxtail + adx * cdytail)) + bdztail * (cdx * ady - cdy * adx)) + (cdz * (adx * bdytail + bdy * adxtail - (ady * bdxtail + bdx * adytail)) + cdztail * (adx * bdy - ady * bdx)));
                 (if det >= errbound || -det >= errbound {
                      det
                  } else {
                      finnow = fin1.as_mut_ptr();
                      finother = fin2.as_mut_ptr();
                      if adxtail == 0.0f64 {
                          if adytail == 0.0f64 {
                              at_b[0usize] = 0.0f64;
                              at_blen = 1i32;
                              at_c[0usize] = 0.0f64;
                              at_clen = 1i32;
                          } else {
                              negate = -adytail;
                              at_blarge = negate * bdx;
                              c = self.splitter * negate;
                              abig = c - negate;
                              ahi = c - abig;
                              alo = negate - ahi;
                              c = self.splitter * bdx;
                              abig = c - bdx;
                              bhi = c - abig;
                              blo = bdx - bhi;
                              err1 = at_blarge - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              at_b[0usize] = alo * blo - err3;
                              at_b[1usize] = at_blarge;
                              at_blen = 2i32;
                              at_clarge = adytail * cdx;
                              c = self.splitter * adytail;
                              abig = c - adytail;
                              ahi = c - abig;
                              alo = adytail - ahi;
                              c = self.splitter * cdx;
                              abig = c - cdx;
                              bhi = c - abig;
                              blo = cdx - bhi;
                              err1 = at_clarge - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              at_c[0usize] = alo * blo - err3;
                              at_c[1usize] = at_clarge;
                              at_clen = 2i32;
                          }
                      } else if adytail == 0.0f64 {
                          at_blarge = adxtail * bdy;
                          c = self.splitter * adxtail;
                          abig = c - adxtail;
                          ahi = c - abig;
                          alo = adxtail - ahi;
                          c = self.splitter * bdy;
                          abig = c - bdy;
                          bhi = c - abig;
                          blo = bdy - bhi;
                          err1 = at_blarge - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          at_b[0usize] = alo * blo - err3;
                          at_b[1usize] = at_blarge;
                          at_blen = 2i32;
                          negate = -adxtail;
                          at_clarge = negate * cdy;
                          c = self.splitter * negate;
                          abig = c - negate;
                          ahi = c - abig;
                          alo = negate - ahi;
                          c = self.splitter * cdy;
                          abig = c - cdy;
                          bhi = c - abig;
                          blo = cdy - bhi;
                          err1 = at_clarge - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          at_c[0usize] = alo * blo - err3;
                          at_c[1usize] = at_clarge;
                          at_clen = 2i32;
                      } else {
                          adxt_bdy1 = adxtail * bdy;
                          c = self.splitter * adxtail;
                          abig = c - adxtail;
                          ahi = c - abig;
                          alo = adxtail - ahi;
                          c = self.splitter * bdy;
                          abig = c - bdy;
                          bhi = c - abig;
                          blo = bdy - bhi;
                          err1 = adxt_bdy1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          adxt_bdy0 = alo * blo - err3;
                          adyt_bdx1 = adytail * bdx;
                          c = self.splitter * adytail;
                          abig = c - adytail;
                          ahi = c - abig;
                          alo = adytail - ahi;
                          c = self.splitter * bdx;
                          abig = c - bdx;
                          bhi = c - abig;
                          blo = bdx - bhi;
                          err1 = adyt_bdx1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          adyt_bdx0 = alo * blo - err3;
                          _i = adxt_bdy0 - adyt_bdx0;
                          bvirt = adxt_bdy0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - adyt_bdx0;
                          around = adxt_bdy0 - avirt;
                          at_b[0usize] = around + bround;
                          _j = adxt_bdy1 + _i;
                          bvirt = _j - adxt_bdy1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = adxt_bdy1 - avirt;
                          _0 = around + bround;
                          _i = _0 - adyt_bdx1;
                          bvirt = _0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - adyt_bdx1;
                          around = _0 - avirt;
                          at_b[1usize] = around + bround;
                          at_blarge = _j + _i;
                          bvirt = at_blarge - _j;
                          avirt = at_blarge - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          at_b[2usize] = around + bround;
                          at_b[3usize] = at_blarge;
                          at_blen = 4i32;
                          adyt_cdx1 = adytail * cdx;
                          c = self.splitter * adytail;
                          abig = c - adytail;
                          ahi = c - abig;
                          alo = adytail - ahi;
                          c = self.splitter * cdx;
                          abig = c - cdx;
                          bhi = c - abig;
                          blo = cdx - bhi;
                          err1 = adyt_cdx1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          adyt_cdx0 = alo * blo - err3;
                          adxt_cdy1 = adxtail * cdy;
                          c = self.splitter * adxtail;
                          abig = c - adxtail;
                          ahi = c - abig;
                          alo = adxtail - ahi;
                          c = self.splitter * cdy;
                          abig = c - cdy;
                          bhi = c - abig;
                          blo = cdy - bhi;
                          err1 = adxt_cdy1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          adxt_cdy0 = alo * blo - err3;
                          _i = adyt_cdx0 - adxt_cdy0;
                          bvirt = adyt_cdx0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - adxt_cdy0;
                          around = adyt_cdx0 - avirt;
                          at_c[0usize] = around + bround;
                          _j = adyt_cdx1 + _i;
                          bvirt = _j - adyt_cdx1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = adyt_cdx1 - avirt;
                          _0 = around + bround;
                          _i = _0 - adxt_cdy1;
                          bvirt = _0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - adxt_cdy1;
                          around = _0 - avirt;
                          at_c[1usize] = around + bround;
                          at_clarge = _j + _i;
                          bvirt = at_clarge - _j;
                          avirt = at_clarge - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          at_c[2usize] = around + bround;
                          at_c[3usize] = at_clarge;
                          at_clen = 4i32;
                      }
                      if bdxtail == 0.0f64 {
                          if bdytail == 0.0f64 {
                              bt_c[0usize] = 0.0f64;
                              bt_clen = 1i32;
                              bt_a[0usize] = 0.0f64;
                              bt_alen = 1i32;
                          } else {
                              negate = -bdytail;
                              bt_clarge = negate * cdx;
                              c = self.splitter * negate;
                              abig = c - negate;
                              ahi = c - abig;
                              alo = negate - ahi;
                              c = self.splitter * cdx;
                              abig = c - cdx;
                              bhi = c - abig;
                              blo = cdx - bhi;
                              err1 = bt_clarge - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              bt_c[0usize] = alo * blo - err3;
                              bt_c[1usize] = bt_clarge;
                              bt_clen = 2i32;
                              bt_alarge = bdytail * adx;
                              c = self.splitter * bdytail;
                              abig = c - bdytail;
                              ahi = c - abig;
                              alo = bdytail - ahi;
                              c = self.splitter * adx;
                              abig = c - adx;
                              bhi = c - abig;
                              blo = adx - bhi;
                              err1 = bt_alarge - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              bt_a[0usize] = alo * blo - err3;
                              bt_a[1usize] = bt_alarge;
                              bt_alen = 2i32;
                          }
                      } else if bdytail == 0.0f64 {
                          bt_clarge = bdxtail * cdy;
                          c = self.splitter * bdxtail;
                          abig = c - bdxtail;
                          ahi = c - abig;
                          alo = bdxtail - ahi;
                          c = self.splitter * cdy;
                          abig = c - cdy;
                          bhi = c - abig;
                          blo = cdy - bhi;
                          err1 = bt_clarge - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          bt_c[0usize] = alo * blo - err3;
                          bt_c[1usize] = bt_clarge;
                          bt_clen = 2i32;
                          negate = -bdxtail;
                          bt_alarge = negate * ady;
                          c = self.splitter * negate;
                          abig = c - negate;
                          ahi = c - abig;
                          alo = negate - ahi;
                          c = self.splitter * ady;
                          abig = c - ady;
                          bhi = c - abig;
                          blo = ady - bhi;
                          err1 = bt_alarge - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          bt_a[0usize] = alo * blo - err3;
                          bt_a[1usize] = bt_alarge;
                          bt_alen = 2i32;
                      } else {
                          bdxt_cdy1 = bdxtail * cdy;
                          c = self.splitter * bdxtail;
                          abig = c - bdxtail;
                          ahi = c - abig;
                          alo = bdxtail - ahi;
                          c = self.splitter * cdy;
                          abig = c - cdy;
                          bhi = c - abig;
                          blo = cdy - bhi;
                          err1 = bdxt_cdy1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          bdxt_cdy0 = alo * blo - err3;
                          bdyt_cdx1 = bdytail * cdx;
                          c = self.splitter * bdytail;
                          abig = c - bdytail;
                          ahi = c - abig;
                          alo = bdytail - ahi;
                          c = self.splitter * cdx;
                          abig = c - cdx;
                          bhi = c - abig;
                          blo = cdx - bhi;
                          err1 = bdyt_cdx1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          bdyt_cdx0 = alo * blo - err3;
                          _i = bdxt_cdy0 - bdyt_cdx0;
                          bvirt = bdxt_cdy0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - bdyt_cdx0;
                          around = bdxt_cdy0 - avirt;
                          bt_c[0usize] = around + bround;
                          _j = bdxt_cdy1 + _i;
                          bvirt = _j - bdxt_cdy1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = bdxt_cdy1 - avirt;
                          _0 = around + bround;
                          _i = _0 - bdyt_cdx1;
                          bvirt = _0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - bdyt_cdx1;
                          around = _0 - avirt;
                          bt_c[1usize] = around + bround;
                          bt_clarge = _j + _i;
                          bvirt = bt_clarge - _j;
                          avirt = bt_clarge - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          bt_c[2usize] = around + bround;
                          bt_c[3usize] = bt_clarge;
                          bt_clen = 4i32;
                          bdyt_adx1 = bdytail * adx;
                          c = self.splitter * bdytail;
                          abig = c - bdytail;
                          ahi = c - abig;
                          alo = bdytail - ahi;
                          c = self.splitter * adx;
                          abig = c - adx;
                          bhi = c - abig;
                          blo = adx - bhi;
                          err1 = bdyt_adx1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          bdyt_adx0 = alo * blo - err3;
                          bdxt_ady1 = bdxtail * ady;
                          c = self.splitter * bdxtail;
                          abig = c - bdxtail;
                          ahi = c - abig;
                          alo = bdxtail - ahi;
                          c = self.splitter * ady;
                          abig = c - ady;
                          bhi = c - abig;
                          blo = ady - bhi;
                          err1 = bdxt_ady1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          bdxt_ady0 = alo * blo - err3;
                          _i = bdyt_adx0 - bdxt_ady0;
                          bvirt = bdyt_adx0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - bdxt_ady0;
                          around = bdyt_adx0 - avirt;
                          bt_a[0usize] = around + bround;
                          _j = bdyt_adx1 + _i;
                          bvirt = _j - bdyt_adx1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = bdyt_adx1 - avirt;
                          _0 = around + bround;
                          _i = _0 - bdxt_ady1;
                          bvirt = _0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - bdxt_ady1;
                          around = _0 - avirt;
                          bt_a[1usize] = around + bround;
                          bt_alarge = _j + _i;
                          bvirt = bt_alarge - _j;
                          avirt = bt_alarge - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          bt_a[2usize] = around + bround;
                          bt_a[3usize] = bt_alarge;
                          bt_alen = 4i32;
                      }
                      if cdxtail == 0.0f64 {
                          if cdytail == 0.0f64 {
                              ct_a[0usize] = 0.0f64;
                              ct_alen = 1i32;
                              ct_b[0usize] = 0.0f64;
                              ct_blen = 1i32;
                          } else {
                              negate = -cdytail;
                              ct_alarge = negate * adx;
                              c = self.splitter * negate;
                              abig = c - negate;
                              ahi = c - abig;
                              alo = negate - ahi;
                              c = self.splitter * adx;
                              abig = c - adx;
                              bhi = c - abig;
                              blo = adx - bhi;
                              err1 = ct_alarge - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ct_a[0usize] = alo * blo - err3;
                              ct_a[1usize] = ct_alarge;
                              ct_alen = 2i32;
                              ct_blarge = cdytail * bdx;
                              c = self.splitter * cdytail;
                              abig = c - cdytail;
                              ahi = c - abig;
                              alo = cdytail - ahi;
                              c = self.splitter * bdx;
                              abig = c - bdx;
                              bhi = c - abig;
                              blo = bdx - bhi;
                              err1 = ct_blarge - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ct_b[0usize] = alo * blo - err3;
                              ct_b[1usize] = ct_blarge;
                              ct_blen = 2i32;
                          }
                      } else if cdytail == 0.0f64 {
                          ct_alarge = cdxtail * ady;
                          c = self.splitter * cdxtail;
                          abig = c - cdxtail;
                          ahi = c - abig;
                          alo = cdxtail - ahi;
                          c = self.splitter * ady;
                          abig = c - ady;
                          bhi = c - abig;
                          blo = ady - bhi;
                          err1 = ct_alarge - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          ct_a[0usize] = alo * blo - err3;
                          ct_a[1usize] = ct_alarge;
                          ct_alen = 2i32;
                          negate = -cdxtail;
                          ct_blarge = negate * bdy;
                          c = self.splitter * negate;
                          abig = c - negate;
                          ahi = c - abig;
                          alo = negate - ahi;
                          c = self.splitter * bdy;
                          abig = c - bdy;
                          bhi = c - abig;
                          blo = bdy - bhi;
                          err1 = ct_blarge - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          ct_b[0usize] = alo * blo - err3;
                          ct_b[1usize] = ct_blarge;
                          ct_blen = 2i32;
                      } else {
                          cdxt_ady1 = cdxtail * ady;
                          c = self.splitter * cdxtail;
                          abig = c - cdxtail;
                          ahi = c - abig;
                          alo = cdxtail - ahi;
                          c = self.splitter * ady;
                          abig = c - ady;
                          bhi = c - abig;
                          blo = ady - bhi;
                          err1 = cdxt_ady1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          cdxt_ady0 = alo * blo - err3;
                          cdyt_adx1 = cdytail * adx;
                          c = self.splitter * cdytail;
                          abig = c - cdytail;
                          ahi = c - abig;
                          alo = cdytail - ahi;
                          c = self.splitter * adx;
                          abig = c - adx;
                          bhi = c - abig;
                          blo = adx - bhi;
                          err1 = cdyt_adx1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          cdyt_adx0 = alo * blo - err3;
                          _i = cdxt_ady0 - cdyt_adx0;
                          bvirt = cdxt_ady0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - cdyt_adx0;
                          around = cdxt_ady0 - avirt;
                          ct_a[0usize] = around + bround;
                          _j = cdxt_ady1 + _i;
                          bvirt = _j - cdxt_ady1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = cdxt_ady1 - avirt;
                          _0 = around + bround;
                          _i = _0 - cdyt_adx1;
                          bvirt = _0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - cdyt_adx1;
                          around = _0 - avirt;
                          ct_a[1usize] = around + bround;
                          ct_alarge = _j + _i;
                          bvirt = ct_alarge - _j;
                          avirt = ct_alarge - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          ct_a[2usize] = around + bround;
                          ct_a[3usize] = ct_alarge;
                          ct_alen = 4i32;
                          cdyt_bdx1 = cdytail * bdx;
                          c = self.splitter * cdytail;
                          abig = c - cdytail;
                          ahi = c - abig;
                          alo = cdytail - ahi;
                          c = self.splitter * bdx;
                          abig = c - bdx;
                          bhi = c - abig;
                          blo = bdx - bhi;
                          err1 = cdyt_bdx1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          cdyt_bdx0 = alo * blo - err3;
                          cdxt_bdy1 = cdxtail * bdy;
                          c = self.splitter * cdxtail;
                          abig = c - cdxtail;
                          ahi = c - abig;
                          alo = cdxtail - ahi;
                          c = self.splitter * bdy;
                          abig = c - bdy;
                          bhi = c - abig;
                          blo = bdy - bhi;
                          err1 = cdxt_bdy1 - ahi * bhi;
                          err2 = err1 - alo * bhi;
                          err3 = err2 - ahi * blo;
                          cdxt_bdy0 = alo * blo - err3;
                          _i = cdyt_bdx0 - cdxt_bdy0;
                          bvirt = cdyt_bdx0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - cdxt_bdy0;
                          around = cdyt_bdx0 - avirt;
                          ct_b[0usize] = around + bround;
                          _j = cdyt_bdx1 + _i;
                          bvirt = _j - cdyt_bdx1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = cdyt_bdx1 - avirt;
                          _0 = around + bround;
                          _i = _0 - cdxt_bdy1;
                          bvirt = _0 - _i;
                          avirt = _i + bvirt;
                          bround = bvirt - cdxt_bdy1;
                          around = _0 - avirt;
                          ct_b[1usize] = around + bround;
                          ct_blarge = _j + _i;
                          bvirt = ct_blarge - _j;
                          avirt = ct_blarge - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          ct_b[2usize] = around + bround;
                          ct_b[3usize] = ct_blarge;
                          ct_blen = 4i32;
                      }
                      bctlen = fast_expansion_sum_zeroelim(
                                   bt_clen,
                                   bt_c.as_ptr(),
                                   ct_blen,
                                   ct_b.as_ptr(),
                                   bct.as_mut_ptr()
                               );
                      wlength = self.scale_expansion_zeroelim(
                                    bctlen,
                                    bct.as_ptr(),
                                    adz,
                                    w.as_mut_ptr()
                                );
                      finlength = fast_expansion_sum_zeroelim(
                                      finlength,
                                      finnow as (*const f64),
                                      wlength,
                                      w.as_ptr(),
                                      finother
                                  );
                      finswap = finnow;
                      finnow = finother;
                      finother = finswap;
                      catlen = fast_expansion_sum_zeroelim(
                                   ct_alen,
                                   ct_a.as_ptr(),
                                   at_clen,
                                   at_c.as_ptr(),
                                   cat.as_mut_ptr()
                               );
                      wlength = self.scale_expansion_zeroelim(
                                    catlen,
                                    cat.as_ptr(),
                                    bdz,
                                    w.as_mut_ptr()
                                );
                      finlength = fast_expansion_sum_zeroelim(
                                      finlength,
                                      finnow as (*const f64),
                                      wlength,
                                      w.as_ptr(),
                                      finother
                                  );
                      finswap = finnow;
                      finnow = finother;
                      finother = finswap;
                      abtlen = fast_expansion_sum_zeroelim(
                                   at_blen,
                                   at_b.as_ptr(),
                                   bt_alen,
                                   bt_a.as_ptr(),
                                   abt.as_mut_ptr()
                               );
                      wlength = self.scale_expansion_zeroelim(
                                    abtlen,
                                    abt.as_ptr(),
                                    cdz,
                                    w.as_mut_ptr()
                                );
                      finlength = fast_expansion_sum_zeroelim(
                                      finlength,
                                      finnow as (*const f64),
                                      wlength,
                                      w.as_ptr(),
                                      finother
                                  );
                      finswap = finnow;
                      finnow = finother;
                      finother = finswap;
                      if adztail != 0.0f64 {
                          vlength = self.scale_expansion_zeroelim(
                                        4i32,
                                        bc.as_ptr(),
                                        adztail,
                                        v.as_mut_ptr()
                                    );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          vlength,
                                          v.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if bdztail != 0.0f64 {
                          vlength = self.scale_expansion_zeroelim(
                                        4i32,
                                        ca.as_ptr(),
                                        bdztail,
                                        v.as_mut_ptr()
                                    );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          vlength,
                                          v.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if cdztail != 0.0f64 {
                          vlength = self.scale_expansion_zeroelim(
                                        4i32,
                                        ab.as_ptr(),
                                        cdztail,
                                        v.as_mut_ptr()
                                    );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          vlength,
                                          v.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if adxtail != 0.0f64 {
                          if bdytail != 0.0f64 {
                              adxt_bdyt1 = adxtail * bdytail;
                              c = self.splitter * adxtail;
                              abig = c - adxtail;
                              ahi = c - abig;
                              alo = adxtail - ahi;
                              c = self.splitter * bdytail;
                              abig = c - bdytail;
                              bhi = c - abig;
                              blo = bdytail - bhi;
                              err1 = adxt_bdyt1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              adxt_bdyt0 = alo * blo - err3;
                              c = self.splitter * cdz;
                              abig = c - cdz;
                              bhi = c - abig;
                              blo = cdz - bhi;
                              _i = adxt_bdyt0 * cdz;
                              c = self.splitter * adxt_bdyt0;
                              abig = c - adxt_bdyt0;
                              ahi = c - abig;
                              alo = adxt_bdyt0 - ahi;
                              err1 = _i - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              u[0usize] = alo * blo - err3;
                              _j = adxt_bdyt1 * cdz;
                              c = self.splitter * adxt_bdyt1;
                              abig = c - adxt_bdyt1;
                              ahi = c - abig;
                              alo = adxt_bdyt1 - ahi;
                              err1 = _j - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              _0 = alo * blo - err3;
                              _k = _i + _0;
                              bvirt = _k - _i;
                              avirt = _k - bvirt;
                              bround = _0 - bvirt;
                              around = _i - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _k;
                              bvirt = u3 - _j;
                              u[2usize] = _k - bvirt;
                              u[3usize] = u3;
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              4i32,
                                              u.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if cdztail != 0.0f64 {
                                  c = self.splitter * cdztail;
                                  abig = c - cdztail;
                                  bhi = c - abig;
                                  blo = cdztail - bhi;
                                  _i = adxt_bdyt0 * cdztail;
                                  c = self.splitter * adxt_bdyt0;
                                  abig = c - adxt_bdyt0;
                                  ahi = c - abig;
                                  alo = adxt_bdyt0 - ahi;
                                  err1 = _i - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  u[0usize] = alo * blo - err3;
                                  _j = adxt_bdyt1 * cdztail;
                                  c = self.splitter * adxt_bdyt1;
                                  abig = c - adxt_bdyt1;
                                  ahi = c - abig;
                                  alo = adxt_bdyt1 - ahi;
                                  err1 = _j - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  _0 = alo * blo - err3;
                                  _k = _i + _0;
                                  bvirt = _k - _i;
                                  avirt = _k - bvirt;
                                  bround = _0 - bvirt;
                                  around = _i - avirt;
                                  u[1usize] = around + bround;
                                  u3 = _j + _k;
                                  bvirt = u3 - _j;
                                  u[2usize] = _k - bvirt;
                                  u[3usize] = u3;
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  4i32,
                                                  u.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                          }
                          if cdytail != 0.0f64 {
                              negate = -adxtail;
                              adxt_cdyt1 = negate * cdytail;
                              c = self.splitter * negate;
                              abig = c - negate;
                              ahi = c - abig;
                              alo = negate - ahi;
                              c = self.splitter * cdytail;
                              abig = c - cdytail;
                              bhi = c - abig;
                              blo = cdytail - bhi;
                              err1 = adxt_cdyt1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              adxt_cdyt0 = alo * blo - err3;
                              c = self.splitter * bdz;
                              abig = c - bdz;
                              bhi = c - abig;
                              blo = bdz - bhi;
                              _i = adxt_cdyt0 * bdz;
                              c = self.splitter * adxt_cdyt0;
                              abig = c - adxt_cdyt0;
                              ahi = c - abig;
                              alo = adxt_cdyt0 - ahi;
                              err1 = _i - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              u[0usize] = alo * blo - err3;
                              _j = adxt_cdyt1 * bdz;
                              c = self.splitter * adxt_cdyt1;
                              abig = c - adxt_cdyt1;
                              ahi = c - abig;
                              alo = adxt_cdyt1 - ahi;
                              err1 = _j - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              _0 = alo * blo - err3;
                              _k = _i + _0;
                              bvirt = _k - _i;
                              avirt = _k - bvirt;
                              bround = _0 - bvirt;
                              around = _i - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _k;
                              bvirt = u3 - _j;
                              u[2usize] = _k - bvirt;
                              u[3usize] = u3;
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              4i32,
                                              u.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if bdztail != 0.0f64 {
                                  c = self.splitter * bdztail;
                                  abig = c - bdztail;
                                  bhi = c - abig;
                                  blo = bdztail - bhi;
                                  _i = adxt_cdyt0 * bdztail;
                                  c = self.splitter * adxt_cdyt0;
                                  abig = c - adxt_cdyt0;
                                  ahi = c - abig;
                                  alo = adxt_cdyt0 - ahi;
                                  err1 = _i - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  u[0usize] = alo * blo - err3;
                                  _j = adxt_cdyt1 * bdztail;
                                  c = self.splitter * adxt_cdyt1;
                                  abig = c - adxt_cdyt1;
                                  ahi = c - abig;
                                  alo = adxt_cdyt1 - ahi;
                                  err1 = _j - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  _0 = alo * blo - err3;
                                  _k = _i + _0;
                                  bvirt = _k - _i;
                                  avirt = _k - bvirt;
                                  bround = _0 - bvirt;
                                  around = _i - avirt;
                                  u[1usize] = around + bround;
                                  u3 = _j + _k;
                                  bvirt = u3 - _j;
                                  u[2usize] = _k - bvirt;
                                  u[3usize] = u3;
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  4i32,
                                                  u.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                          }
                      }
                      if bdxtail != 0.0f64 {
                          if cdytail != 0.0f64 {
                              bdxt_cdyt1 = bdxtail * cdytail;
                              c = self.splitter * bdxtail;
                              abig = c - bdxtail;
                              ahi = c - abig;
                              alo = bdxtail - ahi;
                              c = self.splitter * cdytail;
                              abig = c - cdytail;
                              bhi = c - abig;
                              blo = cdytail - bhi;
                              err1 = bdxt_cdyt1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              bdxt_cdyt0 = alo * blo - err3;
                              c = self.splitter * adz;
                              abig = c - adz;
                              bhi = c - abig;
                              blo = adz - bhi;
                              _i = bdxt_cdyt0 * adz;
                              c = self.splitter * bdxt_cdyt0;
                              abig = c - bdxt_cdyt0;
                              ahi = c - abig;
                              alo = bdxt_cdyt0 - ahi;
                              err1 = _i - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              u[0usize] = alo * blo - err3;
                              _j = bdxt_cdyt1 * adz;
                              c = self.splitter * bdxt_cdyt1;
                              abig = c - bdxt_cdyt1;
                              ahi = c - abig;
                              alo = bdxt_cdyt1 - ahi;
                              err1 = _j - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              _0 = alo * blo - err3;
                              _k = _i + _0;
                              bvirt = _k - _i;
                              avirt = _k - bvirt;
                              bround = _0 - bvirt;
                              around = _i - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _k;
                              bvirt = u3 - _j;
                              u[2usize] = _k - bvirt;
                              u[3usize] = u3;
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              4i32,
                                              u.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if adztail != 0.0f64 {
                                  c = self.splitter * adztail;
                                  abig = c - adztail;
                                  bhi = c - abig;
                                  blo = adztail - bhi;
                                  _i = bdxt_cdyt0 * adztail;
                                  c = self.splitter * bdxt_cdyt0;
                                  abig = c - bdxt_cdyt0;
                                  ahi = c - abig;
                                  alo = bdxt_cdyt0 - ahi;
                                  err1 = _i - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  u[0usize] = alo * blo - err3;
                                  _j = bdxt_cdyt1 * adztail;
                                  c = self.splitter * bdxt_cdyt1;
                                  abig = c - bdxt_cdyt1;
                                  ahi = c - abig;
                                  alo = bdxt_cdyt1 - ahi;
                                  err1 = _j - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  _0 = alo * blo - err3;
                                  _k = _i + _0;
                                  bvirt = _k - _i;
                                  avirt = _k - bvirt;
                                  bround = _0 - bvirt;
                                  around = _i - avirt;
                                  u[1usize] = around + bround;
                                  u3 = _j + _k;
                                  bvirt = u3 - _j;
                                  u[2usize] = _k - bvirt;
                                  u[3usize] = u3;
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  4i32,
                                                  u.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                          }
                          if adytail != 0.0f64 {
                              negate = -bdxtail;
                              bdxt_adyt1 = negate * adytail;
                              c = self.splitter * negate;
                              abig = c - negate;
                              ahi = c - abig;
                              alo = negate - ahi;
                              c = self.splitter * adytail;
                              abig = c - adytail;
                              bhi = c - abig;
                              blo = adytail - bhi;
                              err1 = bdxt_adyt1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              bdxt_adyt0 = alo * blo - err3;
                              c = self.splitter * cdz;
                              abig = c - cdz;
                              bhi = c - abig;
                              blo = cdz - bhi;
                              _i = bdxt_adyt0 * cdz;
                              c = self.splitter * bdxt_adyt0;
                              abig = c - bdxt_adyt0;
                              ahi = c - abig;
                              alo = bdxt_adyt0 - ahi;
                              err1 = _i - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              u[0usize] = alo * blo - err3;
                              _j = bdxt_adyt1 * cdz;
                              c = self.splitter * bdxt_adyt1;
                              abig = c - bdxt_adyt1;
                              ahi = c - abig;
                              alo = bdxt_adyt1 - ahi;
                              err1 = _j - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              _0 = alo * blo - err3;
                              _k = _i + _0;
                              bvirt = _k - _i;
                              avirt = _k - bvirt;
                              bround = _0 - bvirt;
                              around = _i - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _k;
                              bvirt = u3 - _j;
                              u[2usize] = _k - bvirt;
                              u[3usize] = u3;
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              4i32,
                                              u.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if cdztail != 0.0f64 {
                                  c = self.splitter * cdztail;
                                  abig = c - cdztail;
                                  bhi = c - abig;
                                  blo = cdztail - bhi;
                                  _i = bdxt_adyt0 * cdztail;
                                  c = self.splitter * bdxt_adyt0;
                                  abig = c - bdxt_adyt0;
                                  ahi = c - abig;
                                  alo = bdxt_adyt0 - ahi;
                                  err1 = _i - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  u[0usize] = alo * blo - err3;
                                  _j = bdxt_adyt1 * cdztail;
                                  c = self.splitter * bdxt_adyt1;
                                  abig = c - bdxt_adyt1;
                                  ahi = c - abig;
                                  alo = bdxt_adyt1 - ahi;
                                  err1 = _j - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  _0 = alo * blo - err3;
                                  _k = _i + _0;
                                  bvirt = _k - _i;
                                  avirt = _k - bvirt;
                                  bround = _0 - bvirt;
                                  around = _i - avirt;
                                  u[1usize] = around + bround;
                                  u3 = _j + _k;
                                  bvirt = u3 - _j;
                                  u[2usize] = _k - bvirt;
                                  u[3usize] = u3;
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  4i32,
                                                  u.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                          }
                      }
                      if cdxtail != 0.0f64 {
                          if adytail != 0.0f64 {
                              cdxt_adyt1 = cdxtail * adytail;
                              c = self.splitter * cdxtail;
                              abig = c - cdxtail;
                              ahi = c - abig;
                              alo = cdxtail - ahi;
                              c = self.splitter * adytail;
                              abig = c - adytail;
                              bhi = c - abig;
                              blo = adytail - bhi;
                              err1 = cdxt_adyt1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              cdxt_adyt0 = alo * blo - err3;
                              c = self.splitter * bdz;
                              abig = c - bdz;
                              bhi = c - abig;
                              blo = bdz - bhi;
                              _i = cdxt_adyt0 * bdz;
                              c = self.splitter * cdxt_adyt0;
                              abig = c - cdxt_adyt0;
                              ahi = c - abig;
                              alo = cdxt_adyt0 - ahi;
                              err1 = _i - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              u[0usize] = alo * blo - err3;
                              _j = cdxt_adyt1 * bdz;
                              c = self.splitter * cdxt_adyt1;
                              abig = c - cdxt_adyt1;
                              ahi = c - abig;
                              alo = cdxt_adyt1 - ahi;
                              err1 = _j - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              _0 = alo * blo - err3;
                              _k = _i + _0;
                              bvirt = _k - _i;
                              avirt = _k - bvirt;
                              bround = _0 - bvirt;
                              around = _i - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _k;
                              bvirt = u3 - _j;
                              u[2usize] = _k - bvirt;
                              u[3usize] = u3;
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              4i32,
                                              u.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if bdztail != 0.0f64 {
                                  c = self.splitter * bdztail;
                                  abig = c - bdztail;
                                  bhi = c - abig;
                                  blo = bdztail - bhi;
                                  _i = cdxt_adyt0 * bdztail;
                                  c = self.splitter * cdxt_adyt0;
                                  abig = c - cdxt_adyt0;
                                  ahi = c - abig;
                                  alo = cdxt_adyt0 - ahi;
                                  err1 = _i - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  u[0usize] = alo * blo - err3;
                                  _j = cdxt_adyt1 * bdztail;
                                  c = self.splitter * cdxt_adyt1;
                                  abig = c - cdxt_adyt1;
                                  ahi = c - abig;
                                  alo = cdxt_adyt1 - ahi;
                                  err1 = _j - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  _0 = alo * blo - err3;
                                  _k = _i + _0;
                                  bvirt = _k - _i;
                                  avirt = _k - bvirt;
                                  bround = _0 - bvirt;
                                  around = _i - avirt;
                                  u[1usize] = around + bround;
                                  u3 = _j + _k;
                                  bvirt = u3 - _j;
                                  u[2usize] = _k - bvirt;
                                  u[3usize] = u3;
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  4i32,
                                                  u.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                          }
                          if bdytail != 0.0f64 {
                              negate = -cdxtail;
                              cdxt_bdyt1 = negate * bdytail;
                              c = self.splitter * negate;
                              abig = c - negate;
                              ahi = c - abig;
                              alo = negate - ahi;
                              c = self.splitter * bdytail;
                              abig = c - bdytail;
                              bhi = c - abig;
                              blo = bdytail - bhi;
                              err1 = cdxt_bdyt1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              cdxt_bdyt0 = alo * blo - err3;
                              c = self.splitter * adz;
                              abig = c - adz;
                              bhi = c - abig;
                              blo = adz - bhi;
                              _i = cdxt_bdyt0 * adz;
                              c = self.splitter * cdxt_bdyt0;
                              abig = c - cdxt_bdyt0;
                              ahi = c - abig;
                              alo = cdxt_bdyt0 - ahi;
                              err1 = _i - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              u[0usize] = alo * blo - err3;
                              _j = cdxt_bdyt1 * adz;
                              c = self.splitter * cdxt_bdyt1;
                              abig = c - cdxt_bdyt1;
                              ahi = c - abig;
                              alo = cdxt_bdyt1 - ahi;
                              err1 = _j - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              _0 = alo * blo - err3;
                              _k = _i + _0;
                              bvirt = _k - _i;
                              avirt = _k - bvirt;
                              bround = _0 - bvirt;
                              around = _i - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _k;
                              bvirt = u3 - _j;
                              u[2usize] = _k - bvirt;
                              u[3usize] = u3;
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              4i32,
                                              u.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if adztail != 0.0f64 {
                                  c = self.splitter * adztail;
                                  abig = c - adztail;
                                  bhi = c - abig;
                                  blo = adztail - bhi;
                                  _i = cdxt_bdyt0 * adztail;
                                  c = self.splitter * cdxt_bdyt0;
                                  abig = c - cdxt_bdyt0;
                                  ahi = c - abig;
                                  alo = cdxt_bdyt0 - ahi;
                                  err1 = _i - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  u[0usize] = alo * blo - err3;
                                  _j = cdxt_bdyt1 * adztail;
                                  c = self.splitter * cdxt_bdyt1;
                                  abig = c - cdxt_bdyt1;
                                  ahi = c - abig;
                                  alo = cdxt_bdyt1 - ahi;
                                  err1 = _j - ahi * bhi;
                                  err2 = err1 - alo * bhi;
                                  err3 = err2 - ahi * blo;
                                  _0 = alo * blo - err3;
                                  _k = _i + _0;
                                  bvirt = _k - _i;
                                  avirt = _k - bvirt;
                                  bround = _0 - bvirt;
                                  around = _i - avirt;
                                  u[1usize] = around + bround;
                                  u3 = _j + _k;
                                  bvirt = u3 - _j;
                                  u[2usize] = _k - bvirt;
                                  u[3usize] = u3;
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  4i32,
                                                  u.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                          }
                      }
                      if adztail != 0.0f64 {
                          wlength = self.scale_expansion_zeroelim(
                                        bctlen,
                                        bct.as_ptr(),
                                        adztail,
                                        w.as_mut_ptr()
                                    );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          wlength,
                                          w.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if bdztail != 0.0f64 {
                          wlength = self.scale_expansion_zeroelim(
                                        catlen,
                                        cat.as_ptr(),
                                        bdztail,
                                        w.as_mut_ptr()
                                    );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          wlength,
                                          w.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if cdztail != 0.0f64 {
                          wlength = self.scale_expansion_zeroelim(
                                        abtlen,
                                        abt.as_ptr(),
                                        cdztail,
                                        w.as_mut_ptr()
                                    );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          wlength,
                                          w.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      *finnow.offset((finlength - 1i32) as (isize))
                  })
             })
        }
    }

    
    pub unsafe fn orient3d(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64
    ) -> f64 {
        let mut adx : f64;
        let mut bdx : f64;
        let mut cdx : f64;
        let mut ady : f64;
        let mut bdy : f64;
        let mut cdy : f64;
        let mut adz : f64;
        let mut bdz : f64;
        let mut cdz : f64;
        let mut bdxcdy : f64;
        let mut cdxbdy : f64;
        let mut cdxady : f64;
        let mut adxcdy : f64;
        let mut adxbdy : f64;
        let mut bdxady : f64;
        let mut det : f64;
        let mut permanent : f64;
        let mut errbound : f64;
        adx = *pa.offset(0isize) - *pd.offset(0isize);
        bdx = *pb.offset(0isize) - *pd.offset(0isize);
        cdx = *pc.offset(0isize) - *pd.offset(0isize);
        ady = *pa.offset(1isize) - *pd.offset(1isize);
        bdy = *pb.offset(1isize) - *pd.offset(1isize);
        cdy = *pc.offset(1isize) - *pd.offset(1isize);
        adz = *pa.offset(2isize) - *pd.offset(2isize);
        bdz = *pb.offset(2isize) - *pd.offset(2isize);
        cdz = *pc.offset(2isize) - *pd.offset(2isize);
        bdxcdy = bdx * cdy;
        cdxbdy = cdx * bdy;
        cdxady = cdx * ady;
        adxcdy = adx * cdy;
        adxbdy = adx * bdy;
        bdxady = bdx * ady;
        det = adz * (bdxcdy - cdxbdy) + bdz * (cdxady - adxcdy) + cdz * (adxbdy - bdxady);
        permanent = (Absolute(bdxcdy) + Absolute(cdxbdy)) * Absolute(
                                                                adz
                                                            ) + (Absolute(cdxady) + Absolute(
                                                                                        adxcdy
                                                                                    )) * Absolute(
                                                                                             bdz
                                                                                         ) + (Absolute(
                                                                                                  adxbdy
                                                                                              ) + Absolute(
                                                                                                      bdxady
                                                                                                  )) * Absolute(
                                                                                                           cdz
                                                                                                       );
        errbound = self.o3derrboundA * permanent;
        if det > errbound || -det > errbound {
            det
        } else {
            self.orient3dadapt(pa,pb,pc,pd,permanent)
        }
    }

    
    pub unsafe fn incircleexact(&self,
        mut pa : *mut f64,
        mut pb : *mut f64,
        mut pc : *mut f64,
        mut pd : *mut f64
    ) -> f64 {
        let mut axby1 : f64;
        let mut bxcy1 : f64;
        let mut cxdy1 : f64;
        let mut dxay1 : f64;
        let mut axcy1 : f64;
        let mut bxdy1 : f64;
        let mut bxay1 : f64;
        let mut cxby1 : f64;
        let mut dxcy1 : f64;
        let mut axdy1 : f64;
        let mut cxay1 : f64;
        let mut dxby1 : f64;
        let mut axby0 : f64;
        let mut bxcy0 : f64;
        let mut cxdy0 : f64;
        let mut dxay0 : f64;
        let mut axcy0 : f64;
        let mut bxdy0 : f64;
        let mut bxay0 : f64;
        let mut cxby0 : f64;
        let mut dxcy0 : f64;
        let mut axdy0 : f64;
        let mut cxay0 : f64;
        let mut dxby0 : f64;
        let mut ab : [f64; 4];
        let mut bc : [f64; 4];
        let mut cd : [f64; 4];
        let mut da : [f64; 4];
        let mut ac : [f64; 4];
        let mut bd : [f64; 4];
        let mut temp8 : [f64; 8];
        let mut templen : i32;
        let mut abc : [f64; 12];
        let mut bcd : [f64; 12];
        let mut cda : [f64; 12];
        let mut dab : [f64; 12];
        let mut abclen : i32;
        let mut bcdlen : i32;
        let mut cdalen : i32;
        let mut dablen : i32;
        let mut det24x : [f64; 24];
        let mut det24y : [f64; 24];
        let mut det48x : [f64; 48];
        let mut det48y : [f64; 48];
        let mut xlen : i32;
        let mut ylen : i32;
        let mut adet : [f64; 96];
        let mut bdet : [f64; 96];
        let mut cdet : [f64; 96];
        let mut ddet : [f64; 96];
        let mut alen : i32;
        let mut blen : i32;
        let mut clen : i32;
        let mut dlen : i32;
        let mut abdet : [f64; 192];
        let mut cddet : [f64; 192];
        let mut ablen : i32;
        let mut cdlen : i32;
        let mut deter : [f64; 384];
        let mut deterlen : i32;
        let mut i : i32;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _0 : f64;

        axby1 = uninitialized();
        bxcy1 = uninitialized();
        cxdy1 = uninitialized();
        dxay1 = uninitialized();
        axcy1 = uninitialized();
        bxdy1 = uninitialized();
        bxay1 = uninitialized();
        cxby1 = uninitialized();
        dxcy1 = uninitialized();
        axdy1 = uninitialized();
        cxay1 = uninitialized();
        dxby1 = uninitialized();
        axby0 = uninitialized();
        bxcy0 = uninitialized();
        cxdy0 = uninitialized();
        dxay0 = uninitialized();
        axcy0 = uninitialized();
        bxdy0 = uninitialized();
        bxay0 = uninitialized();
        cxby0 = uninitialized();
        dxcy0 = uninitialized();
        axdy0 = uninitialized();
        cxay0 = uninitialized();
        dxby0 = uninitialized();
        ab = uninitialized();
        bc = uninitialized();
        cd = uninitialized();
        da = uninitialized();
        ac = uninitialized();
        bd = uninitialized();
        temp8 = uninitialized();
        templen = uninitialized();
        abc = uninitialized();
        bcd = uninitialized();
        cda = uninitialized();
        dab = uninitialized();
        abclen = uninitialized();
        bcdlen = uninitialized();
        cdalen = uninitialized();
        dablen = uninitialized();
        det24x = uninitialized();
        det24y = uninitialized();
        det48x = uninitialized();
        det48y = uninitialized();
        xlen = uninitialized();
        ylen = uninitialized();
        adet = uninitialized();
        bdet = uninitialized();
        cdet = uninitialized();
        ddet = uninitialized();
        alen = uninitialized();
        blen = uninitialized();
        clen = uninitialized();
        dlen = uninitialized();
        abdet = uninitialized();
        cddet = uninitialized();
        ablen = uninitialized();
        cdlen = uninitialized();
        deter = uninitialized();
        deterlen = uninitialized();
        i = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        ahi = uninitialized();
        alo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _0 = uninitialized();
        axby1 = *pa.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = axby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axby0 = alo * blo - err3;
        bxay1 = *pb.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = bxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxay0 = alo * blo - err3;
        _i = axby0 - bxay0;
        bvirt = axby0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxay0;
        around = axby0 - avirt;
        ab[0usize] = around + bround;
        _j = axby1 + _i;
        bvirt = _j - axby1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = axby1 - avirt;
        _0 = around + bround;
        _i = _0 - bxay1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxay1;
        around = _0 - avirt;
        ab[1usize] = around + bround;
        ab[3usize] = _j + _i;
        bvirt = ab[3usize] - _j;
        avirt = ab[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ab[2usize] = around + bround;
        bxcy1 = *pb.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = bxcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxcy0 = alo * blo - err3;
        cxby1 = *pc.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = cxby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxby0 = alo * blo - err3;
        _i = bxcy0 - cxby0;
        bvirt = bxcy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxby0;
        around = bxcy0 - avirt;
        bc[0usize] = around + bround;
        _j = bxcy1 + _i;
        bvirt = _j - bxcy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bxcy1 - avirt;
        _0 = around + bround;
        _i = _0 - cxby1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxby1;
        around = _0 - avirt;
        bc[1usize] = around + bround;
        bc[3usize] = _j + _i;
        bvirt = bc[3usize] - _j;
        avirt = bc[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bc[2usize] = around + bround;
        cxdy1 = *pc.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = cxdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxdy0 = alo * blo - err3;
        dxcy1 = *pd.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = dxcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxcy0 = alo * blo - err3;
        _i = cxdy0 - dxcy0;
        bvirt = cxdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxcy0;
        around = cxdy0 - avirt;
        cd[0usize] = around + bround;
        _j = cxdy1 + _i;
        bvirt = _j - cxdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = cxdy1 - avirt;
        _0 = around + bround;
        _i = _0 - dxcy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxcy1;
        around = _0 - avirt;
        cd[1usize] = around + bround;
        cd[3usize] = _j + _i;
        bvirt = cd[3usize] - _j;
        avirt = cd[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        cd[2usize] = around + bround;
        dxay1 = *pd.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = dxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxay0 = alo * blo - err3;
        axdy1 = *pa.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = axdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axdy0 = alo * blo - err3;
        _i = dxay0 - axdy0;
        bvirt = dxay0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axdy0;
        around = dxay0 - avirt;
        da[0usize] = around + bround;
        _j = dxay1 + _i;
        bvirt = _j - dxay1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = dxay1 - avirt;
        _0 = around + bround;
        _i = _0 - axdy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axdy1;
        around = _0 - avirt;
        da[1usize] = around + bround;
        da[3usize] = _j + _i;
        bvirt = da[3usize] - _j;
        avirt = da[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        da[2usize] = around + bround;
        axcy1 = *pa.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = axcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axcy0 = alo * blo - err3;
        cxay1 = *pc.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = cxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxay0 = alo * blo - err3;
        _i = axcy0 - cxay0;
        bvirt = axcy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxay0;
        around = axcy0 - avirt;
        ac[0usize] = around + bround;
        _j = axcy1 + _i;
        bvirt = _j - axcy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = axcy1 - avirt;
        _0 = around + bround;
        _i = _0 - cxay1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxay1;
        around = _0 - avirt;
        ac[1usize] = around + bround;
        ac[3usize] = _j + _i;
        bvirt = ac[3usize] - _j;
        avirt = ac[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ac[2usize] = around + bround;
        bxdy1 = *pb.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = bxdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxdy0 = alo * blo - err3;
        dxby1 = *pd.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = dxby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxby0 = alo * blo - err3;
        _i = bxdy0 - dxby0;
        bvirt = bxdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxby0;
        around = bxdy0 - avirt;
        bd[0usize] = around + bround;
        _j = bxdy1 + _i;
        bvirt = _j - bxdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bxdy1 - avirt;
        _0 = around + bround;
        _i = _0 - dxby1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxby1;
        around = _0 - avirt;
        bd[1usize] = around + bround;
        bd[3usize] = _j + _i;
        bvirt = bd[3usize] - _j;
        avirt = bd[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bd[2usize] = around + bround;
        templen = fast_expansion_sum_zeroelim(
                      4i32,
                      cd.as_ptr(),
                      4i32,
                      da.as_ptr(),
                      temp8.as_mut_ptr()
                  );
        cdalen = fast_expansion_sum_zeroelim(
                     templen,
                     temp8.as_ptr(),
                     4i32,
                     ac.as_ptr(),
                     cda.as_mut_ptr()
                 );
        templen = fast_expansion_sum_zeroelim(
                      4i32,
                      da.as_ptr(),
                      4i32,
                      ab.as_ptr(),
                      temp8.as_mut_ptr()
                  );
        dablen = fast_expansion_sum_zeroelim(
                     templen,
                     temp8.as_ptr(),
                     4i32,
                     bd.as_ptr(),
                     dab.as_mut_ptr()
                 );
        i = 0i32;
        'loop1: loop {
            if !(i < 4i32) {
                break;
            }
            bd[i as (usize)] = -bd[i as (usize)];
            ac[i as (usize)] = -ac[i as (usize)];
            i = i + 1;
        }
        templen = fast_expansion_sum_zeroelim(
                      4i32,
                      ab.as_ptr(),
                      4i32,
                      bc.as_ptr(),
                      temp8.as_mut_ptr()
                  );
        abclen = fast_expansion_sum_zeroelim(
                     templen,
                     temp8.as_ptr(),
                     4i32,
                     ac.as_ptr(),
                     abc.as_mut_ptr()
                 );
        templen = fast_expansion_sum_zeroelim(
                      4i32,
                      bc.as_ptr(),
                      4i32,
                      cd.as_ptr(),
                      temp8.as_mut_ptr()
                  );
        bcdlen = fast_expansion_sum_zeroelim(
                     templen,
                     temp8.as_ptr(),
                     4i32,
                     bd.as_ptr(),
                     bcd.as_mut_ptr()
                 );
        xlen = self.scale_expansion_zeroelim(
                   bcdlen,
                   bcd.as_ptr(),
                   *pa.offset(0isize),
                   det24x.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   det24x.as_ptr(),
                   *pa.offset(0isize),
                   det48x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   bcdlen,
                   bcd.as_ptr(),
                   *pa.offset(1isize),
                   det24y.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   det24y.as_ptr(),
                   *pa.offset(1isize),
                   det48y.as_mut_ptr()
               );
        alen = fast_expansion_sum_zeroelim(
                   xlen,
                   det48x.as_ptr(),
                   ylen,
                   det48y.as_ptr(),
                   adet.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   cdalen,
                   cda.as_ptr(),
                   *pb.offset(0isize),
                   det24x.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   det24x.as_ptr(),
                   -*pb.offset(0isize),
                   det48x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   cdalen,
                   cda.as_ptr(),
                   *pb.offset(1isize),
                   det24y.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   det24y.as_ptr(),
                   -*pb.offset(1isize),
                   det48y.as_mut_ptr()
               );
        blen = fast_expansion_sum_zeroelim(
                   xlen,
                   det48x.as_ptr(),
                   ylen,
                   det48y.as_ptr(),
                   bdet.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   dablen,
                   dab.as_ptr(),
                   *pc.offset(0isize),
                   det24x.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   det24x.as_ptr(),
                   *pc.offset(0isize),
                   det48x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   dablen,
                   dab.as_ptr(),
                   *pc.offset(1isize),
                   det24y.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   det24y.as_ptr(),
                   *pc.offset(1isize),
                   det48y.as_mut_ptr()
               );
        clen = fast_expansion_sum_zeroelim(
                   xlen,
                   det48x.as_ptr(),
                   ylen,
                   det48y.as_ptr(),
                   cdet.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   abclen,
                   abc.as_ptr(),
                   *pd.offset(0isize),
                   det24x.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   det24x.as_ptr(),
                   -*pd.offset(0isize),
                   det48x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   abclen,
                   abc.as_ptr(),
                   *pd.offset(1isize),
                   det24y.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   det24y.as_ptr(),
                   -*pd.offset(1isize),
                   det48y.as_mut_ptr()
               );
        dlen = fast_expansion_sum_zeroelim(
                   xlen,
                   det48x.as_ptr(),
                   ylen,
                   det48y.as_ptr(),
                   ddet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        cdlen = fast_expansion_sum_zeroelim(
                    clen,
                    cdet.as_ptr(),
                    dlen,
                    ddet.as_ptr(),
                    cddet.as_mut_ptr()
                );
        deterlen = fast_expansion_sum_zeroelim(
                       ablen,
                       abdet.as_ptr(),
                       cdlen,
                       cddet.as_ptr(),
                       deter.as_mut_ptr()
                   );
        deter[(deterlen - 1i32) as (usize)]
    }

    
    pub unsafe fn incircleslow(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64
    ) -> f64 {
        let mut adx : f64;
        let mut bdx : f64;
        let mut cdx : f64;
        let mut ady : f64;
        let mut bdy : f64;
        let mut cdy : f64;
        let mut adxtail : f64;
        let mut bdxtail : f64;
        let mut cdxtail : f64;
        let mut adytail : f64;
        let mut bdytail : f64;
        let mut cdytail : f64;
        let mut negate : f64;
        let mut negatetail : f64;
        let mut axby7 : f64;
        let mut bxcy7 : f64;
        let mut axcy7 : f64;
        let mut bxay7 : f64;
        let mut cxby7 : f64;
        let mut cxay7 : f64;
        let mut axby : [f64; 8];
        let mut bxcy : [f64; 8];
        let mut axcy : [f64; 8];
        let mut bxay : [f64; 8];
        let mut cxby : [f64; 8];
        let mut cxay : [f64; 8];
        let mut temp16 : [f64; 16];
        let mut temp16len : i32;
        let mut detx : [f64; 32];
        let mut detxx : [f64; 64];
        let mut detxt : [f64; 32];
        let mut detxxt : [f64; 64];
        let mut detxtxt : [f64; 64];
        let mut xlen : i32;
        let mut xxlen : i32;
        let mut xtlen : i32;
        let mut xxtlen : i32;
        let mut xtxtlen : i32;
        let mut x1 : [f64; 128];
        let mut x2 : [f64; 192];
        let mut x1len : i32;
        let mut x2len : i32;
        let mut dety : [f64; 32];
        let mut detyy : [f64; 64];
        let mut detyt : [f64; 32];
        let mut detyyt : [f64; 64];
        let mut detytyt : [f64; 64];
        let mut ylen : i32;
        let mut yylen : i32;
        let mut ytlen : i32;
        let mut yytlen : i32;
        let mut ytytlen : i32;
        let mut y1 : [f64; 128];
        let mut y2 : [f64; 192];
        let mut y1len : i32;
        let mut y2len : i32;
        let mut adet : [f64; 384];
        let mut bdet : [f64; 384];
        let mut cdet : [f64; 384];
        let mut abdet : [f64; 768];
        let mut deter : [f64; 1152];
        let mut alen : i32;
        let mut blen : i32;
        let mut clen : i32;
        let mut ablen : i32;
        let mut deterlen : i32;
        let mut i : i32;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut a0hi : f64;
        let mut a0lo : f64;
        let mut a1hi : f64;
        let mut a1lo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _k : f64;
        let mut _l : f64;
        let mut _m : f64;
        let mut _n : f64;
        let mut _0 : f64;
        let mut _1 : f64;
        let mut _2 : f64;

        adx = uninitialized();
        bdx = uninitialized();
        cdx = uninitialized();
        ady = uninitialized();
        bdy = uninitialized();
        cdy = uninitialized();
        adxtail = uninitialized();
        bdxtail = uninitialized();
        cdxtail = uninitialized();
        adytail = uninitialized();
        bdytail = uninitialized();
        cdytail = uninitialized();
        negate = uninitialized();
        negatetail = uninitialized();
        axby7 = uninitialized();
        bxcy7 = uninitialized();
        axcy7 = uninitialized();
        bxay7 = uninitialized();
        cxby7 = uninitialized();
        cxay7 = uninitialized();
        axby = uninitialized();
        bxcy = uninitialized();
        axcy = uninitialized();
        bxay = uninitialized();
        cxby = uninitialized();
        cxay = uninitialized();
        temp16 = uninitialized();
        temp16len = uninitialized();
        detx = uninitialized();
        detxx = uninitialized();
        detxt = uninitialized();
        detxxt = uninitialized();
        detxtxt = uninitialized();
        xlen = uninitialized();
        xxlen = uninitialized();
        xtlen = uninitialized();
        xxtlen = uninitialized();
        xtxtlen = uninitialized();
        x1 = uninitialized();
        x2 = uninitialized();
        x1len = uninitialized();
        x2len = uninitialized();
        dety = uninitialized();
        detyy = uninitialized();
        detyt = uninitialized();
        detyyt = uninitialized();
        detytyt = uninitialized();
        ylen = uninitialized();
        yylen = uninitialized();
        ytlen = uninitialized();
        yytlen = uninitialized();
        ytytlen = uninitialized();
        y1 = uninitialized();
        y2 = uninitialized();
        y1len = uninitialized();
        y2len = uninitialized();
        adet = uninitialized();
        bdet = uninitialized();
        cdet = uninitialized();
        abdet = uninitialized();
        deter = uninitialized();
        alen = uninitialized();
        blen = uninitialized();
        clen = uninitialized();
        ablen = uninitialized();
        deterlen = uninitialized();
        i = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        a0hi = uninitialized();
        a0lo = uninitialized();
        a1hi = uninitialized();
        a1lo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _k = uninitialized();
        _l = uninitialized();
        _m = uninitialized();
        _n = uninitialized();
        _0 = uninitialized();
        _1 = uninitialized();
        _2 = uninitialized();
        adx = *pa.offset(0isize) - *pd.offset(0isize);
        bvirt = *pa.offset(0isize) - adx;
        avirt = adx + bvirt;
        bround = bvirt - *pd.offset(0isize);
        around = *pa.offset(0isize) - avirt;
        adxtail = around + bround;
        ady = *pa.offset(1isize) - *pd.offset(1isize);
        bvirt = *pa.offset(1isize) - ady;
        avirt = ady + bvirt;
        bround = bvirt - *pd.offset(1isize);
        around = *pa.offset(1isize) - avirt;
        adytail = around + bround;
        bdx = *pb.offset(0isize) - *pd.offset(0isize);
        bvirt = *pb.offset(0isize) - bdx;
        avirt = bdx + bvirt;
        bround = bvirt - *pd.offset(0isize);
        around = *pb.offset(0isize) - avirt;
        bdxtail = around + bround;
        bdy = *pb.offset(1isize) - *pd.offset(1isize);
        bvirt = *pb.offset(1isize) - bdy;
        avirt = bdy + bvirt;
        bround = bvirt - *pd.offset(1isize);
        around = *pb.offset(1isize) - avirt;
        bdytail = around + bround;
        cdx = *pc.offset(0isize) - *pd.offset(0isize);
        bvirt = *pc.offset(0isize) - cdx;
        avirt = cdx + bvirt;
        bround = bvirt - *pd.offset(0isize);
        around = *pc.offset(0isize) - avirt;
        cdxtail = around + bround;
        cdy = *pc.offset(1isize) - *pd.offset(1isize);
        bvirt = *pc.offset(1isize) - cdy;
        avirt = cdy + bvirt;
        bround = bvirt - *pd.offset(1isize);
        around = *pc.offset(1isize) - avirt;
        cdytail = around + bround;
        c = self.splitter * adxtail;
        abig = c - adxtail;
        a0hi = c - abig;
        a0lo = adxtail - a0hi;
        c = self.splitter * bdytail;
        abig = c - bdytail;
        bhi = c - abig;
        blo = bdytail - bhi;
        _i = adxtail * bdytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        axby[0usize] = a0lo * blo - err3;
        c = self.splitter * adx;
        abig = c - adx;
        a1hi = c - abig;
        a1lo = adx - a1hi;
        _j = adx * bdytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * bdy;
        abig = c - bdy;
        bhi = c - abig;
        blo = bdy - bhi;
        _i = adxtail * bdy;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = adx * bdy;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        axby[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        axby[5usize] = around + bround;
        axby7 = _m + _k;
        bvirt = axby7 - _m;
        avirt = axby7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        axby[6usize] = around + bround;
        axby[7usize] = axby7;
        negate = -ady;
        negatetail = -adytail;
        c = self.splitter * bdxtail;
        abig = c - bdxtail;
        a0hi = c - abig;
        a0lo = bdxtail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = bdxtail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        bxay[0usize] = a0lo * blo - err3;
        c = self.splitter * bdx;
        abig = c - bdx;
        a1hi = c - abig;
        a1lo = bdx - a1hi;
        _j = bdx * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = bdxtail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = bdx * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        bxay[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        bxay[5usize] = around + bround;
        bxay7 = _m + _k;
        bvirt = bxay7 - _m;
        avirt = bxay7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        bxay[6usize] = around + bround;
        bxay[7usize] = bxay7;
        c = self.splitter * bdxtail;
        abig = c - bdxtail;
        a0hi = c - abig;
        a0lo = bdxtail - a0hi;
        c = self.splitter * cdytail;
        abig = c - cdytail;
        bhi = c - abig;
        blo = cdytail - bhi;
        _i = bdxtail * cdytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        bxcy[0usize] = a0lo * blo - err3;
        c = self.splitter * bdx;
        abig = c - bdx;
        a1hi = c - abig;
        a1lo = bdx - a1hi;
        _j = bdx * cdytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * cdy;
        abig = c - cdy;
        bhi = c - abig;
        blo = cdy - bhi;
        _i = bdxtail * cdy;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = bdx * cdy;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        bxcy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        bxcy[5usize] = around + bround;
        bxcy7 = _m + _k;
        bvirt = bxcy7 - _m;
        avirt = bxcy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        bxcy[6usize] = around + bround;
        bxcy[7usize] = bxcy7;
        negate = -bdy;
        negatetail = -bdytail;
        c = self.splitter * cdxtail;
        abig = c - cdxtail;
        a0hi = c - abig;
        a0lo = cdxtail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = cdxtail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        cxby[0usize] = a0lo * blo - err3;
        c = self.splitter * cdx;
        abig = c - cdx;
        a1hi = c - abig;
        a1lo = cdx - a1hi;
        _j = cdx * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = cdxtail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = cdx * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        cxby[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        cxby[5usize] = around + bround;
        cxby7 = _m + _k;
        bvirt = cxby7 - _m;
        avirt = cxby7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        cxby[6usize] = around + bround;
        cxby[7usize] = cxby7;
        c = self.splitter * cdxtail;
        abig = c - cdxtail;
        a0hi = c - abig;
        a0lo = cdxtail - a0hi;
        c = self.splitter * adytail;
        abig = c - adytail;
        bhi = c - abig;
        blo = adytail - bhi;
        _i = cdxtail * adytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        cxay[0usize] = a0lo * blo - err3;
        c = self.splitter * cdx;
        abig = c - cdx;
        a1hi = c - abig;
        a1lo = cdx - a1hi;
        _j = cdx * adytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * ady;
        abig = c - ady;
        bhi = c - abig;
        blo = ady - bhi;
        _i = cdxtail * ady;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = cdx * ady;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        cxay[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        cxay[5usize] = around + bround;
        cxay7 = _m + _k;
        bvirt = cxay7 - _m;
        avirt = cxay7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        cxay[6usize] = around + bround;
        cxay[7usize] = cxay7;
        negate = -cdy;
        negatetail = -cdytail;
        c = self.splitter * adxtail;
        abig = c - adxtail;
        a0hi = c - abig;
        a0lo = adxtail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = adxtail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        axcy[0usize] = a0lo * blo - err3;
        c = self.splitter * adx;
        abig = c - adx;
        a1hi = c - abig;
        a1lo = adx - a1hi;
        _j = adx * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = adxtail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = adx * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        axcy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        axcy[5usize] = around + bround;
        axcy7 = _m + _k;
        bvirt = axcy7 - _m;
        avirt = axcy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        axcy[6usize] = around + bround;
        axcy[7usize] = axcy7;
        temp16len = fast_expansion_sum_zeroelim(
                        8i32,
                        bxcy.as_ptr(),
                        8i32,
                        cxby.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        xlen = self.scale_expansion_zeroelim(
                   temp16len,
                   temp16.as_ptr(),
                   adx,
                   detx.as_mut_ptr()
               );
        xxlen = self.scale_expansion_zeroelim(
                    xlen,
                    detx.as_ptr(),
                    adx,
                    detxx.as_mut_ptr()
                );
        xtlen = self.scale_expansion_zeroelim(
                    temp16len,
                    temp16.as_ptr(),
                    adxtail,
                    detxt.as_mut_ptr()
                );
        xxtlen = self.scale_expansion_zeroelim(
                     xtlen,
                     detxt.as_ptr(),
                     adx,
                     detxxt.as_mut_ptr()
                 );
        i = 0i32;
        'loop1: loop {
            if !(i < xxtlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detxxt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        xtxtlen = self.scale_expansion_zeroelim(
                      xtlen,
                      detxt.as_ptr(),
                      adxtail,
                      detxtxt.as_mut_ptr()
                  );
        x1len = fast_expansion_sum_zeroelim(
                    xxlen,
                    detxx.as_ptr(),
                    xxtlen,
                    detxxt.as_ptr(),
                    x1.as_mut_ptr()
                );
        x2len = fast_expansion_sum_zeroelim(
                    x1len,
                    x1.as_ptr(),
                    xtxtlen,
                    detxtxt.as_ptr(),
                    x2.as_mut_ptr()
                );
        ylen = self.scale_expansion_zeroelim(
                   temp16len,
                   temp16.as_ptr(),
                   ady,
                   dety.as_mut_ptr()
               );
        yylen = self.scale_expansion_zeroelim(
                    ylen,
                    dety.as_ptr(),
                    ady,
                    detyy.as_mut_ptr()
                );
        ytlen = self.scale_expansion_zeroelim(
                    temp16len,
                    temp16.as_ptr(),
                    adytail,
                    detyt.as_mut_ptr()
                );
        yytlen = self.scale_expansion_zeroelim(
                     ytlen,
                     detyt.as_ptr(),
                     ady,
                     detyyt.as_mut_ptr()
                 );
        i = 0i32;
        'loop3: loop {
            if !(i < yytlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detyyt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ytytlen = self.scale_expansion_zeroelim(
                      ytlen,
                      detyt.as_ptr(),
                      adytail,
                      detytyt.as_mut_ptr()
                  );
        y1len = fast_expansion_sum_zeroelim(
                    yylen,
                    detyy.as_ptr(),
                    yytlen,
                    detyyt.as_ptr(),
                    y1.as_mut_ptr()
                );
        y2len = fast_expansion_sum_zeroelim(
                    y1len,
                    y1.as_ptr(),
                    ytytlen,
                    detytyt.as_ptr(),
                    y2.as_mut_ptr()
                );
        alen = fast_expansion_sum_zeroelim(
                   x2len,
                   x2.as_ptr(),
                   y2len,
                   y2.as_ptr(),
                   adet.as_mut_ptr()
               );
        temp16len = fast_expansion_sum_zeroelim(
                        8i32,
                        cxay.as_ptr(),
                        8i32,
                        axcy.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        xlen = self.scale_expansion_zeroelim(
                   temp16len,
                   temp16.as_ptr(),
                   bdx,
                   detx.as_mut_ptr()
               );
        xxlen = self.scale_expansion_zeroelim(
                    xlen,
                    detx.as_ptr(),
                    bdx,
                    detxx.as_mut_ptr()
                );
        xtlen = self.scale_expansion_zeroelim(
                    temp16len,
                    temp16.as_ptr(),
                    bdxtail,
                    detxt.as_mut_ptr()
                );
        xxtlen = self.scale_expansion_zeroelim(
                     xtlen,
                     detxt.as_ptr(),
                     bdx,
                     detxxt.as_mut_ptr()
                 );
        i = 0i32;
        'loop5: loop {
            if !(i < xxtlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detxxt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        xtxtlen = self.scale_expansion_zeroelim(
                      xtlen,
                      detxt.as_ptr(),
                      bdxtail,
                      detxtxt.as_mut_ptr()
                  );
        x1len = fast_expansion_sum_zeroelim(
                    xxlen,
                    detxx.as_ptr(),
                    xxtlen,
                    detxxt.as_ptr(),
                    x1.as_mut_ptr()
                );
        x2len = fast_expansion_sum_zeroelim(
                    x1len,
                    x1.as_ptr(),
                    xtxtlen,
                    detxtxt.as_ptr(),
                    x2.as_mut_ptr()
                );
        ylen = self.scale_expansion_zeroelim(
                   temp16len,
                   temp16.as_ptr(),
                   bdy,
                   dety.as_mut_ptr()
               );
        yylen = self.scale_expansion_zeroelim(
                    ylen,
                    dety.as_ptr(),
                    bdy,
                    detyy.as_mut_ptr()
                );
        ytlen = self.scale_expansion_zeroelim(
                    temp16len,
                    temp16.as_ptr(),
                    bdytail,
                    detyt.as_mut_ptr()
                );
        yytlen = self.scale_expansion_zeroelim(
                     ytlen,
                     detyt.as_ptr(),
                     bdy,
                     detyyt.as_mut_ptr()
                 );
        i = 0i32;
        'loop7: loop {
            if !(i < yytlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detyyt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ytytlen = self.scale_expansion_zeroelim(
                      ytlen,
                      detyt.as_ptr(),
                      bdytail,
                      detytyt.as_mut_ptr()
                  );
        y1len = fast_expansion_sum_zeroelim(
                    yylen,
                    detyy.as_ptr(),
                    yytlen,
                    detyyt.as_ptr(),
                    y1.as_mut_ptr()
                );
        y2len = fast_expansion_sum_zeroelim(
                    y1len,
                    y1.as_ptr(),
                    ytytlen,
                    detytyt.as_ptr(),
                    y2.as_mut_ptr()
                );
        blen = fast_expansion_sum_zeroelim(
                   x2len,
                   x2.as_ptr(),
                   y2len,
                   y2.as_ptr(),
                   bdet.as_mut_ptr()
               );
        temp16len = fast_expansion_sum_zeroelim(
                        8i32,
                        axby.as_ptr(),
                        8i32,
                        bxay.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        xlen = self.scale_expansion_zeroelim(
                   temp16len,
                   temp16.as_ptr(),
                   cdx,
                   detx.as_mut_ptr()
               );
        xxlen = self.scale_expansion_zeroelim(
                    xlen,
                    detx.as_ptr(),
                    cdx,
                    detxx.as_mut_ptr()
                );
        xtlen = self.scale_expansion_zeroelim(
                    temp16len,
                    temp16.as_ptr(),
                    cdxtail,
                    detxt.as_mut_ptr()
                );
        xxtlen = self.scale_expansion_zeroelim(
                     xtlen,
                     detxt.as_ptr(),
                     cdx,
                     detxxt.as_mut_ptr()
                 );
        i = 0i32;
        'loop9: loop {
            if !(i < xxtlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detxxt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        xtxtlen = self.scale_expansion_zeroelim(
                      xtlen,
                      detxt.as_ptr(),
                      cdxtail,
                      detxtxt.as_mut_ptr()
                  );
        x1len = fast_expansion_sum_zeroelim(
                    xxlen,
                    detxx.as_ptr(),
                    xxtlen,
                    detxxt.as_ptr(),
                    x1.as_mut_ptr()
                );
        x2len = fast_expansion_sum_zeroelim(
                    x1len,
                    x1.as_ptr(),
                    xtxtlen,
                    detxtxt.as_ptr(),
                    x2.as_mut_ptr()
                );
        ylen = self.scale_expansion_zeroelim(
                   temp16len,
                   temp16.as_ptr(),
                   cdy,
                   dety.as_mut_ptr()
               );
        yylen = self.scale_expansion_zeroelim(
                    ylen,
                    dety.as_ptr(),
                    cdy,
                    detyy.as_mut_ptr()
                );
        ytlen = self.scale_expansion_zeroelim(
                    temp16len,
                    temp16.as_ptr(),
                    cdytail,
                    detyt.as_mut_ptr()
                );
        yytlen = self.scale_expansion_zeroelim(
                     ytlen,
                     detyt.as_ptr(),
                     cdy,
                     detyyt.as_mut_ptr()
                 );
        i = 0i32;
        'loop11: loop {
            if !(i < yytlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detyyt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ytytlen = self.scale_expansion_zeroelim(
                      ytlen,
                      detyt.as_ptr(),
                      cdytail,
                      detytyt.as_mut_ptr()
                  );
        y1len = fast_expansion_sum_zeroelim(
                    yylen,
                    detyy.as_ptr(),
                    yytlen,
                    detyyt.as_ptr(),
                    y1.as_mut_ptr()
                );
        y2len = fast_expansion_sum_zeroelim(
                    y1len,
                    y1.as_ptr(),
                    ytytlen,
                    detytyt.as_ptr(),
                    y2.as_mut_ptr()
                );
        clen = fast_expansion_sum_zeroelim(
                   x2len,
                   x2.as_ptr(),
                   y2len,
                   y2.as_ptr(),
                   cdet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        deterlen = fast_expansion_sum_zeroelim(
                       ablen,
                       abdet.as_ptr(),
                       clen,
                       cdet.as_ptr(),
                       deter.as_mut_ptr()
                   );
        deter[(deterlen - 1i32) as (usize)]
    }

    
    pub unsafe fn incircleadapt(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64,
        mut permanent : f64
    ) -> f64 {
        let mut adx : f64;
        let mut bdx : f64;
        let mut cdx : f64;
        let mut ady : f64;
        let mut bdy : f64;
        let mut cdy : f64;
        let mut det : f64;
        let mut errbound : f64;
        let mut bdxcdy1 : f64;
        let mut cdxbdy1 : f64;
        let mut cdxady1 : f64;
        let mut adxcdy1 : f64;
        let mut adxbdy1 : f64;
        let mut bdxady1 : f64;
        let mut bdxcdy0 : f64;
        let mut cdxbdy0 : f64;
        let mut cdxady0 : f64;
        let mut adxcdy0 : f64;
        let mut adxbdy0 : f64;
        let mut bdxady0 : f64;
        let mut bc : [f64; 4];
        let mut ca : [f64; 4];
        let mut ab : [f64; 4];
        let mut bc3 : f64;
        let mut ca3 : f64;
        let mut ab3 : f64;
        let mut axbc : [f64; 8];
        let mut axxbc : [f64; 16];
        let mut aybc : [f64; 8];
        let mut ayybc : [f64; 16];
        let mut adet : [f64; 32];
        let mut axbclen : i32;
        let mut axxbclen : i32;
        let mut aybclen : i32;
        let mut ayybclen : i32;
        let mut alen : i32;
        let mut bxca : [f64; 8];
        let mut bxxca : [f64; 16];
        let mut byca : [f64; 8];
        let mut byyca : [f64; 16];
        let mut bdet : [f64; 32];
        let mut bxcalen : i32;
        let mut bxxcalen : i32;
        let mut bycalen : i32;
        let mut byycalen : i32;
        let mut blen : i32;
        let mut cxab : [f64; 8];
        let mut cxxab : [f64; 16];
        let mut cyab : [f64; 8];
        let mut cyyab : [f64; 16];
        let mut cdet : [f64; 32];
        let mut cxablen : i32;
        let mut cxxablen : i32;
        let mut cyablen : i32;
        let mut cyyablen : i32;
        let mut clen : i32;
        let mut abdet : [f64; 64];
        let mut ablen : i32;
        let mut fin1 : [f64; 1152];
        let mut fin2 : [f64; 1152];
        let mut finnow : *mut f64;
        let mut finother : *mut f64;
        let mut finswap : *mut f64;
        let mut finlength : i32;
        let mut adxtail : f64;
        let mut bdxtail : f64;
        let mut cdxtail : f64;
        let mut adytail : f64;
        let mut bdytail : f64;
        let mut cdytail : f64;
        let mut adxadx1 : f64;
        let mut adyady1 : f64;
        let mut bdxbdx1 : f64;
        let mut bdybdy1 : f64;
        let mut cdxcdx1 : f64;
        let mut cdycdy1 : f64;
        let mut adxadx0 : f64;
        let mut adyady0 : f64;
        let mut bdxbdx0 : f64;
        let mut bdybdy0 : f64;
        let mut cdxcdx0 : f64;
        let mut cdycdy0 : f64;
        let mut aa : [f64; 4];
        let mut bb : [f64; 4];
        let mut cc : [f64; 4];
        let mut aa3 : f64;
        let mut bb3 : f64;
        let mut cc3 : f64;
        let mut ti1 : f64;
        let mut tj1 : f64;
        let mut ti0 : f64;
        let mut tj0 : f64;
        let mut u : [f64; 4];
        let mut v : [f64; 4];
        let mut u3 : f64;
        let mut v3 : f64;
        let mut temp8 : [f64; 8];
        let mut temp16a : [f64; 16];
        let mut temp16b : [f64; 16];
        let mut temp16c : [f64; 16];
        let mut temp32a : [f64; 32];
        let mut temp32b : [f64; 32];
        let mut temp48 : [f64; 48];
        let mut temp64 : [f64; 64];
        let mut temp8len : i32;
        let mut temp16alen : i32;
        let mut temp16blen : i32;
        let mut temp16clen : i32;
        let mut temp32alen : i32;
        let mut temp32blen : i32;
        let mut temp48len : i32;
        let mut temp64len : i32;
        let mut axtbb : [f64; 8];
        let mut axtcc : [f64; 8];
        let mut aytbb : [f64; 8];
        let mut aytcc : [f64; 8];
        let mut axtbblen : i32;
        let mut axtcclen : i32;
        let mut aytbblen : i32;
        let mut aytcclen : i32;
        let mut bxtaa : [f64; 8];
        let mut bxtcc : [f64; 8];
        let mut bytaa : [f64; 8];
        let mut bytcc : [f64; 8];
        let mut bxtaalen : i32;
        let mut bxtcclen : i32;
        let mut bytaalen : i32;
        let mut bytcclen : i32;
        let mut cxtaa : [f64; 8];
        let mut cxtbb : [f64; 8];
        let mut cytaa : [f64; 8];
        let mut cytbb : [f64; 8];
        let mut cxtaalen : i32;
        let mut cxtbblen : i32;
        let mut cytaalen : i32;
        let mut cytbblen : i32;
        let mut axtbc : [f64; 8];
        let mut aytbc : [f64; 8];
        let mut bxtca : [f64; 8];
        let mut bytca : [f64; 8];
        let mut cxtab : [f64; 8];
        let mut cytab : [f64; 8];
        let mut axtbclen : i32;
        let mut aytbclen : i32;
        let mut bxtcalen : i32;
        let mut bytcalen : i32;
        let mut cxtablen : i32;
        let mut cytablen : i32;
        let mut axtbct : [f64; 16];
        let mut aytbct : [f64; 16];
        let mut bxtcat : [f64; 16];
        let mut bytcat : [f64; 16];
        let mut cxtabt : [f64; 16];
        let mut cytabt : [f64; 16];
        let mut axtbctlen : i32;
        let mut aytbctlen : i32;
        let mut bxtcatlen : i32;
        let mut bytcatlen : i32;
        let mut cxtabtlen : i32;
        let mut cytabtlen : i32;
        let mut axtbctt : [f64; 8];
        let mut aytbctt : [f64; 8];
        let mut bxtcatt : [f64; 8];
        let mut bytcatt : [f64; 8];
        let mut cxtabtt : [f64; 8];
        let mut cytabtt : [f64; 8];
        let mut axtbcttlen : i32;
        let mut aytbcttlen : i32;
        let mut bxtcattlen : i32;
        let mut bytcattlen : i32;
        let mut cxtabttlen : i32;
        let mut cytabttlen : i32;
        let mut abt : [f64; 8];
        let mut bct : [f64; 8];
        let mut cat : [f64; 8];
        let mut abtlen : i32;
        let mut bctlen : i32;
        let mut catlen : i32;
        let mut abtt : [f64; 4];
        let mut bctt : [f64; 4];
        let mut catt : [f64; 4];
        let mut abttlen : i32;
        let mut bcttlen : i32;
        let mut cattlen : i32;
        let mut abtt3 : f64;
        let mut bctt3 : f64;
        let mut catt3 : f64;
        let mut negate : f64;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _0 : f64;

        adx = uninitialized();
        bdx = uninitialized();
        cdx = uninitialized();
        ady = uninitialized();
        bdy = uninitialized();
        cdy = uninitialized();
        det = uninitialized();
        errbound = uninitialized();
        bdxcdy1 = uninitialized();
        cdxbdy1 = uninitialized();
        cdxady1 = uninitialized();
        adxcdy1 = uninitialized();
        adxbdy1 = uninitialized();
        bdxady1 = uninitialized();
        bdxcdy0 = uninitialized();
        cdxbdy0 = uninitialized();
        cdxady0 = uninitialized();
        adxcdy0 = uninitialized();
        adxbdy0 = uninitialized();
        bdxady0 = uninitialized();
        bc = uninitialized();
        ca = uninitialized();
        ab = uninitialized();
        bc3 = uninitialized();
        ca3 = uninitialized();
        ab3 = uninitialized();
        axbc = uninitialized();
        axxbc = uninitialized();
        aybc = uninitialized();
        ayybc = uninitialized();
        adet = uninitialized();
        axbclen = uninitialized();
        axxbclen = uninitialized();
        aybclen = uninitialized();
        ayybclen = uninitialized();
        alen = uninitialized();
        bxca = uninitialized();
        bxxca = uninitialized();
        byca = uninitialized();
        byyca = uninitialized();
        bdet = uninitialized();
        bxcalen = uninitialized();
        bxxcalen = uninitialized();
        bycalen = uninitialized();
        byycalen = uninitialized();
        blen = uninitialized();
        cxab = uninitialized();
        cxxab = uninitialized();
        cyab = uninitialized();
        cyyab = uninitialized();
        cdet = uninitialized();
        cxablen = uninitialized();
        cxxablen = uninitialized();
        cyablen = uninitialized();
        cyyablen = uninitialized();
        clen = uninitialized();
        abdet = uninitialized();
        ablen = uninitialized();
        fin1 = uninitialized();
        fin2 = uninitialized();
        finnow = uninitialized();
        finother = uninitialized();
        finswap = uninitialized();
        finlength = uninitialized();
        adxtail = uninitialized();
        bdxtail = uninitialized();
        cdxtail = uninitialized();
        adytail = uninitialized();
        bdytail = uninitialized();
        cdytail = uninitialized();
        adxadx1 = uninitialized();
        adyady1 = uninitialized();
        bdxbdx1 = uninitialized();
        bdybdy1 = uninitialized();
        cdxcdx1 = uninitialized();
        cdycdy1 = uninitialized();
        adxadx0 = uninitialized();
        adyady0 = uninitialized();
        bdxbdx0 = uninitialized();
        bdybdy0 = uninitialized();
        cdxcdx0 = uninitialized();
        cdycdy0 = uninitialized();
        aa = uninitialized();
        bb = uninitialized();
        cc = uninitialized();
        aa3 = uninitialized();
        bb3 = uninitialized();
        cc3 = uninitialized();
        ti1 = uninitialized();
        tj1 = uninitialized();
        ti0 = uninitialized();
        tj0 = uninitialized();
        u = uninitialized();
        v = uninitialized();
        u3 = uninitialized();
        v3 = uninitialized();
        temp8 = uninitialized();
        temp16a = uninitialized();
        temp16b = uninitialized();
        temp16c = uninitialized();
        temp32a = uninitialized();
        temp32b = uninitialized();
        temp48 = uninitialized();
        temp64 = uninitialized();
        temp8len = uninitialized();
        temp16alen = uninitialized();
        temp16blen = uninitialized();
        temp16clen = uninitialized();
        temp32alen = uninitialized();
        temp32blen = uninitialized();
        temp48len = uninitialized();
        temp64len = uninitialized();
        axtbb = uninitialized();
        axtcc = uninitialized();
        aytbb = uninitialized();
        aytcc = uninitialized();
        axtbblen = uninitialized();
        axtcclen = uninitialized();
        aytbblen = uninitialized();
        aytcclen = uninitialized();
        bxtaa = uninitialized();
        bxtcc = uninitialized();
        bytaa = uninitialized();
        bytcc = uninitialized();
        bxtaalen = uninitialized();
        bxtcclen = uninitialized();
        bytaalen = uninitialized();
        bytcclen = uninitialized();
        cxtaa = uninitialized();
        cxtbb = uninitialized();
        cytaa = uninitialized();
        cytbb = uninitialized();
        cxtaalen = uninitialized();
        cxtbblen = uninitialized();
        cytaalen = uninitialized();
        cytbblen = uninitialized();
        axtbc = uninitialized();
        aytbc = uninitialized();
        bxtca = uninitialized();
        bytca = uninitialized();
        cxtab = uninitialized();
        cytab = uninitialized();
        axtbclen = uninitialized();
        aytbclen = uninitialized();
        bxtcalen = uninitialized();
        bytcalen = uninitialized();
        cxtablen = uninitialized();
        cytablen = uninitialized();
        axtbct = uninitialized();
        aytbct = uninitialized();
        bxtcat = uninitialized();
        bytcat = uninitialized();
        cxtabt = uninitialized();
        cytabt = uninitialized();
        axtbctlen = uninitialized();
        aytbctlen = uninitialized();
        bxtcatlen = uninitialized();
        bytcatlen = uninitialized();
        cxtabtlen = uninitialized();
        cytabtlen = uninitialized();
        axtbctt = uninitialized();
        aytbctt = uninitialized();
        bxtcatt = uninitialized();
        bytcatt = uninitialized();
        cxtabtt = uninitialized();
        cytabtt = uninitialized();
        axtbcttlen = uninitialized();
        aytbcttlen = uninitialized();
        bxtcattlen = uninitialized();
        bytcattlen = uninitialized();
        cxtabttlen = uninitialized();
        cytabttlen = uninitialized();
        abt = uninitialized();
        bct = uninitialized();
        cat = uninitialized();
        abtlen = uninitialized();
        bctlen = uninitialized();
        catlen = uninitialized();
        abtt = uninitialized();
        bctt = uninitialized();
        catt = uninitialized();
        abttlen = uninitialized();
        bcttlen = uninitialized();
        cattlen = uninitialized();
        abtt3 = uninitialized();
        bctt3 = uninitialized();
        catt3 = uninitialized();
        negate = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        ahi = uninitialized();
        alo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _0 = uninitialized();
        adx = *pa.offset(0isize) - *pd.offset(0isize);
        bdx = *pb.offset(0isize) - *pd.offset(0isize);
        cdx = *pc.offset(0isize) - *pd.offset(0isize);
        ady = *pa.offset(1isize) - *pd.offset(1isize);
        bdy = *pb.offset(1isize) - *pd.offset(1isize);
        cdy = *pc.offset(1isize) - *pd.offset(1isize);
        bdxcdy1 = bdx * cdy;
        c = self.splitter * bdx;
        abig = c - bdx;
        ahi = c - abig;
        alo = bdx - ahi;
        c = self.splitter * cdy;
        abig = c - cdy;
        bhi = c - abig;
        blo = cdy - bhi;
        err1 = bdxcdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bdxcdy0 = alo * blo - err3;
        cdxbdy1 = cdx * bdy;
        c = self.splitter * cdx;
        abig = c - cdx;
        ahi = c - abig;
        alo = cdx - ahi;
        c = self.splitter * bdy;
        abig = c - bdy;
        bhi = c - abig;
        blo = bdy - bhi;
        err1 = cdxbdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cdxbdy0 = alo * blo - err3;
        _i = bdxcdy0 - cdxbdy0;
        bvirt = bdxcdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cdxbdy0;
        around = bdxcdy0 - avirt;
        bc[0usize] = around + bround;
        _j = bdxcdy1 + _i;
        bvirt = _j - bdxcdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bdxcdy1 - avirt;
        _0 = around + bround;
        _i = _0 - cdxbdy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cdxbdy1;
        around = _0 - avirt;
        bc[1usize] = around + bround;
        bc3 = _j + _i;
        bvirt = bc3 - _j;
        avirt = bc3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bc[2usize] = around + bround;
        bc[3usize] = bc3;
        axbclen = self.scale_expansion_zeroelim(
                      4i32,
                      bc.as_ptr(),
                      adx,
                      axbc.as_mut_ptr()
                  );
        axxbclen = self.scale_expansion_zeroelim(
                       axbclen,
                       axbc.as_ptr(),
                       adx,
                       axxbc.as_mut_ptr()
                   );
        aybclen = self.scale_expansion_zeroelim(
                      4i32,
                      bc.as_ptr(),
                      ady,
                      aybc.as_mut_ptr()
                  );
        ayybclen = self.scale_expansion_zeroelim(
                       aybclen,
                       aybc.as_ptr(),
                       ady,
                       ayybc.as_mut_ptr()
                   );
        alen = fast_expansion_sum_zeroelim(
                   axxbclen,
                   axxbc.as_ptr(),
                   ayybclen,
                   ayybc.as_ptr(),
                   adet.as_mut_ptr()
               );
        cdxady1 = cdx * ady;
        c = self.splitter * cdx;
        abig = c - cdx;
        ahi = c - abig;
        alo = cdx - ahi;
        c = self.splitter * ady;
        abig = c - ady;
        bhi = c - abig;
        blo = ady - bhi;
        err1 = cdxady1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cdxady0 = alo * blo - err3;
        adxcdy1 = adx * cdy;
        c = self.splitter * adx;
        abig = c - adx;
        ahi = c - abig;
        alo = adx - ahi;
        c = self.splitter * cdy;
        abig = c - cdy;
        bhi = c - abig;
        blo = cdy - bhi;
        err1 = adxcdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        adxcdy0 = alo * blo - err3;
        _i = cdxady0 - adxcdy0;
        bvirt = cdxady0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - adxcdy0;
        around = cdxady0 - avirt;
        ca[0usize] = around + bround;
        _j = cdxady1 + _i;
        bvirt = _j - cdxady1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = cdxady1 - avirt;
        _0 = around + bround;
        _i = _0 - adxcdy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - adxcdy1;
        around = _0 - avirt;
        ca[1usize] = around + bround;
        ca3 = _j + _i;
        bvirt = ca3 - _j;
        avirt = ca3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ca[2usize] = around + bround;
        ca[3usize] = ca3;
        bxcalen = self.scale_expansion_zeroelim(
                      4i32,
                      ca.as_ptr(),
                      bdx,
                      bxca.as_mut_ptr()
                  );
        bxxcalen = self.scale_expansion_zeroelim(
                       bxcalen,
                       bxca.as_ptr(),
                       bdx,
                       bxxca.as_mut_ptr()
                   );
        bycalen = self.scale_expansion_zeroelim(
                      4i32,
                      ca.as_ptr(),
                      bdy,
                      byca.as_mut_ptr()
                  );
        byycalen = self.scale_expansion_zeroelim(
                       bycalen,
                       byca.as_ptr(),
                       bdy,
                       byyca.as_mut_ptr()
                   );
        blen = fast_expansion_sum_zeroelim(
                   bxxcalen,
                   bxxca.as_ptr(),
                   byycalen,
                   byyca.as_ptr(),
                   bdet.as_mut_ptr()
               );
        adxbdy1 = adx * bdy;
        c = self.splitter * adx;
        abig = c - adx;
        ahi = c - abig;
        alo = adx - ahi;
        c = self.splitter * bdy;
        abig = c - bdy;
        bhi = c - abig;
        blo = bdy - bhi;
        err1 = adxbdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        adxbdy0 = alo * blo - err3;
        bdxady1 = bdx * ady;
        c = self.splitter * bdx;
        abig = c - bdx;
        ahi = c - abig;
        alo = bdx - ahi;
        c = self.splitter * ady;
        abig = c - ady;
        bhi = c - abig;
        blo = ady - bhi;
        err1 = bdxady1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bdxady0 = alo * blo - err3;
        _i = adxbdy0 - bdxady0;
        bvirt = adxbdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bdxady0;
        around = adxbdy0 - avirt;
        ab[0usize] = around + bround;
        _j = adxbdy1 + _i;
        bvirt = _j - adxbdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = adxbdy1 - avirt;
        _0 = around + bround;
        _i = _0 - bdxady1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bdxady1;
        around = _0 - avirt;
        ab[1usize] = around + bround;
        ab3 = _j + _i;
        bvirt = ab3 - _j;
        avirt = ab3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ab[2usize] = around + bround;
        ab[3usize] = ab3;
        cxablen = self.scale_expansion_zeroelim(
                      4i32,
                      ab.as_ptr(),
                      cdx,
                      cxab.as_mut_ptr()
                  );
        cxxablen = self.scale_expansion_zeroelim(
                       cxablen,
                       cxab.as_ptr(),
                       cdx,
                       cxxab.as_mut_ptr()
                   );
        cyablen = self.scale_expansion_zeroelim(
                      4i32,
                      ab.as_ptr(),
                      cdy,
                      cyab.as_mut_ptr()
                  );
        cyyablen = self.scale_expansion_zeroelim(
                       cyablen,
                       cyab.as_ptr(),
                       cdy,
                       cyyab.as_mut_ptr()
                   );
        clen = fast_expansion_sum_zeroelim(
                   cxxablen,
                   cxxab.as_ptr(),
                   cyyablen,
                   cyyab.as_ptr(),
                   cdet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        finlength = fast_expansion_sum_zeroelim(
                        ablen,
                        abdet.as_ptr(),
                        clen,
                        cdet.as_ptr(),
                        fin1.as_mut_ptr()
                    );
        det = estimate(finlength,fin1.as_ptr());
        errbound = self.iccerrboundB * permanent;
        if det >= errbound || -det >= errbound {
            det
        } else {
            bvirt = *pa.offset(0isize) - adx;
            avirt = adx + bvirt;
            bround = bvirt - *pd.offset(0isize);
            around = *pa.offset(0isize) - avirt;
            adxtail = around + bround;
            bvirt = *pa.offset(1isize) - ady;
            avirt = ady + bvirt;
            bround = bvirt - *pd.offset(1isize);
            around = *pa.offset(1isize) - avirt;
            adytail = around + bround;
            bvirt = *pb.offset(0isize) - bdx;
            avirt = bdx + bvirt;
            bround = bvirt - *pd.offset(0isize);
            around = *pb.offset(0isize) - avirt;
            bdxtail = around + bround;
            bvirt = *pb.offset(1isize) - bdy;
            avirt = bdy + bvirt;
            bround = bvirt - *pd.offset(1isize);
            around = *pb.offset(1isize) - avirt;
            bdytail = around + bround;
            bvirt = *pc.offset(0isize) - cdx;
            avirt = cdx + bvirt;
            bround = bvirt - *pd.offset(0isize);
            around = *pc.offset(0isize) - avirt;
            cdxtail = around + bround;
            bvirt = *pc.offset(1isize) - cdy;
            avirt = cdy + bvirt;
            bround = bvirt - *pd.offset(1isize);
            around = *pc.offset(1isize) - avirt;
            cdytail = around + bround;
            (if adxtail == 0.0f64 && (bdxtail == 0.0f64) && (cdxtail == 0.0f64) && (adytail == 0.0f64) && (bdytail == 0.0f64) && (cdytail == 0.0f64) {
                 det
             } else {
                 errbound = self.iccerrboundC * permanent + self.resulterrbound * Absolute(
                                                                            det
                                                                        );
                 det = det + ((adx * adx + ady * ady) * (bdx * cdytail + cdy * bdxtail - (bdy * cdxtail + cdx * bdytail)) + 2.0f64 * (adx * adxtail + ady * adytail) * (bdx * cdy - bdy * cdx) + ((bdx * bdx + bdy * bdy) * (cdx * adytail + ady * cdxtail - (cdy * adxtail + adx * cdytail)) + 2.0f64 * (bdx * bdxtail + bdy * bdytail) * (cdx * ady - cdy * adx)) + ((cdx * cdx + cdy * cdy) * (adx * bdytail + bdy * adxtail - (ady * bdxtail + bdx * adytail)) + 2.0f64 * (cdx * cdxtail + cdy * cdytail) * (adx * bdy - ady * bdx)));
                 (if det >= errbound || -det >= errbound {
                      det
                  } else {
                      finnow = fin1.as_mut_ptr();
                      finother = fin2.as_mut_ptr();
                      if bdxtail != 0.0f64 || bdytail != 0.0f64 || cdxtail != 0.0f64 || cdytail != 0.0f64 {
                          adxadx1 = adx * adx;
                          c = self.splitter * adx;
                          abig = c - adx;
                          ahi = c - abig;
                          alo = adx - ahi;
                          err1 = adxadx1 - ahi * ahi;
                          err3 = err1 - (ahi + ahi) * alo;
                          adxadx0 = alo * alo - err3;
                          adyady1 = ady * ady;
                          c = self.splitter * ady;
                          abig = c - ady;
                          ahi = c - abig;
                          alo = ady - ahi;
                          err1 = adyady1 - ahi * ahi;
                          err3 = err1 - (ahi + ahi) * alo;
                          adyady0 = alo * alo - err3;
                          _i = adxadx0 + adyady0;
                          bvirt = _i - adxadx0;
                          avirt = _i - bvirt;
                          bround = adyady0 - bvirt;
                          around = adxadx0 - avirt;
                          aa[0usize] = around + bround;
                          _j = adxadx1 + _i;
                          bvirt = _j - adxadx1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = adxadx1 - avirt;
                          _0 = around + bround;
                          _i = _0 + adyady1;
                          bvirt = _i - _0;
                          avirt = _i - bvirt;
                          bround = adyady1 - bvirt;
                          around = _0 - avirt;
                          aa[1usize] = around + bround;
                          aa3 = _j + _i;
                          bvirt = aa3 - _j;
                          avirt = aa3 - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          aa[2usize] = around + bround;
                          aa[3usize] = aa3;
                      }
                      if cdxtail != 0.0f64 || cdytail != 0.0f64 || adxtail != 0.0f64 || adytail != 0.0f64 {
                          bdxbdx1 = bdx * bdx;
                          c = self.splitter * bdx;
                          abig = c - bdx;
                          ahi = c - abig;
                          alo = bdx - ahi;
                          err1 = bdxbdx1 - ahi * ahi;
                          err3 = err1 - (ahi + ahi) * alo;
                          bdxbdx0 = alo * alo - err3;
                          bdybdy1 = bdy * bdy;
                          c = self.splitter * bdy;
                          abig = c - bdy;
                          ahi = c - abig;
                          alo = bdy - ahi;
                          err1 = bdybdy1 - ahi * ahi;
                          err3 = err1 - (ahi + ahi) * alo;
                          bdybdy0 = alo * alo - err3;
                          _i = bdxbdx0 + bdybdy0;
                          bvirt = _i - bdxbdx0;
                          avirt = _i - bvirt;
                          bround = bdybdy0 - bvirt;
                          around = bdxbdx0 - avirt;
                          bb[0usize] = around + bround;
                          _j = bdxbdx1 + _i;
                          bvirt = _j - bdxbdx1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = bdxbdx1 - avirt;
                          _0 = around + bround;
                          _i = _0 + bdybdy1;
                          bvirt = _i - _0;
                          avirt = _i - bvirt;
                          bround = bdybdy1 - bvirt;
                          around = _0 - avirt;
                          bb[1usize] = around + bround;
                          bb3 = _j + _i;
                          bvirt = bb3 - _j;
                          avirt = bb3 - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          bb[2usize] = around + bround;
                          bb[3usize] = bb3;
                      }
                      if adxtail != 0.0f64 || adytail != 0.0f64 || bdxtail != 0.0f64 || bdytail != 0.0f64 {
                          cdxcdx1 = cdx * cdx;
                          c = self.splitter * cdx;
                          abig = c - cdx;
                          ahi = c - abig;
                          alo = cdx - ahi;
                          err1 = cdxcdx1 - ahi * ahi;
                          err3 = err1 - (ahi + ahi) * alo;
                          cdxcdx0 = alo * alo - err3;
                          cdycdy1 = cdy * cdy;
                          c = self.splitter * cdy;
                          abig = c - cdy;
                          ahi = c - abig;
                          alo = cdy - ahi;
                          err1 = cdycdy1 - ahi * ahi;
                          err3 = err1 - (ahi + ahi) * alo;
                          cdycdy0 = alo * alo - err3;
                          _i = cdxcdx0 + cdycdy0;
                          bvirt = _i - cdxcdx0;
                          avirt = _i - bvirt;
                          bround = cdycdy0 - bvirt;
                          around = cdxcdx0 - avirt;
                          cc[0usize] = around + bround;
                          _j = cdxcdx1 + _i;
                          bvirt = _j - cdxcdx1;
                          avirt = _j - bvirt;
                          bround = _i - bvirt;
                          around = cdxcdx1 - avirt;
                          _0 = around + bround;
                          _i = _0 + cdycdy1;
                          bvirt = _i - _0;
                          avirt = _i - bvirt;
                          bround = cdycdy1 - bvirt;
                          around = _0 - avirt;
                          cc[1usize] = around + bround;
                          cc3 = _j + _i;
                          bvirt = cc3 - _j;
                          avirt = cc3 - bvirt;
                          bround = _i - bvirt;
                          around = _j - avirt;
                          cc[2usize] = around + bround;
                          cc[3usize] = cc3;
                      }
                      if adxtail != 0.0f64 {
                          axtbclen = self.scale_expansion_zeroelim(
                                         4i32,
                                         bc.as_ptr(),
                                         adxtail,
                                         axtbc.as_mut_ptr()
                                     );
                          temp16alen = self.scale_expansion_zeroelim(
                                           axtbclen,
                                           axtbc.as_ptr(),
                                           2.0f64 * adx,
                                           temp16a.as_mut_ptr()
                                       );
                          axtcclen = self.scale_expansion_zeroelim(
                                         4i32,
                                         cc.as_ptr(),
                                         adxtail,
                                         axtcc.as_mut_ptr()
                                     );
                          temp16blen = self.scale_expansion_zeroelim(
                                           axtcclen,
                                           axtcc.as_ptr(),
                                           bdy,
                                           temp16b.as_mut_ptr()
                                       );
                          axtbblen = self.scale_expansion_zeroelim(
                                         4i32,
                                         bb.as_ptr(),
                                         adxtail,
                                         axtbb.as_mut_ptr()
                                     );
                          temp16clen = self.scale_expansion_zeroelim(
                                           axtbblen,
                                           axtbb.as_ptr(),
                                           -cdy,
                                           temp16c.as_mut_ptr()
                                       );
                          temp32alen = fast_expansion_sum_zeroelim(
                                           temp16alen,
                                           temp16a.as_ptr(),
                                           temp16blen,
                                           temp16b.as_ptr(),
                                           temp32a.as_mut_ptr()
                                       );
                          temp48len = fast_expansion_sum_zeroelim(
                                          temp16clen,
                                          temp16c.as_ptr(),
                                          temp32alen,
                                          temp32a.as_ptr(),
                                          temp48.as_mut_ptr()
                                      );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          temp48len,
                                          temp48.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if adytail != 0.0f64 {
                          aytbclen = self.scale_expansion_zeroelim(
                                         4i32,
                                         bc.as_ptr(),
                                         adytail,
                                         aytbc.as_mut_ptr()
                                     );
                          temp16alen = self.scale_expansion_zeroelim(
                                           aytbclen,
                                           aytbc.as_ptr(),
                                           2.0f64 * ady,
                                           temp16a.as_mut_ptr()
                                       );
                          aytbblen = self.scale_expansion_zeroelim(
                                         4i32,
                                         bb.as_ptr(),
                                         adytail,
                                         aytbb.as_mut_ptr()
                                     );
                          temp16blen = self.scale_expansion_zeroelim(
                                           aytbblen,
                                           aytbb.as_ptr(),
                                           cdx,
                                           temp16b.as_mut_ptr()
                                       );
                          aytcclen = self.scale_expansion_zeroelim(
                                         4i32,
                                         cc.as_ptr(),
                                         adytail,
                                         aytcc.as_mut_ptr()
                                     );
                          temp16clen = self.scale_expansion_zeroelim(
                                           aytcclen,
                                           aytcc.as_ptr(),
                                           -bdx,
                                           temp16c.as_mut_ptr()
                                       );
                          temp32alen = fast_expansion_sum_zeroelim(
                                           temp16alen,
                                           temp16a.as_ptr(),
                                           temp16blen,
                                           temp16b.as_ptr(),
                                           temp32a.as_mut_ptr()
                                       );
                          temp48len = fast_expansion_sum_zeroelim(
                                          temp16clen,
                                          temp16c.as_ptr(),
                                          temp32alen,
                                          temp32a.as_ptr(),
                                          temp48.as_mut_ptr()
                                      );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          temp48len,
                                          temp48.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if bdxtail != 0.0f64 {
                          bxtcalen = self.scale_expansion_zeroelim(
                                         4i32,
                                         ca.as_ptr(),
                                         bdxtail,
                                         bxtca.as_mut_ptr()
                                     );
                          temp16alen = self.scale_expansion_zeroelim(
                                           bxtcalen,
                                           bxtca.as_ptr(),
                                           2.0f64 * bdx,
                                           temp16a.as_mut_ptr()
                                       );
                          bxtaalen = self.scale_expansion_zeroelim(
                                         4i32,
                                         aa.as_ptr(),
                                         bdxtail,
                                         bxtaa.as_mut_ptr()
                                     );
                          temp16blen = self.scale_expansion_zeroelim(
                                           bxtaalen,
                                           bxtaa.as_ptr(),
                                           cdy,
                                           temp16b.as_mut_ptr()
                                       );
                          bxtcclen = self.scale_expansion_zeroelim(
                                         4i32,
                                         cc.as_ptr(),
                                         bdxtail,
                                         bxtcc.as_mut_ptr()
                                     );
                          temp16clen = self.scale_expansion_zeroelim(
                                           bxtcclen,
                                           bxtcc.as_ptr(),
                                           -ady,
                                           temp16c.as_mut_ptr()
                                       );
                          temp32alen = fast_expansion_sum_zeroelim(
                                           temp16alen,
                                           temp16a.as_ptr(),
                                           temp16blen,
                                           temp16b.as_ptr(),
                                           temp32a.as_mut_ptr()
                                       );
                          temp48len = fast_expansion_sum_zeroelim(
                                          temp16clen,
                                          temp16c.as_ptr(),
                                          temp32alen,
                                          temp32a.as_ptr(),
                                          temp48.as_mut_ptr()
                                      );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          temp48len,
                                          temp48.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if bdytail != 0.0f64 {
                          bytcalen = self.scale_expansion_zeroelim(
                                         4i32,
                                         ca.as_ptr(),
                                         bdytail,
                                         bytca.as_mut_ptr()
                                     );
                          temp16alen = self.scale_expansion_zeroelim(
                                           bytcalen,
                                           bytca.as_ptr(),
                                           2.0f64 * bdy,
                                           temp16a.as_mut_ptr()
                                       );
                          bytcclen = self.scale_expansion_zeroelim(
                                         4i32,
                                         cc.as_ptr(),
                                         bdytail,
                                         bytcc.as_mut_ptr()
                                     );
                          temp16blen = self.scale_expansion_zeroelim(
                                           bytcclen,
                                           bytcc.as_ptr(),
                                           adx,
                                           temp16b.as_mut_ptr()
                                       );
                          bytaalen = self.scale_expansion_zeroelim(
                                         4i32,
                                         aa.as_ptr(),
                                         bdytail,
                                         bytaa.as_mut_ptr()
                                     );
                          temp16clen = self.scale_expansion_zeroelim(
                                           bytaalen,
                                           bytaa.as_ptr(),
                                           -cdx,
                                           temp16c.as_mut_ptr()
                                       );
                          temp32alen = fast_expansion_sum_zeroelim(
                                           temp16alen,
                                           temp16a.as_ptr(),
                                           temp16blen,
                                           temp16b.as_ptr(),
                                           temp32a.as_mut_ptr()
                                       );
                          temp48len = fast_expansion_sum_zeroelim(
                                          temp16clen,
                                          temp16c.as_ptr(),
                                          temp32alen,
                                          temp32a.as_ptr(),
                                          temp48.as_mut_ptr()
                                      );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          temp48len,
                                          temp48.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if cdxtail != 0.0f64 {
                          cxtablen = self.scale_expansion_zeroelim(
                                         4i32,
                                         ab.as_ptr(),
                                         cdxtail,
                                         cxtab.as_mut_ptr()
                                     );
                          temp16alen = self.scale_expansion_zeroelim(
                                           cxtablen,
                                           cxtab.as_ptr(),
                                           2.0f64 * cdx,
                                           temp16a.as_mut_ptr()
                                       );
                          cxtbblen = self.scale_expansion_zeroelim(
                                         4i32,
                                         bb.as_ptr(),
                                         cdxtail,
                                         cxtbb.as_mut_ptr()
                                     );
                          temp16blen = self.scale_expansion_zeroelim(
                                           cxtbblen,
                                           cxtbb.as_ptr(),
                                           ady,
                                           temp16b.as_mut_ptr()
                                       );
                          cxtaalen = self.scale_expansion_zeroelim(
                                         4i32,
                                         aa.as_ptr(),
                                         cdxtail,
                                         cxtaa.as_mut_ptr()
                                     );
                          temp16clen = self.scale_expansion_zeroelim(
                                           cxtaalen,
                                           cxtaa.as_ptr(),
                                           -bdy,
                                           temp16c.as_mut_ptr()
                                       );
                          temp32alen = fast_expansion_sum_zeroelim(
                                           temp16alen,
                                           temp16a.as_ptr(),
                                           temp16blen,
                                           temp16b.as_ptr(),
                                           temp32a.as_mut_ptr()
                                       );
                          temp48len = fast_expansion_sum_zeroelim(
                                          temp16clen,
                                          temp16c.as_ptr(),
                                          temp32alen,
                                          temp32a.as_ptr(),
                                          temp48.as_mut_ptr()
                                      );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          temp48len,
                                          temp48.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if cdytail != 0.0f64 {
                          cytablen = self.scale_expansion_zeroelim(
                                         4i32,
                                         ab.as_ptr(),
                                         cdytail,
                                         cytab.as_mut_ptr()
                                     );
                          temp16alen = self.scale_expansion_zeroelim(
                                           cytablen,
                                           cytab.as_ptr(),
                                           2.0f64 * cdy,
                                           temp16a.as_mut_ptr()
                                       );
                          cytaalen = self.scale_expansion_zeroelim(
                                         4i32,
                                         aa.as_ptr(),
                                         cdytail,
                                         cytaa.as_mut_ptr()
                                     );
                          temp16blen = self.scale_expansion_zeroelim(
                                           cytaalen,
                                           cytaa.as_ptr(),
                                           bdx,
                                           temp16b.as_mut_ptr()
                                       );
                          cytbblen = self.scale_expansion_zeroelim(
                                         4i32,
                                         bb.as_ptr(),
                                         cdytail,
                                         cytbb.as_mut_ptr()
                                     );
                          temp16clen = self.scale_expansion_zeroelim(
                                           cytbblen,
                                           cytbb.as_ptr(),
                                           -adx,
                                           temp16c.as_mut_ptr()
                                       );
                          temp32alen = fast_expansion_sum_zeroelim(
                                           temp16alen,
                                           temp16a.as_ptr(),
                                           temp16blen,
                                           temp16b.as_ptr(),
                                           temp32a.as_mut_ptr()
                                       );
                          temp48len = fast_expansion_sum_zeroelim(
                                          temp16clen,
                                          temp16c.as_ptr(),
                                          temp32alen,
                                          temp32a.as_ptr(),
                                          temp48.as_mut_ptr()
                                      );
                          finlength = fast_expansion_sum_zeroelim(
                                          finlength,
                                          finnow as (*const f64),
                                          temp48len,
                                          temp48.as_ptr(),
                                          finother
                                      );
                          finswap = finnow;
                          finnow = finother;
                          finother = finswap;
                      }
                      if adxtail != 0.0f64 || adytail != 0.0f64 {
                          if bdxtail != 0.0f64 || bdytail != 0.0f64 || cdxtail != 0.0f64 || cdytail != 0.0f64 {
                              ti1 = bdxtail * cdy;
                              c = self.splitter * bdxtail;
                              abig = c - bdxtail;
                              ahi = c - abig;
                              alo = bdxtail - ahi;
                              c = self.splitter * cdy;
                              abig = c - cdy;
                              bhi = c - abig;
                              blo = cdy - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              tj1 = bdx * cdytail;
                              c = self.splitter * bdx;
                              abig = c - bdx;
                              ahi = c - abig;
                              alo = bdx - ahi;
                              c = self.splitter * cdytail;
                              abig = c - cdytail;
                              bhi = c - abig;
                              blo = cdytail - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 + tj0;
                              bvirt = _i - ti0;
                              avirt = _i - bvirt;
                              bround = tj0 - bvirt;
                              around = ti0 - avirt;
                              u[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 + tj1;
                              bvirt = _i - _0;
                              avirt = _i - bvirt;
                              bround = tj1 - bvirt;
                              around = _0 - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _i;
                              bvirt = u3 - _j;
                              avirt = u3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              u[2usize] = around + bround;
                              u[3usize] = u3;
                              negate = -bdy;
                              ti1 = cdxtail * negate;
                              c = self.splitter * cdxtail;
                              abig = c - cdxtail;
                              ahi = c - abig;
                              alo = cdxtail - ahi;
                              c = self.splitter * negate;
                              abig = c - negate;
                              bhi = c - abig;
                              blo = negate - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              negate = -bdytail;
                              tj1 = cdx * negate;
                              c = self.splitter * cdx;
                              abig = c - cdx;
                              ahi = c - abig;
                              alo = cdx - ahi;
                              c = self.splitter * negate;
                              abig = c - negate;
                              bhi = c - abig;
                              blo = negate - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 + tj0;
                              bvirt = _i - ti0;
                              avirt = _i - bvirt;
                              bround = tj0 - bvirt;
                              around = ti0 - avirt;
                              v[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 + tj1;
                              bvirt = _i - _0;
                              avirt = _i - bvirt;
                              bround = tj1 - bvirt;
                              around = _0 - avirt;
                              v[1usize] = around + bround;
                              v3 = _j + _i;
                              bvirt = v3 - _j;
                              avirt = v3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              v[2usize] = around + bround;
                              v[3usize] = v3;
                              bctlen = fast_expansion_sum_zeroelim(
                                           4i32,
                                           u.as_ptr(),
                                           4i32,
                                           v.as_ptr(),
                                           bct.as_mut_ptr()
                                       );
                              ti1 = bdxtail * cdytail;
                              c = self.splitter * bdxtail;
                              abig = c - bdxtail;
                              ahi = c - abig;
                              alo = bdxtail - ahi;
                              c = self.splitter * cdytail;
                              abig = c - cdytail;
                              bhi = c - abig;
                              blo = cdytail - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              tj1 = cdxtail * bdytail;
                              c = self.splitter * cdxtail;
                              abig = c - cdxtail;
                              ahi = c - abig;
                              alo = cdxtail - ahi;
                              c = self.splitter * bdytail;
                              abig = c - bdytail;
                              bhi = c - abig;
                              blo = bdytail - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 - tj0;
                              bvirt = ti0 - _i;
                              avirt = _i + bvirt;
                              bround = bvirt - tj0;
                              around = ti0 - avirt;
                              bctt[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 - tj1;
                              bvirt = _0 - _i;
                              avirt = _i + bvirt;
                              bround = bvirt - tj1;
                              around = _0 - avirt;
                              bctt[1usize] = around + bround;
                              bctt3 = _j + _i;
                              bvirt = bctt3 - _j;
                              avirt = bctt3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              bctt[2usize] = around + bround;
                              bctt[3usize] = bctt3;
                              bcttlen = 4i32;
                          } else {
                              bct[0usize] = 0.0f64;
                              bctlen = 1i32;
                              bctt[0usize] = 0.0f64;
                              bcttlen = 1i32;
                          }
                          if adxtail != 0.0f64 {
                              temp16alen = self.scale_expansion_zeroelim(
                                               axtbclen,
                                               axtbc.as_ptr(),
                                               adxtail,
                                               temp16a.as_mut_ptr()
                                           );
                              axtbctlen = self.scale_expansion_zeroelim(
                                              bctlen,
                                              bct.as_ptr(),
                                              adxtail,
                                              axtbct.as_mut_ptr()
                                          );
                              temp32alen = self.scale_expansion_zeroelim(
                                               axtbctlen,
                                               axtbct.as_ptr(),
                                               2.0f64 * adx,
                                               temp32a.as_mut_ptr()
                                           );
                              temp48len = fast_expansion_sum_zeroelim(
                                              temp16alen,
                                              temp16a.as_ptr(),
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp48.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp48len,
                                              temp48.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if bdytail != 0.0f64 {
                                  temp8len = self.scale_expansion_zeroelim(
                                                 4i32,
                                                 cc.as_ptr(),
                                                 adxtail,
                                                 temp8.as_mut_ptr()
                                             );
                                  temp16alen = self.scale_expansion_zeroelim(
                                                   temp8len,
                                                   temp8.as_ptr(),
                                                   bdytail,
                                                   temp16a.as_mut_ptr()
                                               );
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  temp16alen,
                                                  temp16a.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                              if cdytail != 0.0f64 {
                                  temp8len = self.scale_expansion_zeroelim(
                                                 4i32,
                                                 bb.as_ptr(),
                                                 -adxtail,
                                                 temp8.as_mut_ptr()
                                             );
                                  temp16alen = self.scale_expansion_zeroelim(
                                                   temp8len,
                                                   temp8.as_ptr(),
                                                   cdytail,
                                                   temp16a.as_mut_ptr()
                                               );
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  temp16alen,
                                                  temp16a.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                              temp32alen = self.scale_expansion_zeroelim(
                                               axtbctlen,
                                               axtbct.as_ptr(),
                                               adxtail,
                                               temp32a.as_mut_ptr()
                                           );
                              axtbcttlen = self.scale_expansion_zeroelim(
                                               bcttlen,
                                               bctt.as_ptr(),
                                               adxtail,
                                               axtbctt.as_mut_ptr()
                                           );
                              temp16alen = self.scale_expansion_zeroelim(
                                               axtbcttlen,
                                               axtbctt.as_ptr(),
                                               2.0f64 * adx,
                                               temp16a.as_mut_ptr()
                                           );
                              temp16blen = self.scale_expansion_zeroelim(
                                               axtbcttlen,
                                               axtbctt.as_ptr(),
                                               adxtail,
                                               temp16b.as_mut_ptr()
                                           );
                              temp32blen = fast_expansion_sum_zeroelim(
                                               temp16alen,
                                               temp16a.as_ptr(),
                                               temp16blen,
                                               temp16b.as_ptr(),
                                               temp32b.as_mut_ptr()
                                           );
                              temp64len = fast_expansion_sum_zeroelim(
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp32blen,
                                              temp32b.as_ptr(),
                                              temp64.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp64len,
                                              temp64.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                          }
                          if adytail != 0.0f64 {
                              temp16alen = self.scale_expansion_zeroelim(
                                               aytbclen,
                                               aytbc.as_ptr(),
                                               adytail,
                                               temp16a.as_mut_ptr()
                                           );
                              aytbctlen = self.scale_expansion_zeroelim(
                                              bctlen,
                                              bct.as_ptr(),
                                              adytail,
                                              aytbct.as_mut_ptr()
                                          );
                              temp32alen = self.scale_expansion_zeroelim(
                                               aytbctlen,
                                               aytbct.as_ptr(),
                                               2.0f64 * ady,
                                               temp32a.as_mut_ptr()
                                           );
                              temp48len = fast_expansion_sum_zeroelim(
                                              temp16alen,
                                              temp16a.as_ptr(),
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp48.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp48len,
                                              temp48.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              temp32alen = self.scale_expansion_zeroelim(
                                               aytbctlen,
                                               aytbct.as_ptr(),
                                               adytail,
                                               temp32a.as_mut_ptr()
                                           );
                              aytbcttlen = self.scale_expansion_zeroelim(
                                               bcttlen,
                                               bctt.as_ptr(),
                                               adytail,
                                               aytbctt.as_mut_ptr()
                                           );
                              temp16alen = self.scale_expansion_zeroelim(
                                               aytbcttlen,
                                               aytbctt.as_ptr(),
                                               2.0f64 * ady,
                                               temp16a.as_mut_ptr()
                                           );
                              temp16blen = self.scale_expansion_zeroelim(
                                               aytbcttlen,
                                               aytbctt.as_ptr(),
                                               adytail,
                                               temp16b.as_mut_ptr()
                                           );
                              temp32blen = fast_expansion_sum_zeroelim(
                                               temp16alen,
                                               temp16a.as_ptr(),
                                               temp16blen,
                                               temp16b.as_ptr(),
                                               temp32b.as_mut_ptr()
                                           );
                              temp64len = fast_expansion_sum_zeroelim(
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp32blen,
                                              temp32b.as_ptr(),
                                              temp64.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp64len,
                                              temp64.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                          }
                      }
                      if bdxtail != 0.0f64 || bdytail != 0.0f64 {
                          if cdxtail != 0.0f64 || cdytail != 0.0f64 || adxtail != 0.0f64 || adytail != 0.0f64 {
                              ti1 = cdxtail * ady;
                              c = self.splitter * cdxtail;
                              abig = c - cdxtail;
                              ahi = c - abig;
                              alo = cdxtail - ahi;
                              c = self.splitter * ady;
                              abig = c - ady;
                              bhi = c - abig;
                              blo = ady - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              tj1 = cdx * adytail;
                              c = self.splitter * cdx;
                              abig = c - cdx;
                              ahi = c - abig;
                              alo = cdx - ahi;
                              c = self.splitter * adytail;
                              abig = c - adytail;
                              bhi = c - abig;
                              blo = adytail - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 + tj0;
                              bvirt = _i - ti0;
                              avirt = _i - bvirt;
                              bround = tj0 - bvirt;
                              around = ti0 - avirt;
                              u[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 + tj1;
                              bvirt = _i - _0;
                              avirt = _i - bvirt;
                              bround = tj1 - bvirt;
                              around = _0 - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _i;
                              bvirt = u3 - _j;
                              avirt = u3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              u[2usize] = around + bround;
                              u[3usize] = u3;
                              negate = -cdy;
                              ti1 = adxtail * negate;
                              c = self.splitter * adxtail;
                              abig = c - adxtail;
                              ahi = c - abig;
                              alo = adxtail - ahi;
                              c = self.splitter * negate;
                              abig = c - negate;
                              bhi = c - abig;
                              blo = negate - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              negate = -cdytail;
                              tj1 = adx * negate;
                              c = self.splitter * adx;
                              abig = c - adx;
                              ahi = c - abig;
                              alo = adx - ahi;
                              c = self.splitter * negate;
                              abig = c - negate;
                              bhi = c - abig;
                              blo = negate - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 + tj0;
                              bvirt = _i - ti0;
                              avirt = _i - bvirt;
                              bround = tj0 - bvirt;
                              around = ti0 - avirt;
                              v[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 + tj1;
                              bvirt = _i - _0;
                              avirt = _i - bvirt;
                              bround = tj1 - bvirt;
                              around = _0 - avirt;
                              v[1usize] = around + bround;
                              v3 = _j + _i;
                              bvirt = v3 - _j;
                              avirt = v3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              v[2usize] = around + bround;
                              v[3usize] = v3;
                              catlen = fast_expansion_sum_zeroelim(
                                           4i32,
                                           u.as_ptr(),
                                           4i32,
                                           v.as_ptr(),
                                           cat.as_mut_ptr()
                                       );
                              ti1 = cdxtail * adytail;
                              c = self.splitter * cdxtail;
                              abig = c - cdxtail;
                              ahi = c - abig;
                              alo = cdxtail - ahi;
                              c = self.splitter * adytail;
                              abig = c - adytail;
                              bhi = c - abig;
                              blo = adytail - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              tj1 = adxtail * cdytail;
                              c = self.splitter * adxtail;
                              abig = c - adxtail;
                              ahi = c - abig;
                              alo = adxtail - ahi;
                              c = self.splitter * cdytail;
                              abig = c - cdytail;
                              bhi = c - abig;
                              blo = cdytail - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 - tj0;
                              bvirt = ti0 - _i;
                              avirt = _i + bvirt;
                              bround = bvirt - tj0;
                              around = ti0 - avirt;
                              catt[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 - tj1;
                              bvirt = _0 - _i;
                              avirt = _i + bvirt;
                              bround = bvirt - tj1;
                              around = _0 - avirt;
                              catt[1usize] = around + bround;
                              catt3 = _j + _i;
                              bvirt = catt3 - _j;
                              avirt = catt3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              catt[2usize] = around + bround;
                              catt[3usize] = catt3;
                              cattlen = 4i32;
                          } else {
                              cat[0usize] = 0.0f64;
                              catlen = 1i32;
                              catt[0usize] = 0.0f64;
                              cattlen = 1i32;
                          }
                          if bdxtail != 0.0f64 {
                              temp16alen = self.scale_expansion_zeroelim(
                                               bxtcalen,
                                               bxtca.as_ptr(),
                                               bdxtail,
                                               temp16a.as_mut_ptr()
                                           );
                              bxtcatlen = self.scale_expansion_zeroelim(
                                              catlen,
                                              cat.as_ptr(),
                                              bdxtail,
                                              bxtcat.as_mut_ptr()
                                          );
                              temp32alen = self.scale_expansion_zeroelim(
                                               bxtcatlen,
                                               bxtcat.as_ptr(),
                                               2.0f64 * bdx,
                                               temp32a.as_mut_ptr()
                                           );
                              temp48len = fast_expansion_sum_zeroelim(
                                              temp16alen,
                                              temp16a.as_ptr(),
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp48.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp48len,
                                              temp48.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if cdytail != 0.0f64 {
                                  temp8len = self.scale_expansion_zeroelim(
                                                 4i32,
                                                 aa.as_ptr(),
                                                 bdxtail,
                                                 temp8.as_mut_ptr()
                                             );
                                  temp16alen = self.scale_expansion_zeroelim(
                                                   temp8len,
                                                   temp8.as_ptr(),
                                                   cdytail,
                                                   temp16a.as_mut_ptr()
                                               );
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  temp16alen,
                                                  temp16a.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                              if adytail != 0.0f64 {
                                  temp8len = self.scale_expansion_zeroelim(
                                                 4i32,
                                                 cc.as_ptr(),
                                                 -bdxtail,
                                                 temp8.as_mut_ptr()
                                             );
                                  temp16alen = self.scale_expansion_zeroelim(
                                                   temp8len,
                                                   temp8.as_ptr(),
                                                   adytail,
                                                   temp16a.as_mut_ptr()
                                               );
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  temp16alen,
                                                  temp16a.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                              temp32alen = self.scale_expansion_zeroelim(
                                               bxtcatlen,
                                               bxtcat.as_ptr(),
                                               bdxtail,
                                               temp32a.as_mut_ptr()
                                           );
                              bxtcattlen = self.scale_expansion_zeroelim(
                                               cattlen,
                                               catt.as_ptr(),
                                               bdxtail,
                                               bxtcatt.as_mut_ptr()
                                           );
                              temp16alen = self.scale_expansion_zeroelim(
                                               bxtcattlen,
                                               bxtcatt.as_ptr(),
                                               2.0f64 * bdx,
                                               temp16a.as_mut_ptr()
                                           );
                              temp16blen = self.scale_expansion_zeroelim(
                                               bxtcattlen,
                                               bxtcatt.as_ptr(),
                                               bdxtail,
                                               temp16b.as_mut_ptr()
                                           );
                              temp32blen = fast_expansion_sum_zeroelim(
                                               temp16alen,
                                               temp16a.as_ptr(),
                                               temp16blen,
                                               temp16b.as_ptr(),
                                               temp32b.as_mut_ptr()
                                           );
                              temp64len = fast_expansion_sum_zeroelim(
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp32blen,
                                              temp32b.as_ptr(),
                                              temp64.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp64len,
                                              temp64.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                          }
                          if bdytail != 0.0f64 {
                              temp16alen = self.scale_expansion_zeroelim(
                                               bytcalen,
                                               bytca.as_ptr(),
                                               bdytail,
                                               temp16a.as_mut_ptr()
                                           );
                              bytcatlen = self.scale_expansion_zeroelim(
                                              catlen,
                                              cat.as_ptr(),
                                              bdytail,
                                              bytcat.as_mut_ptr()
                                          );
                              temp32alen = self.scale_expansion_zeroelim(
                                               bytcatlen,
                                               bytcat.as_ptr(),
                                               2.0f64 * bdy,
                                               temp32a.as_mut_ptr()
                                           );
                              temp48len = fast_expansion_sum_zeroelim(
                                              temp16alen,
                                              temp16a.as_ptr(),
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp48.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp48len,
                                              temp48.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              temp32alen = self.scale_expansion_zeroelim(
                                               bytcatlen,
                                               bytcat.as_ptr(),
                                               bdytail,
                                               temp32a.as_mut_ptr()
                                           );
                              bytcattlen = self.scale_expansion_zeroelim(
                                               cattlen,
                                               catt.as_ptr(),
                                               bdytail,
                                               bytcatt.as_mut_ptr()
                                           );
                              temp16alen = self.scale_expansion_zeroelim(
                                               bytcattlen,
                                               bytcatt.as_ptr(),
                                               2.0f64 * bdy,
                                               temp16a.as_mut_ptr()
                                           );
                              temp16blen = self.scale_expansion_zeroelim(
                                               bytcattlen,
                                               bytcatt.as_ptr(),
                                               bdytail,
                                               temp16b.as_mut_ptr()
                                           );
                              temp32blen = fast_expansion_sum_zeroelim(
                                               temp16alen,
                                               temp16a.as_ptr(),
                                               temp16blen,
                                               temp16b.as_ptr(),
                                               temp32b.as_mut_ptr()
                                           );
                              temp64len = fast_expansion_sum_zeroelim(
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp32blen,
                                              temp32b.as_ptr(),
                                              temp64.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp64len,
                                              temp64.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                          }
                      }
                      if cdxtail != 0.0f64 || cdytail != 0.0f64 {
                          if adxtail != 0.0f64 || adytail != 0.0f64 || bdxtail != 0.0f64 || bdytail != 0.0f64 {
                              ti1 = adxtail * bdy;
                              c = self.splitter * adxtail;
                              abig = c - adxtail;
                              ahi = c - abig;
                              alo = adxtail - ahi;
                              c = self.splitter * bdy;
                              abig = c - bdy;
                              bhi = c - abig;
                              blo = bdy - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              tj1 = adx * bdytail;
                              c = self.splitter * adx;
                              abig = c - adx;
                              ahi = c - abig;
                              alo = adx - ahi;
                              c = self.splitter * bdytail;
                              abig = c - bdytail;
                              bhi = c - abig;
                              blo = bdytail - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 + tj0;
                              bvirt = _i - ti0;
                              avirt = _i - bvirt;
                              bround = tj0 - bvirt;
                              around = ti0 - avirt;
                              u[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 + tj1;
                              bvirt = _i - _0;
                              avirt = _i - bvirt;
                              bround = tj1 - bvirt;
                              around = _0 - avirt;
                              u[1usize] = around + bround;
                              u3 = _j + _i;
                              bvirt = u3 - _j;
                              avirt = u3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              u[2usize] = around + bround;
                              u[3usize] = u3;
                              negate = -ady;
                              ti1 = bdxtail * negate;
                              c = self.splitter * bdxtail;
                              abig = c - bdxtail;
                              ahi = c - abig;
                              alo = bdxtail - ahi;
                              c = self.splitter * negate;
                              abig = c - negate;
                              bhi = c - abig;
                              blo = negate - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              negate = -adytail;
                              tj1 = bdx * negate;
                              c = self.splitter * bdx;
                              abig = c - bdx;
                              ahi = c - abig;
                              alo = bdx - ahi;
                              c = self.splitter * negate;
                              abig = c - negate;
                              bhi = c - abig;
                              blo = negate - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 + tj0;
                              bvirt = _i - ti0;
                              avirt = _i - bvirt;
                              bround = tj0 - bvirt;
                              around = ti0 - avirt;
                              v[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 + tj1;
                              bvirt = _i - _0;
                              avirt = _i - bvirt;
                              bround = tj1 - bvirt;
                              around = _0 - avirt;
                              v[1usize] = around + bround;
                              v3 = _j + _i;
                              bvirt = v3 - _j;
                              avirt = v3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              v[2usize] = around + bround;
                              v[3usize] = v3;
                              abtlen = fast_expansion_sum_zeroelim(
                                           4i32,
                                           u.as_ptr(),
                                           4i32,
                                           v.as_ptr(),
                                           abt.as_mut_ptr()
                                       );
                              ti1 = adxtail * bdytail;
                              c = self.splitter * adxtail;
                              abig = c - adxtail;
                              ahi = c - abig;
                              alo = adxtail - ahi;
                              c = self.splitter * bdytail;
                              abig = c - bdytail;
                              bhi = c - abig;
                              blo = bdytail - bhi;
                              err1 = ti1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              ti0 = alo * blo - err3;
                              tj1 = bdxtail * adytail;
                              c = self.splitter * bdxtail;
                              abig = c - bdxtail;
                              ahi = c - abig;
                              alo = bdxtail - ahi;
                              c = self.splitter * adytail;
                              abig = c - adytail;
                              bhi = c - abig;
                              blo = adytail - bhi;
                              err1 = tj1 - ahi * bhi;
                              err2 = err1 - alo * bhi;
                              err3 = err2 - ahi * blo;
                              tj0 = alo * blo - err3;
                              _i = ti0 - tj0;
                              bvirt = ti0 - _i;
                              avirt = _i + bvirt;
                              bround = bvirt - tj0;
                              around = ti0 - avirt;
                              abtt[0usize] = around + bround;
                              _j = ti1 + _i;
                              bvirt = _j - ti1;
                              avirt = _j - bvirt;
                              bround = _i - bvirt;
                              around = ti1 - avirt;
                              _0 = around + bround;
                              _i = _0 - tj1;
                              bvirt = _0 - _i;
                              avirt = _i + bvirt;
                              bround = bvirt - tj1;
                              around = _0 - avirt;
                              abtt[1usize] = around + bround;
                              abtt3 = _j + _i;
                              bvirt = abtt3 - _j;
                              avirt = abtt3 - bvirt;
                              bround = _i - bvirt;
                              around = _j - avirt;
                              abtt[2usize] = around + bround;
                              abtt[3usize] = abtt3;
                              abttlen = 4i32;
                          } else {
                              abt[0usize] = 0.0f64;
                              abtlen = 1i32;
                              abtt[0usize] = 0.0f64;
                              abttlen = 1i32;
                          }
                          if cdxtail != 0.0f64 {
                              temp16alen = self.scale_expansion_zeroelim(
                                               cxtablen,
                                               cxtab.as_ptr(),
                                               cdxtail,
                                               temp16a.as_mut_ptr()
                                           );
                              cxtabtlen = self.scale_expansion_zeroelim(
                                              abtlen,
                                              abt.as_ptr(),
                                              cdxtail,
                                              cxtabt.as_mut_ptr()
                                          );
                              temp32alen = self.scale_expansion_zeroelim(
                                               cxtabtlen,
                                               cxtabt.as_ptr(),
                                               2.0f64 * cdx,
                                               temp32a.as_mut_ptr()
                                           );
                              temp48len = fast_expansion_sum_zeroelim(
                                              temp16alen,
                                              temp16a.as_ptr(),
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp48.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp48len,
                                              temp48.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              if adytail != 0.0f64 {
                                  temp8len = self.scale_expansion_zeroelim(
                                                 4i32,
                                                 bb.as_ptr(),
                                                 cdxtail,
                                                 temp8.as_mut_ptr()
                                             );
                                  temp16alen = self.scale_expansion_zeroelim(
                                                   temp8len,
                                                   temp8.as_ptr(),
                                                   adytail,
                                                   temp16a.as_mut_ptr()
                                               );
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  temp16alen,
                                                  temp16a.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                              if bdytail != 0.0f64 {
                                  temp8len = self.scale_expansion_zeroelim(
                                                 4i32,
                                                 aa.as_ptr(),
                                                 -cdxtail,
                                                 temp8.as_mut_ptr()
                                             );
                                  temp16alen = self.scale_expansion_zeroelim(
                                                   temp8len,
                                                   temp8.as_ptr(),
                                                   bdytail,
                                                   temp16a.as_mut_ptr()
                                               );
                                  finlength = fast_expansion_sum_zeroelim(
                                                  finlength,
                                                  finnow as (*const f64),
                                                  temp16alen,
                                                  temp16a.as_ptr(),
                                                  finother
                                              );
                                  finswap = finnow;
                                  finnow = finother;
                                  finother = finswap;
                              }
                              temp32alen = self.scale_expansion_zeroelim(
                                               cxtabtlen,
                                               cxtabt.as_ptr(),
                                               cdxtail,
                                               temp32a.as_mut_ptr()
                                           );
                              cxtabttlen = self.scale_expansion_zeroelim(
                                               abttlen,
                                               abtt.as_ptr(),
                                               cdxtail,
                                               cxtabtt.as_mut_ptr()
                                           );
                              temp16alen = self.scale_expansion_zeroelim(
                                               cxtabttlen,
                                               cxtabtt.as_ptr(),
                                               2.0f64 * cdx,
                                               temp16a.as_mut_ptr()
                                           );
                              temp16blen = self.scale_expansion_zeroelim(
                                               cxtabttlen,
                                               cxtabtt.as_ptr(),
                                               cdxtail,
                                               temp16b.as_mut_ptr()
                                           );
                              temp32blen = fast_expansion_sum_zeroelim(
                                               temp16alen,
                                               temp16a.as_ptr(),
                                               temp16blen,
                                               temp16b.as_ptr(),
                                               temp32b.as_mut_ptr()
                                           );
                              temp64len = fast_expansion_sum_zeroelim(
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp32blen,
                                              temp32b.as_ptr(),
                                              temp64.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp64len,
                                              temp64.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                          }
                          if cdytail != 0.0f64 {
                              temp16alen = self.scale_expansion_zeroelim(
                                               cytablen,
                                               cytab.as_ptr(),
                                               cdytail,
                                               temp16a.as_mut_ptr()
                                           );
                              cytabtlen = self.scale_expansion_zeroelim(
                                              abtlen,
                                              abt.as_ptr(),
                                              cdytail,
                                              cytabt.as_mut_ptr()
                                          );
                              temp32alen = self.scale_expansion_zeroelim(
                                               cytabtlen,
                                               cytabt.as_ptr(),
                                               2.0f64 * cdy,
                                               temp32a.as_mut_ptr()
                                           );
                              temp48len = fast_expansion_sum_zeroelim(
                                              temp16alen,
                                              temp16a.as_ptr(),
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp48.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp48len,
                                              temp48.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                              temp32alen = self.scale_expansion_zeroelim(
                                               cytabtlen,
                                               cytabt.as_ptr(),
                                               cdytail,
                                               temp32a.as_mut_ptr()
                                           );
                              cytabttlen = self.scale_expansion_zeroelim(
                                               abttlen,
                                               abtt.as_ptr(),
                                               cdytail,
                                               cytabtt.as_mut_ptr()
                                           );
                              temp16alen = self.scale_expansion_zeroelim(
                                               cytabttlen,
                                               cytabtt.as_ptr(),
                                               2.0f64 * cdy,
                                               temp16a.as_mut_ptr()
                                           );
                              temp16blen = self.scale_expansion_zeroelim(
                                               cytabttlen,
                                               cytabtt.as_ptr(),
                                               cdytail,
                                               temp16b.as_mut_ptr()
                                           );
                              temp32blen = fast_expansion_sum_zeroelim(
                                               temp16alen,
                                               temp16a.as_ptr(),
                                               temp16blen,
                                               temp16b.as_ptr(),
                                               temp32b.as_mut_ptr()
                                           );
                              temp64len = fast_expansion_sum_zeroelim(
                                              temp32alen,
                                              temp32a.as_ptr(),
                                              temp32blen,
                                              temp32b.as_ptr(),
                                              temp64.as_mut_ptr()
                                          );
                              finlength = fast_expansion_sum_zeroelim(
                                              finlength,
                                              finnow as (*const f64),
                                              temp64len,
                                              temp64.as_ptr(),
                                              finother
                                          );
                              finswap = finnow;
                              finnow = finother;
                              finother = finswap;
                          }
                      }
                      *finnow.offset((finlength - 1i32) as (isize))
                  })
             })
        }
    }

    
    pub unsafe fn incircle(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64
    ) -> f64 {
        let mut adx : f64;
        let mut bdx : f64;
        let mut cdx : f64;
        let mut ady : f64;
        let mut bdy : f64;
        let mut cdy : f64;
        let mut bdxcdy : f64;
        let mut cdxbdy : f64;
        let mut cdxady : f64;
        let mut adxcdy : f64;
        let mut adxbdy : f64;
        let mut bdxady : f64;
        let mut alift : f64;
        let mut blift : f64;
        let mut clift : f64;
        let mut det : f64;
        let mut permanent : f64;
        let mut errbound : f64;
        adx = *pa.offset(0isize) - *pd.offset(0isize);
        bdx = *pb.offset(0isize) - *pd.offset(0isize);
        cdx = *pc.offset(0isize) - *pd.offset(0isize);
        ady = *pa.offset(1isize) - *pd.offset(1isize);
        bdy = *pb.offset(1isize) - *pd.offset(1isize);
        cdy = *pc.offset(1isize) - *pd.offset(1isize);
        bdxcdy = bdx * cdy;
        cdxbdy = cdx * bdy;
        alift = adx * adx + ady * ady;
        cdxady = cdx * ady;
        adxcdy = adx * cdy;
        blift = bdx * bdx + bdy * bdy;
        adxbdy = adx * bdy;
        bdxady = bdx * ady;
        clift = cdx * cdx + cdy * cdy;
        det = alift * (bdxcdy - cdxbdy) + blift * (cdxady - adxcdy) + clift * (adxbdy - bdxady);
        permanent = (Absolute(bdxcdy) + Absolute(
                                            cdxbdy
                                        )) * alift + (Absolute(cdxady) + Absolute(
                                                                             adxcdy
                                                                         )) * blift + (Absolute(
                                                                                           adxbdy
                                                                                       ) + Absolute(
                                                                                               bdxady
                                                                                           )) * clift;
        errbound = self.iccerrboundA * permanent;
        if det > errbound || -det > errbound {
            det
        } else {
            self.incircleadapt(pa,pb,pc,pd,permanent)
        }
    }

    
    pub unsafe fn insphereexact(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64,
        mut pe : *const f64
    ) -> f64 {
        let mut axby1 : f64;
        let mut bxcy1 : f64;
        let mut cxdy1 : f64;
        let mut dxey1 : f64;
        let mut exay1 : f64;
        let mut bxay1 : f64;
        let mut cxby1 : f64;
        let mut dxcy1 : f64;
        let mut exdy1 : f64;
        let mut axey1 : f64;
        let mut axcy1 : f64;
        let mut bxdy1 : f64;
        let mut cxey1 : f64;
        let mut dxay1 : f64;
        let mut exby1 : f64;
        let mut cxay1 : f64;
        let mut dxby1 : f64;
        let mut excy1 : f64;
        let mut axdy1 : f64;
        let mut bxey1 : f64;
        let mut axby0 : f64;
        let mut bxcy0 : f64;
        let mut cxdy0 : f64;
        let mut dxey0 : f64;
        let mut exay0 : f64;
        let mut bxay0 : f64;
        let mut cxby0 : f64;
        let mut dxcy0 : f64;
        let mut exdy0 : f64;
        let mut axey0 : f64;
        let mut axcy0 : f64;
        let mut bxdy0 : f64;
        let mut cxey0 : f64;
        let mut dxay0 : f64;
        let mut exby0 : f64;
        let mut cxay0 : f64;
        let mut dxby0 : f64;
        let mut excy0 : f64;
        let mut axdy0 : f64;
        let mut bxey0 : f64;
        let mut ab : [f64; 4];
        let mut bc : [f64; 4];
        let mut cd : [f64; 4];
        let mut de : [f64; 4];
        let mut ea : [f64; 4];
        let mut ac : [f64; 4];
        let mut bd : [f64; 4];
        let mut ce : [f64; 4];
        let mut da : [f64; 4];
        let mut eb : [f64; 4];
        let mut temp8a : [f64; 8];
        let mut temp8b : [f64; 8];
        let mut temp16 : [f64; 16];
        let mut temp8alen : i32;
        let mut temp8blen : i32;
        let mut temp16len : i32;
        let mut abc : [f64; 24];
        let mut bcd : [f64; 24];
        let mut cde : [f64; 24];
        let mut dea : [f64; 24];
        let mut eab : [f64; 24];
        let mut abd : [f64; 24];
        let mut bce : [f64; 24];
        let mut cda : [f64; 24];
        let mut deb : [f64; 24];
        let mut eac : [f64; 24];
        let mut abclen : i32;
        let mut bcdlen : i32;
        let mut cdelen : i32;
        let mut dealen : i32;
        let mut eablen : i32;
        let mut abdlen : i32;
        let mut bcelen : i32;
        let mut cdalen : i32;
        let mut deblen : i32;
        let mut eaclen : i32;
        let mut temp48a : [f64; 48];
        let mut temp48b : [f64; 48];
        let mut temp48alen : i32;
        let mut temp48blen : i32;
        let mut abcd : [f64; 96];
        let mut bcde : [f64; 96];
        let mut cdea : [f64; 96];
        let mut deab : [f64; 96];
        let mut eabc : [f64; 96];
        let mut abcdlen : i32;
        let mut bcdelen : i32;
        let mut cdealen : i32;
        let mut deablen : i32;
        let mut eabclen : i32;
        let mut temp192 : [f64; 192];
        let mut det384x : [f64; 384];
        let mut det384y : [f64; 384];
        let mut det384z : [f64; 384];
        let mut xlen : i32;
        let mut ylen : i32;
        let mut zlen : i32;
        let mut detxy : [f64; 768];
        let mut xylen : i32;
        let mut adet : [f64; 1152];
        let mut bdet : [f64; 1152];
        let mut cdet : [f64; 1152];
        let mut ddet : [f64; 1152];
        let mut edet : [f64; 1152];
        let mut alen : i32;
        let mut blen : i32;
        let mut clen : i32;
        let mut dlen : i32;
        let mut elen : i32;
        let mut abdet : [f64; 2304];
        let mut cddet : [f64; 2304];
        let mut cdedet : [f64; 3456];
        let mut ablen : i32;
        let mut cdlen : i32;
        let mut deter : [f64; 5760];
        let mut deterlen : i32;
        let mut i : i32;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _0 : f64;

        axby1 = uninitialized();
        bxcy1 = uninitialized();
        cxdy1 = uninitialized();
        dxey1 = uninitialized();
        exay1 = uninitialized();
        bxay1 = uninitialized();
        cxby1 = uninitialized();
        dxcy1 = uninitialized();
        exdy1 = uninitialized();
        axey1 = uninitialized();
        axcy1 = uninitialized();
        bxdy1 = uninitialized();
        cxey1 = uninitialized();
        dxay1 = uninitialized();
        exby1 = uninitialized();
        cxay1 = uninitialized();
        dxby1 = uninitialized();
        excy1 = uninitialized();
        axdy1 = uninitialized();
        bxey1 = uninitialized();
        axby0 = uninitialized();
        bxcy0 = uninitialized();
        cxdy0 = uninitialized();
        dxey0 = uninitialized();
        exay0 = uninitialized();
        bxay0 = uninitialized();
        cxby0 = uninitialized();
        dxcy0 = uninitialized();
        exdy0 = uninitialized();
        axey0 = uninitialized();
        axcy0 = uninitialized();
        bxdy0 = uninitialized();
        cxey0 = uninitialized();
        dxay0 = uninitialized();
        exby0 = uninitialized();
        cxay0 = uninitialized();
        dxby0 = uninitialized();
        excy0 = uninitialized();
        axdy0 = uninitialized();
        bxey0 = uninitialized();
        ab = uninitialized();
        bc = uninitialized();
        cd = uninitialized();
        de = uninitialized();
        ea = uninitialized();
        ac = uninitialized();
        bd = uninitialized();
        ce = uninitialized();
        da = uninitialized();
        eb = uninitialized();
        temp8a = uninitialized();
        temp8b = uninitialized();
        temp16 = uninitialized();
        temp8alen = uninitialized();
        temp8blen = uninitialized();
        temp16len = uninitialized();
        abc = uninitialized();
        bcd = uninitialized();
        cde = uninitialized();
        dea = uninitialized();
        eab = uninitialized();
        abd = uninitialized();
        bce = uninitialized();
        cda = uninitialized();
        deb = uninitialized();
        eac = uninitialized();
        abclen = uninitialized();
        bcdlen = uninitialized();
        cdelen = uninitialized();
        dealen = uninitialized();
        eablen = uninitialized();
        abdlen = uninitialized();
        bcelen = uninitialized();
        cdalen = uninitialized();
        deblen = uninitialized();
        eaclen = uninitialized();
        temp48a = uninitialized();
        temp48b = uninitialized();
        temp48alen = uninitialized();
        temp48blen = uninitialized();
        abcd = uninitialized();
        bcde = uninitialized();
        cdea = uninitialized();
        deab = uninitialized();
        eabc = uninitialized();
        abcdlen = uninitialized();
        bcdelen = uninitialized();
        cdealen = uninitialized();
        deablen = uninitialized();
        eabclen = uninitialized();
        temp192 = uninitialized();
        det384x = uninitialized();
        det384y = uninitialized();
        det384z = uninitialized();
        xlen = uninitialized();
        ylen = uninitialized();
        zlen = uninitialized();
        detxy = uninitialized();
        xylen = uninitialized();
        adet = uninitialized();
        bdet = uninitialized();
        cdet = uninitialized();
        ddet = uninitialized();
        edet = uninitialized();
        alen = uninitialized();
        blen = uninitialized();
        clen = uninitialized();
        dlen = uninitialized();
        elen = uninitialized();
        abdet = uninitialized();
        cddet = uninitialized();
        cdedet = uninitialized();
        ablen = uninitialized();
        cdlen = uninitialized();
        deter = uninitialized();
        deterlen = uninitialized();
        i = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        ahi = uninitialized();
        alo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _0 = uninitialized();
        axby1 = *pa.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = axby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axby0 = alo * blo - err3;
        bxay1 = *pb.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = bxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxay0 = alo * blo - err3;
        _i = axby0 - bxay0;
        bvirt = axby0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxay0;
        around = axby0 - avirt;
        ab[0usize] = around + bround;
        _j = axby1 + _i;
        bvirt = _j - axby1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = axby1 - avirt;
        _0 = around + bround;
        _i = _0 - bxay1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxay1;
        around = _0 - avirt;
        ab[1usize] = around + bround;
        ab[3usize] = _j + _i;
        bvirt = ab[3usize] - _j;
        avirt = ab[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ab[2usize] = around + bround;
        bxcy1 = *pb.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = bxcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxcy0 = alo * blo - err3;
        cxby1 = *pc.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = cxby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxby0 = alo * blo - err3;
        _i = bxcy0 - cxby0;
        bvirt = bxcy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxby0;
        around = bxcy0 - avirt;
        bc[0usize] = around + bround;
        _j = bxcy1 + _i;
        bvirt = _j - bxcy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bxcy1 - avirt;
        _0 = around + bround;
        _i = _0 - cxby1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxby1;
        around = _0 - avirt;
        bc[1usize] = around + bround;
        bc[3usize] = _j + _i;
        bvirt = bc[3usize] - _j;
        avirt = bc[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bc[2usize] = around + bround;
        cxdy1 = *pc.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = cxdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxdy0 = alo * blo - err3;
        dxcy1 = *pd.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = dxcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxcy0 = alo * blo - err3;
        _i = cxdy0 - dxcy0;
        bvirt = cxdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxcy0;
        around = cxdy0 - avirt;
        cd[0usize] = around + bround;
        _j = cxdy1 + _i;
        bvirt = _j - cxdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = cxdy1 - avirt;
        _0 = around + bround;
        _i = _0 - dxcy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxcy1;
        around = _0 - avirt;
        cd[1usize] = around + bround;
        cd[3usize] = _j + _i;
        bvirt = cd[3usize] - _j;
        avirt = cd[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        cd[2usize] = around + bround;
        dxey1 = *pd.offset(0isize) * *pe.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pe.offset(1isize);
        abig = c - *pe.offset(1isize);
        bhi = c - abig;
        blo = *pe.offset(1isize) - bhi;
        err1 = dxey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxey0 = alo * blo - err3;
        exdy1 = *pe.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pe.offset(0isize);
        abig = c - *pe.offset(0isize);
        ahi = c - abig;
        alo = *pe.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = exdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        exdy0 = alo * blo - err3;
        _i = dxey0 - exdy0;
        bvirt = dxey0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - exdy0;
        around = dxey0 - avirt;
        de[0usize] = around + bround;
        _j = dxey1 + _i;
        bvirt = _j - dxey1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = dxey1 - avirt;
        _0 = around + bround;
        _i = _0 - exdy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - exdy1;
        around = _0 - avirt;
        de[1usize] = around + bround;
        de[3usize] = _j + _i;
        bvirt = de[3usize] - _j;
        avirt = de[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        de[2usize] = around + bround;
        exay1 = *pe.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pe.offset(0isize);
        abig = c - *pe.offset(0isize);
        ahi = c - abig;
        alo = *pe.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = exay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        exay0 = alo * blo - err3;
        axey1 = *pa.offset(0isize) * *pe.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pe.offset(1isize);
        abig = c - *pe.offset(1isize);
        bhi = c - abig;
        blo = *pe.offset(1isize) - bhi;
        err1 = axey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axey0 = alo * blo - err3;
        _i = exay0 - axey0;
        bvirt = exay0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axey0;
        around = exay0 - avirt;
        ea[0usize] = around + bround;
        _j = exay1 + _i;
        bvirt = _j - exay1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = exay1 - avirt;
        _0 = around + bround;
        _i = _0 - axey1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axey1;
        around = _0 - avirt;
        ea[1usize] = around + bround;
        ea[3usize] = _j + _i;
        bvirt = ea[3usize] - _j;
        avirt = ea[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ea[2usize] = around + bround;
        axcy1 = *pa.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = axcy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axcy0 = alo * blo - err3;
        cxay1 = *pc.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = cxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxay0 = alo * blo - err3;
        _i = axcy0 - cxay0;
        bvirt = axcy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxay0;
        around = axcy0 - avirt;
        ac[0usize] = around + bround;
        _j = axcy1 + _i;
        bvirt = _j - axcy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = axcy1 - avirt;
        _0 = around + bround;
        _i = _0 - cxay1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cxay1;
        around = _0 - avirt;
        ac[1usize] = around + bround;
        ac[3usize] = _j + _i;
        bvirt = ac[3usize] - _j;
        avirt = ac[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ac[2usize] = around + bround;
        bxdy1 = *pb.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = bxdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxdy0 = alo * blo - err3;
        dxby1 = *pd.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = dxby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxby0 = alo * blo - err3;
        _i = bxdy0 - dxby0;
        bvirt = bxdy0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxby0;
        around = bxdy0 - avirt;
        bd[0usize] = around + bround;
        _j = bxdy1 + _i;
        bvirt = _j - bxdy1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bxdy1 - avirt;
        _0 = around + bround;
        _i = _0 - dxby1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dxby1;
        around = _0 - avirt;
        bd[1usize] = around + bround;
        bd[3usize] = _j + _i;
        bvirt = bd[3usize] - _j;
        avirt = bd[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bd[2usize] = around + bround;
        cxey1 = *pc.offset(0isize) * *pe.offset(1isize);
        c = self.splitter * *pc.offset(0isize);
        abig = c - *pc.offset(0isize);
        ahi = c - abig;
        alo = *pc.offset(0isize) - ahi;
        c = self.splitter * *pe.offset(1isize);
        abig = c - *pe.offset(1isize);
        bhi = c - abig;
        blo = *pe.offset(1isize) - bhi;
        err1 = cxey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cxey0 = alo * blo - err3;
        excy1 = *pe.offset(0isize) * *pc.offset(1isize);
        c = self.splitter * *pe.offset(0isize);
        abig = c - *pe.offset(0isize);
        ahi = c - abig;
        alo = *pe.offset(0isize) - ahi;
        c = self.splitter * *pc.offset(1isize);
        abig = c - *pc.offset(1isize);
        bhi = c - abig;
        blo = *pc.offset(1isize) - bhi;
        err1 = excy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        excy0 = alo * blo - err3;
        _i = cxey0 - excy0;
        bvirt = cxey0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - excy0;
        around = cxey0 - avirt;
        ce[0usize] = around + bround;
        _j = cxey1 + _i;
        bvirt = _j - cxey1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = cxey1 - avirt;
        _0 = around + bround;
        _i = _0 - excy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - excy1;
        around = _0 - avirt;
        ce[1usize] = around + bround;
        ce[3usize] = _j + _i;
        bvirt = ce[3usize] - _j;
        avirt = ce[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ce[2usize] = around + bround;
        dxay1 = *pd.offset(0isize) * *pa.offset(1isize);
        c = self.splitter * *pd.offset(0isize);
        abig = c - *pd.offset(0isize);
        ahi = c - abig;
        alo = *pd.offset(0isize) - ahi;
        c = self.splitter * *pa.offset(1isize);
        abig = c - *pa.offset(1isize);
        bhi = c - abig;
        blo = *pa.offset(1isize) - bhi;
        err1 = dxay1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dxay0 = alo * blo - err3;
        axdy1 = *pa.offset(0isize) * *pd.offset(1isize);
        c = self.splitter * *pa.offset(0isize);
        abig = c - *pa.offset(0isize);
        ahi = c - abig;
        alo = *pa.offset(0isize) - ahi;
        c = self.splitter * *pd.offset(1isize);
        abig = c - *pd.offset(1isize);
        bhi = c - abig;
        blo = *pd.offset(1isize) - bhi;
        err1 = axdy1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        axdy0 = alo * blo - err3;
        _i = dxay0 - axdy0;
        bvirt = dxay0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axdy0;
        around = dxay0 - avirt;
        da[0usize] = around + bround;
        _j = dxay1 + _i;
        bvirt = _j - dxay1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = dxay1 - avirt;
        _0 = around + bround;
        _i = _0 - axdy1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - axdy1;
        around = _0 - avirt;
        da[1usize] = around + bround;
        da[3usize] = _j + _i;
        bvirt = da[3usize] - _j;
        avirt = da[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        da[2usize] = around + bround;
        exby1 = *pe.offset(0isize) * *pb.offset(1isize);
        c = self.splitter * *pe.offset(0isize);
        abig = c - *pe.offset(0isize);
        ahi = c - abig;
        alo = *pe.offset(0isize) - ahi;
        c = self.splitter * *pb.offset(1isize);
        abig = c - *pb.offset(1isize);
        bhi = c - abig;
        blo = *pb.offset(1isize) - bhi;
        err1 = exby1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        exby0 = alo * blo - err3;
        bxey1 = *pb.offset(0isize) * *pe.offset(1isize);
        c = self.splitter * *pb.offset(0isize);
        abig = c - *pb.offset(0isize);
        ahi = c - abig;
        alo = *pb.offset(0isize) - ahi;
        c = self.splitter * *pe.offset(1isize);
        abig = c - *pe.offset(1isize);
        bhi = c - abig;
        blo = *pe.offset(1isize) - bhi;
        err1 = bxey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bxey0 = alo * blo - err3;
        _i = exby0 - bxey0;
        bvirt = exby0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxey0;
        around = exby0 - avirt;
        eb[0usize] = around + bround;
        _j = exby1 + _i;
        bvirt = _j - exby1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = exby1 - avirt;
        _0 = around + bround;
        _i = _0 - bxey1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bxey1;
        around = _0 - avirt;
        eb[1usize] = around + bround;
        eb[3usize] = _j + _i;
        bvirt = eb[3usize] - _j;
        avirt = eb[3usize] - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        eb[2usize] = around + bround;
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        bc.as_ptr(),
                        *pa.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        ac.as_ptr(),
                        -*pb.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ab.as_ptr(),
                        *pc.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        abclen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     abc.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        cd.as_ptr(),
                        *pb.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        bd.as_ptr(),
                        -*pc.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        bc.as_ptr(),
                        *pd.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        bcdlen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     bcd.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        de.as_ptr(),
                        *pc.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        ce.as_ptr(),
                        -*pd.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        cd.as_ptr(),
                        *pe.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        cdelen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     cde.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ea.as_ptr(),
                        *pd.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        da.as_ptr(),
                        -*pe.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        de.as_ptr(),
                        *pa.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        dealen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     dea.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ab.as_ptr(),
                        *pe.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        eb.as_ptr(),
                        -*pa.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ea.as_ptr(),
                        *pb.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        eablen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     eab.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        bd.as_ptr(),
                        *pa.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        da.as_ptr(),
                        *pb.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ab.as_ptr(),
                        *pd.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        abdlen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     abd.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ce.as_ptr(),
                        *pb.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        eb.as_ptr(),
                        *pc.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        bc.as_ptr(),
                        *pe.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        bcelen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     bce.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        da.as_ptr(),
                        *pc.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        ac.as_ptr(),
                        *pd.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        cd.as_ptr(),
                        *pa.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        cdalen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     cda.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        eb.as_ptr(),
                        *pd.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        bd.as_ptr(),
                        *pe.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        de.as_ptr(),
                        *pb.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        deblen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     deb.as_mut_ptr()
                 );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ac.as_ptr(),
                        *pe.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        ce.as_ptr(),
                        *pa.offset(2isize),
                        temp8b.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ea.as_ptr(),
                        *pc.offset(2isize),
                        temp8a.as_mut_ptr()
                    );
        eaclen = fast_expansion_sum_zeroelim(
                     temp8alen,
                     temp8a.as_ptr(),
                     temp16len,
                     temp16.as_ptr(),
                     eac.as_mut_ptr()
                 );
        temp48alen = fast_expansion_sum_zeroelim(
                         cdelen,
                         cde.as_ptr(),
                         bcelen,
                         bce.as_ptr(),
                         temp48a.as_mut_ptr()
                     );
        temp48blen = fast_expansion_sum_zeroelim(
                         deblen,
                         deb.as_ptr(),
                         bcdlen,
                         bcd.as_ptr(),
                         temp48b.as_mut_ptr()
                     );
        i = 0i32;
        'loop1: loop {
            if !(i < temp48blen) {
                break;
            }
            temp48b[i as (usize)] = -temp48b[i as (usize)];
            i = i + 1;
        }
        bcdelen = fast_expansion_sum_zeroelim(
                      temp48alen,
                      temp48a.as_ptr(),
                      temp48blen,
                      temp48b.as_ptr(),
                      bcde.as_mut_ptr()
                  );
        xlen = self.scale_expansion_zeroelim(
                   bcdelen,
                   bcde.as_ptr(),
                   *pa.offset(0isize),
                   temp192.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   temp192.as_ptr(),
                   *pa.offset(0isize),
                   det384x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   bcdelen,
                   bcde.as_ptr(),
                   *pa.offset(1isize),
                   temp192.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   temp192.as_ptr(),
                   *pa.offset(1isize),
                   det384y.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   bcdelen,
                   bcde.as_ptr(),
                   *pa.offset(2isize),
                   temp192.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   zlen,
                   temp192.as_ptr(),
                   *pa.offset(2isize),
                   det384z.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    det384x.as_ptr(),
                    ylen,
                    det384y.as_ptr(),
                    detxy.as_mut_ptr()
                );
        alen = fast_expansion_sum_zeroelim(
                   xylen,
                   detxy.as_ptr(),
                   zlen,
                   det384z.as_ptr(),
                   adet.as_mut_ptr()
               );
        temp48alen = fast_expansion_sum_zeroelim(
                         dealen,
                         dea.as_ptr(),
                         cdalen,
                         cda.as_ptr(),
                         temp48a.as_mut_ptr()
                     );
        temp48blen = fast_expansion_sum_zeroelim(
                         eaclen,
                         eac.as_ptr(),
                         cdelen,
                         cde.as_ptr(),
                         temp48b.as_mut_ptr()
                     );
        i = 0i32;
        'loop3: loop {
            if !(i < temp48blen) {
                break;
            }
            temp48b[i as (usize)] = -temp48b[i as (usize)];
            i = i + 1;
        }
        cdealen = fast_expansion_sum_zeroelim(
                      temp48alen,
                      temp48a.as_ptr(),
                      temp48blen,
                      temp48b.as_ptr(),
                      cdea.as_mut_ptr()
                  );
        xlen = self.scale_expansion_zeroelim(
                   cdealen,
                   cdea.as_ptr(),
                   *pb.offset(0isize),
                   temp192.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   temp192.as_ptr(),
                   *pb.offset(0isize),
                   det384x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   cdealen,
                   cdea.as_ptr(),
                   *pb.offset(1isize),
                   temp192.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   temp192.as_ptr(),
                   *pb.offset(1isize),
                   det384y.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   cdealen,
                   cdea.as_ptr(),
                   *pb.offset(2isize),
                   temp192.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   zlen,
                   temp192.as_ptr(),
                   *pb.offset(2isize),
                   det384z.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    det384x.as_ptr(),
                    ylen,
                    det384y.as_ptr(),
                    detxy.as_mut_ptr()
                );
        blen = fast_expansion_sum_zeroelim(
                   xylen,
                   detxy.as_ptr(),
                   zlen,
                   det384z.as_ptr(),
                   bdet.as_mut_ptr()
               );
        temp48alen = fast_expansion_sum_zeroelim(
                         eablen,
                         eab.as_ptr(),
                         deblen,
                         deb.as_ptr(),
                         temp48a.as_mut_ptr()
                     );
        temp48blen = fast_expansion_sum_zeroelim(
                         abdlen,
                         abd.as_ptr(),
                         dealen,
                         dea.as_ptr(),
                         temp48b.as_mut_ptr()
                     );
        i = 0i32;
        'loop5: loop {
            if !(i < temp48blen) {
                break;
            }
            temp48b[i as (usize)] = -temp48b[i as (usize)];
            i = i + 1;
        }
        deablen = fast_expansion_sum_zeroelim(
                      temp48alen,
                      temp48a.as_ptr(),
                      temp48blen,
                      temp48b.as_ptr(),
                      deab.as_mut_ptr()
                  );
        xlen = self.scale_expansion_zeroelim(
                   deablen,
                   deab.as_ptr(),
                   *pc.offset(0isize),
                   temp192.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   temp192.as_ptr(),
                   *pc.offset(0isize),
                   det384x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   deablen,
                   deab.as_ptr(),
                   *pc.offset(1isize),
                   temp192.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   temp192.as_ptr(),
                   *pc.offset(1isize),
                   det384y.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   deablen,
                   deab.as_ptr(),
                   *pc.offset(2isize),
                   temp192.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   zlen,
                   temp192.as_ptr(),
                   *pc.offset(2isize),
                   det384z.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    det384x.as_ptr(),
                    ylen,
                    det384y.as_ptr(),
                    detxy.as_mut_ptr()
                );
        clen = fast_expansion_sum_zeroelim(
                   xylen,
                   detxy.as_ptr(),
                   zlen,
                   det384z.as_ptr(),
                   cdet.as_mut_ptr()
               );
        temp48alen = fast_expansion_sum_zeroelim(
                         abclen,
                         abc.as_ptr(),
                         eaclen,
                         eac.as_ptr(),
                         temp48a.as_mut_ptr()
                     );
        temp48blen = fast_expansion_sum_zeroelim(
                         bcelen,
                         bce.as_ptr(),
                         eablen,
                         eab.as_ptr(),
                         temp48b.as_mut_ptr()
                     );
        i = 0i32;
        'loop7: loop {
            if !(i < temp48blen) {
                break;
            }
            temp48b[i as (usize)] = -temp48b[i as (usize)];
            i = i + 1;
        }
        eabclen = fast_expansion_sum_zeroelim(
                      temp48alen,
                      temp48a.as_ptr(),
                      temp48blen,
                      temp48b.as_ptr(),
                      eabc.as_mut_ptr()
                  );
        xlen = self.scale_expansion_zeroelim(
                   eabclen,
                   eabc.as_ptr(),
                   *pd.offset(0isize),
                   temp192.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   temp192.as_ptr(),
                   *pd.offset(0isize),
                   det384x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   eabclen,
                   eabc.as_ptr(),
                   *pd.offset(1isize),
                   temp192.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   temp192.as_ptr(),
                   *pd.offset(1isize),
                   det384y.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   eabclen,
                   eabc.as_ptr(),
                   *pd.offset(2isize),
                   temp192.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   zlen,
                   temp192.as_ptr(),
                   *pd.offset(2isize),
                   det384z.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    det384x.as_ptr(),
                    ylen,
                    det384y.as_ptr(),
                    detxy.as_mut_ptr()
                );
        dlen = fast_expansion_sum_zeroelim(
                   xylen,
                   detxy.as_ptr(),
                   zlen,
                   det384z.as_ptr(),
                   ddet.as_mut_ptr()
               );
        temp48alen = fast_expansion_sum_zeroelim(
                         bcdlen,
                         bcd.as_ptr(),
                         abdlen,
                         abd.as_ptr(),
                         temp48a.as_mut_ptr()
                     );
        temp48blen = fast_expansion_sum_zeroelim(
                         cdalen,
                         cda.as_ptr(),
                         abclen,
                         abc.as_ptr(),
                         temp48b.as_mut_ptr()
                     );
        i = 0i32;
        'loop9: loop {
            if !(i < temp48blen) {
                break;
            }
            temp48b[i as (usize)] = -temp48b[i as (usize)];
            i = i + 1;
        }
        abcdlen = fast_expansion_sum_zeroelim(
                      temp48alen,
                      temp48a.as_ptr(),
                      temp48blen,
                      temp48b.as_ptr(),
                      abcd.as_mut_ptr()
                  );
        xlen = self.scale_expansion_zeroelim(
                   abcdlen,
                   abcd.as_ptr(),
                   *pe.offset(0isize),
                   temp192.as_mut_ptr()
               );
        xlen = self.scale_expansion_zeroelim(
                   xlen,
                   temp192.as_ptr(),
                   *pe.offset(0isize),
                   det384x.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   abcdlen,
                   abcd.as_ptr(),
                   *pe.offset(1isize),
                   temp192.as_mut_ptr()
               );
        ylen = self.scale_expansion_zeroelim(
                   ylen,
                   temp192.as_ptr(),
                   *pe.offset(1isize),
                   det384y.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   abcdlen,
                   abcd.as_ptr(),
                   *pe.offset(2isize),
                   temp192.as_mut_ptr()
               );
        zlen = self.scale_expansion_zeroelim(
                   zlen,
                   temp192.as_ptr(),
                   *pe.offset(2isize),
                   det384z.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    det384x.as_ptr(),
                    ylen,
                    det384y.as_ptr(),
                    detxy.as_mut_ptr()
                );
        elen = fast_expansion_sum_zeroelim(
                   xylen,
                   detxy.as_ptr(),
                   zlen,
                   det384z.as_ptr(),
                   edet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        cdlen = fast_expansion_sum_zeroelim(
                    clen,
                    cdet.as_ptr(),
                    dlen,
                    ddet.as_ptr(),
                    cddet.as_mut_ptr()
                );
        cdelen = fast_expansion_sum_zeroelim(
                     cdlen,
                     cddet.as_ptr(),
                     elen,
                     edet.as_ptr(),
                     cdedet.as_mut_ptr()
                 );
        deterlen = fast_expansion_sum_zeroelim(
                       ablen,
                       abdet.as_ptr(),
                       cdelen,
                       cdedet.as_ptr(),
                       deter.as_mut_ptr()
                   );
        deter[(deterlen - 1i32) as (usize)]
    }

    
    pub unsafe fn insphereslow(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64,
        mut pe : *const f64
    ) -> f64 {
        let mut aex : f64;
        let mut bex : f64;
        let mut cex : f64;
        let mut dex : f64;
        let mut aey : f64;
        let mut bey : f64;
        let mut cey : f64;
        let mut dey : f64;
        let mut aez : f64;
        let mut bez : f64;
        let mut cez : f64;
        let mut dez : f64;
        let mut aextail : f64;
        let mut bextail : f64;
        let mut cextail : f64;
        let mut dextail : f64;
        let mut aeytail : f64;
        let mut beytail : f64;
        let mut ceytail : f64;
        let mut deytail : f64;
        let mut aeztail : f64;
        let mut beztail : f64;
        let mut ceztail : f64;
        let mut deztail : f64;
        let mut negate : f64;
        let mut negatetail : f64;
        let mut axby7 : f64;
        let mut bxcy7 : f64;
        let mut cxdy7 : f64;
        let mut dxay7 : f64;
        let mut axcy7 : f64;
        let mut bxdy7 : f64;
        let mut bxay7 : f64;
        let mut cxby7 : f64;
        let mut dxcy7 : f64;
        let mut axdy7 : f64;
        let mut cxay7 : f64;
        let mut dxby7 : f64;
        let mut axby : [f64; 8];
        let mut bxcy : [f64; 8];
        let mut cxdy : [f64; 8];
        let mut dxay : [f64; 8];
        let mut axcy : [f64; 8];
        let mut bxdy : [f64; 8];
        let mut bxay : [f64; 8];
        let mut cxby : [f64; 8];
        let mut dxcy : [f64; 8];
        let mut axdy : [f64; 8];
        let mut cxay : [f64; 8];
        let mut dxby : [f64; 8];
        let mut ab : [f64; 16];
        let mut bc : [f64; 16];
        let mut cd : [f64; 16];
        let mut da : [f64; 16];
        let mut ac : [f64; 16];
        let mut bd : [f64; 16];
        let mut ablen : i32;
        let mut bclen : i32;
        let mut cdlen : i32;
        let mut dalen : i32;
        let mut aclen : i32;
        let mut bdlen : i32;
        let mut temp32a : [f64; 32];
        let mut temp32b : [f64; 32];
        let mut temp64a : [f64; 64];
        let mut temp64b : [f64; 64];
        let mut temp64c : [f64; 64];
        let mut temp32alen : i32;
        let mut temp32blen : i32;
        let mut temp64alen : i32;
        let mut temp64blen : i32;
        let mut temp64clen : i32;
        let mut temp128 : [f64; 128];
        let mut temp192 : [f64; 192];
        let mut temp128len : i32;
        let mut temp192len : i32;
        let mut detx : [f64; 384];
        let mut detxx : [f64; 768];
        let mut detxt : [f64; 384];
        let mut detxxt : [f64; 768];
        let mut detxtxt : [f64; 768];
        let mut xlen : i32;
        let mut xxlen : i32;
        let mut xtlen : i32;
        let mut xxtlen : i32;
        let mut xtxtlen : i32;
        let mut x1 : [f64; 1536];
        let mut x2 : [f64; 2304];
        let mut x1len : i32;
        let mut x2len : i32;
        let mut dety : [f64; 384];
        let mut detyy : [f64; 768];
        let mut detyt : [f64; 384];
        let mut detyyt : [f64; 768];
        let mut detytyt : [f64; 768];
        let mut ylen : i32;
        let mut yylen : i32;
        let mut ytlen : i32;
        let mut yytlen : i32;
        let mut ytytlen : i32;
        let mut y1 : [f64; 1536];
        let mut y2 : [f64; 2304];
        let mut y1len : i32;
        let mut y2len : i32;
        let mut detz : [f64; 384];
        let mut detzz : [f64; 768];
        let mut detzt : [f64; 384];
        let mut detzzt : [f64; 768];
        let mut detztzt : [f64; 768];
        let mut zlen : i32;
        let mut zzlen : i32;
        let mut ztlen : i32;
        let mut zztlen : i32;
        let mut ztztlen : i32;
        let mut z1 : [f64; 1536];
        let mut z2 : [f64; 2304];
        let mut z1len : i32;
        let mut z2len : i32;
        let mut detxy : [f64; 4608];
        let mut xylen : i32;
        let mut adet : [f64; 6912];
        let mut bdet : [f64; 6912];
        let mut cdet : [f64; 6912];
        let mut ddet : [f64; 6912];
        let mut alen : i32;
        let mut blen : i32;
        let mut clen : i32;
        let mut dlen : i32;
        let mut abdet : [f64; 13824];
        let mut cddet : [f64; 13824];
        let mut deter : [f64; 27648];
        let mut deterlen : i32;
        let mut i : i32;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut a0hi : f64;
        let mut a0lo : f64;
        let mut a1hi : f64;
        let mut a1lo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _k : f64;
        let mut _l : f64;
        let mut _m : f64;
        let mut _n : f64;
        let mut _0 : f64;
        let mut _1 : f64;
        let mut _2 : f64;

        aex = uninitialized();
        bex = uninitialized();
        cex = uninitialized();
        dex = uninitialized();
        aey = uninitialized();
        bey = uninitialized();
        cey = uninitialized();
        dey = uninitialized();
        aez = uninitialized();
        bez = uninitialized();
        cez = uninitialized();
        dez = uninitialized();
        aextail = uninitialized();
        bextail = uninitialized();
        cextail = uninitialized();
        dextail = uninitialized();
        aeytail = uninitialized();
        beytail = uninitialized();
        ceytail = uninitialized();
        deytail = uninitialized();
        aeztail = uninitialized();
        beztail = uninitialized();
        ceztail = uninitialized();
        deztail = uninitialized();
        negate = uninitialized();
        negatetail = uninitialized();
        axby7 = uninitialized();
        bxcy7 = uninitialized();
        cxdy7 = uninitialized();
        dxay7 = uninitialized();
        axcy7 = uninitialized();
        bxdy7 = uninitialized();
        bxay7 = uninitialized();
        cxby7 = uninitialized();
        dxcy7 = uninitialized();
        axdy7 = uninitialized();
        cxay7 = uninitialized();
        dxby7 = uninitialized();
        axby = uninitialized();
        bxcy = uninitialized();
        cxdy = uninitialized();
        dxay = uninitialized();
        axcy = uninitialized();
        bxdy = uninitialized();
        bxay = uninitialized();
        cxby = uninitialized();
        dxcy = uninitialized();
        axdy = uninitialized();
        cxay = uninitialized();
        dxby = uninitialized();
        ab = uninitialized();
        bc = uninitialized();
        cd = uninitialized();
        da = uninitialized();
        ac = uninitialized();
        bd = uninitialized();
        ablen = uninitialized();
        bclen = uninitialized();
        cdlen = uninitialized();
        dalen = uninitialized();
        aclen = uninitialized();
        bdlen = uninitialized();
        temp32a = uninitialized();
        temp32b = uninitialized();
        temp64a = uninitialized();
        temp64b = uninitialized();
        temp64c = uninitialized();
        temp32alen = uninitialized();
        temp32blen = uninitialized();
        temp64alen = uninitialized();
        temp64blen = uninitialized();
        temp64clen = uninitialized();
        temp128 = uninitialized();
        temp192 = uninitialized();
        temp128len = uninitialized();
        temp192len = uninitialized();
        detx = uninitialized();
        detxx = uninitialized();
        detxt = uninitialized();
        detxxt = uninitialized();
        detxtxt = uninitialized();
        xlen = uninitialized();
        xxlen = uninitialized();
        xtlen = uninitialized();
        xxtlen = uninitialized();
        xtxtlen = uninitialized();
        x1 = uninitialized();
        x2 = uninitialized();
        x1len = uninitialized();
        x2len = uninitialized();
        dety = uninitialized();
        detyy = uninitialized();
        detyt = uninitialized();
        detyyt = uninitialized();
        detytyt = uninitialized();
        ylen = uninitialized();
        yylen = uninitialized();
        ytlen = uninitialized();
        yytlen = uninitialized();
        ytytlen = uninitialized();
        y1 = uninitialized();
        y2 = uninitialized();
        y1len = uninitialized();
        y2len = uninitialized();
        detz = uninitialized();
        detzz = uninitialized();
        detzt = uninitialized();
        detzzt = uninitialized();
        detztzt = uninitialized();
        zlen = uninitialized();
        zzlen = uninitialized();
        ztlen = uninitialized();
        zztlen = uninitialized();
        ztztlen = uninitialized();
        z1 = uninitialized();
        z2 = uninitialized();
        z1len = uninitialized();
        z2len = uninitialized();
        detxy = uninitialized();
        xylen = uninitialized();
        adet = uninitialized();
        bdet = uninitialized();
        cdet = uninitialized();
        ddet = uninitialized();
        alen = uninitialized();
        blen = uninitialized();
        clen = uninitialized();
        dlen = uninitialized();
        abdet = uninitialized();
        cddet = uninitialized();
        deter = uninitialized();
        deterlen = uninitialized();
        i = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        a0hi = uninitialized();
        a0lo = uninitialized();
        a1hi = uninitialized();
        a1lo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _k = uninitialized();
        _l = uninitialized();
        _m = uninitialized();
        _n = uninitialized();
        _0 = uninitialized();
        _1 = uninitialized();
        _2 = uninitialized();
        aex = *pa.offset(0isize) - *pe.offset(0isize);
        bvirt = *pa.offset(0isize) - aex;
        avirt = aex + bvirt;
        bround = bvirt - *pe.offset(0isize);
        around = *pa.offset(0isize) - avirt;
        aextail = around + bround;
        aey = *pa.offset(1isize) - *pe.offset(1isize);
        bvirt = *pa.offset(1isize) - aey;
        avirt = aey + bvirt;
        bround = bvirt - *pe.offset(1isize);
        around = *pa.offset(1isize) - avirt;
        aeytail = around + bround;
        aez = *pa.offset(2isize) - *pe.offset(2isize);
        bvirt = *pa.offset(2isize) - aez;
        avirt = aez + bvirt;
        bround = bvirt - *pe.offset(2isize);
        around = *pa.offset(2isize) - avirt;
        aeztail = around + bround;
        bex = *pb.offset(0isize) - *pe.offset(0isize);
        bvirt = *pb.offset(0isize) - bex;
        avirt = bex + bvirt;
        bround = bvirt - *pe.offset(0isize);
        around = *pb.offset(0isize) - avirt;
        bextail = around + bround;
        bey = *pb.offset(1isize) - *pe.offset(1isize);
        bvirt = *pb.offset(1isize) - bey;
        avirt = bey + bvirt;
        bround = bvirt - *pe.offset(1isize);
        around = *pb.offset(1isize) - avirt;
        beytail = around + bround;
        bez = *pb.offset(2isize) - *pe.offset(2isize);
        bvirt = *pb.offset(2isize) - bez;
        avirt = bez + bvirt;
        bround = bvirt - *pe.offset(2isize);
        around = *pb.offset(2isize) - avirt;
        beztail = around + bround;
        cex = *pc.offset(0isize) - *pe.offset(0isize);
        bvirt = *pc.offset(0isize) - cex;
        avirt = cex + bvirt;
        bround = bvirt - *pe.offset(0isize);
        around = *pc.offset(0isize) - avirt;
        cextail = around + bround;
        cey = *pc.offset(1isize) - *pe.offset(1isize);
        bvirt = *pc.offset(1isize) - cey;
        avirt = cey + bvirt;
        bround = bvirt - *pe.offset(1isize);
        around = *pc.offset(1isize) - avirt;
        ceytail = around + bround;
        cez = *pc.offset(2isize) - *pe.offset(2isize);
        bvirt = *pc.offset(2isize) - cez;
        avirt = cez + bvirt;
        bround = bvirt - *pe.offset(2isize);
        around = *pc.offset(2isize) - avirt;
        ceztail = around + bround;
        dex = *pd.offset(0isize) - *pe.offset(0isize);
        bvirt = *pd.offset(0isize) - dex;
        avirt = dex + bvirt;
        bround = bvirt - *pe.offset(0isize);
        around = *pd.offset(0isize) - avirt;
        dextail = around + bround;
        dey = *pd.offset(1isize) - *pe.offset(1isize);
        bvirt = *pd.offset(1isize) - dey;
        avirt = dey + bvirt;
        bround = bvirt - *pe.offset(1isize);
        around = *pd.offset(1isize) - avirt;
        deytail = around + bround;
        dez = *pd.offset(2isize) - *pe.offset(2isize);
        bvirt = *pd.offset(2isize) - dez;
        avirt = dez + bvirt;
        bround = bvirt - *pe.offset(2isize);
        around = *pd.offset(2isize) - avirt;
        deztail = around + bround;
        c = self.splitter * aextail;
        abig = c - aextail;
        a0hi = c - abig;
        a0lo = aextail - a0hi;
        c = self.splitter * beytail;
        abig = c - beytail;
        bhi = c - abig;
        blo = beytail - bhi;
        _i = aextail * beytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        axby[0usize] = a0lo * blo - err3;
        c = self.splitter * aex;
        abig = c - aex;
        a1hi = c - abig;
        a1lo = aex - a1hi;
        _j = aex * beytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * bey;
        abig = c - bey;
        bhi = c - abig;
        blo = bey - bhi;
        _i = aextail * bey;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = aex * bey;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axby[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        axby[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        axby[5usize] = around + bround;
        axby7 = _m + _k;
        bvirt = axby7 - _m;
        avirt = axby7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        axby[6usize] = around + bround;
        axby[7usize] = axby7;
        negate = -aey;
        negatetail = -aeytail;
        c = self.splitter * bextail;
        abig = c - bextail;
        a0hi = c - abig;
        a0lo = bextail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = bextail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        bxay[0usize] = a0lo * blo - err3;
        c = self.splitter * bex;
        abig = c - bex;
        a1hi = c - abig;
        a1lo = bex - a1hi;
        _j = bex * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = bextail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = bex * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxay[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        bxay[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        bxay[5usize] = around + bround;
        bxay7 = _m + _k;
        bvirt = bxay7 - _m;
        avirt = bxay7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        bxay[6usize] = around + bround;
        bxay[7usize] = bxay7;
        ablen = fast_expansion_sum_zeroelim(
                    8i32,
                    axby.as_ptr(),
                    8i32,
                    bxay.as_ptr(),
                    ab.as_mut_ptr()
                );
        c = self.splitter * bextail;
        abig = c - bextail;
        a0hi = c - abig;
        a0lo = bextail - a0hi;
        c = self.splitter * ceytail;
        abig = c - ceytail;
        bhi = c - abig;
        blo = ceytail - bhi;
        _i = bextail * ceytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        bxcy[0usize] = a0lo * blo - err3;
        c = self.splitter * bex;
        abig = c - bex;
        a1hi = c - abig;
        a1lo = bex - a1hi;
        _j = bex * ceytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * cey;
        abig = c - cey;
        bhi = c - abig;
        blo = cey - bhi;
        _i = bextail * cey;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = bex * cey;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxcy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        bxcy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        bxcy[5usize] = around + bround;
        bxcy7 = _m + _k;
        bvirt = bxcy7 - _m;
        avirt = bxcy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        bxcy[6usize] = around + bround;
        bxcy[7usize] = bxcy7;
        negate = -bey;
        negatetail = -beytail;
        c = self.splitter * cextail;
        abig = c - cextail;
        a0hi = c - abig;
        a0lo = cextail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = cextail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        cxby[0usize] = a0lo * blo - err3;
        c = self.splitter * cex;
        abig = c - cex;
        a1hi = c - abig;
        a1lo = cex - a1hi;
        _j = cex * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = cextail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = cex * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxby[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        cxby[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        cxby[5usize] = around + bround;
        cxby7 = _m + _k;
        bvirt = cxby7 - _m;
        avirt = cxby7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        cxby[6usize] = around + bround;
        cxby[7usize] = cxby7;
        bclen = fast_expansion_sum_zeroelim(
                    8i32,
                    bxcy.as_ptr(),
                    8i32,
                    cxby.as_ptr(),
                    bc.as_mut_ptr()
                );
        c = self.splitter * cextail;
        abig = c - cextail;
        a0hi = c - abig;
        a0lo = cextail - a0hi;
        c = self.splitter * deytail;
        abig = c - deytail;
        bhi = c - abig;
        blo = deytail - bhi;
        _i = cextail * deytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        cxdy[0usize] = a0lo * blo - err3;
        c = self.splitter * cex;
        abig = c - cex;
        a1hi = c - abig;
        a1lo = cex - a1hi;
        _j = cex * deytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * dey;
        abig = c - dey;
        bhi = c - abig;
        blo = dey - bhi;
        _i = cextail * dey;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxdy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = cex * dey;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxdy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxdy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        cxdy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        cxdy[5usize] = around + bround;
        cxdy7 = _m + _k;
        bvirt = cxdy7 - _m;
        avirt = cxdy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        cxdy[6usize] = around + bround;
        cxdy[7usize] = cxdy7;
        negate = -cey;
        negatetail = -ceytail;
        c = self.splitter * dextail;
        abig = c - dextail;
        a0hi = c - abig;
        a0lo = dextail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = dextail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        dxcy[0usize] = a0lo * blo - err3;
        c = self.splitter * dex;
        abig = c - dex;
        a1hi = c - abig;
        a1lo = dex - a1hi;
        _j = dex * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = dextail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxcy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = dex * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxcy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxcy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        dxcy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        dxcy[5usize] = around + bround;
        dxcy7 = _m + _k;
        bvirt = dxcy7 - _m;
        avirt = dxcy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        dxcy[6usize] = around + bround;
        dxcy[7usize] = dxcy7;
        cdlen = fast_expansion_sum_zeroelim(
                    8i32,
                    cxdy.as_ptr(),
                    8i32,
                    dxcy.as_ptr(),
                    cd.as_mut_ptr()
                );
        c = self.splitter * dextail;
        abig = c - dextail;
        a0hi = c - abig;
        a0lo = dextail - a0hi;
        c = self.splitter * aeytail;
        abig = c - aeytail;
        bhi = c - abig;
        blo = aeytail - bhi;
        _i = dextail * aeytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        dxay[0usize] = a0lo * blo - err3;
        c = self.splitter * dex;
        abig = c - dex;
        a1hi = c - abig;
        a1lo = dex - a1hi;
        _j = dex * aeytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * aey;
        abig = c - aey;
        bhi = c - abig;
        blo = aey - bhi;
        _i = dextail * aey;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxay[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = dex * aey;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxay[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxay[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        dxay[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        dxay[5usize] = around + bround;
        dxay7 = _m + _k;
        bvirt = dxay7 - _m;
        avirt = dxay7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        dxay[6usize] = around + bround;
        dxay[7usize] = dxay7;
        negate = -dey;
        negatetail = -deytail;
        c = self.splitter * aextail;
        abig = c - aextail;
        a0hi = c - abig;
        a0lo = aextail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = aextail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        axdy[0usize] = a0lo * blo - err3;
        c = self.splitter * aex;
        abig = c - aex;
        a1hi = c - abig;
        a1lo = aex - a1hi;
        _j = aex * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = aextail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axdy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = aex * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axdy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axdy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        axdy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        axdy[5usize] = around + bround;
        axdy7 = _m + _k;
        bvirt = axdy7 - _m;
        avirt = axdy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        axdy[6usize] = around + bround;
        axdy[7usize] = axdy7;
        dalen = fast_expansion_sum_zeroelim(
                    8i32,
                    dxay.as_ptr(),
                    8i32,
                    axdy.as_ptr(),
                    da.as_mut_ptr()
                );
        c = self.splitter * aextail;
        abig = c - aextail;
        a0hi = c - abig;
        a0lo = aextail - a0hi;
        c = self.splitter * ceytail;
        abig = c - ceytail;
        bhi = c - abig;
        blo = ceytail - bhi;
        _i = aextail * ceytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        axcy[0usize] = a0lo * blo - err3;
        c = self.splitter * aex;
        abig = c - aex;
        a1hi = c - abig;
        a1lo = aex - a1hi;
        _j = aex * ceytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * cey;
        abig = c - cey;
        bhi = c - abig;
        blo = cey - bhi;
        _i = aextail * cey;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = aex * cey;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        axcy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        axcy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        axcy[5usize] = around + bround;
        axcy7 = _m + _k;
        bvirt = axcy7 - _m;
        avirt = axcy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        axcy[6usize] = around + bround;
        axcy[7usize] = axcy7;
        negate = -aey;
        negatetail = -aeytail;
        c = self.splitter * cextail;
        abig = c - cextail;
        a0hi = c - abig;
        a0lo = cextail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = cextail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        cxay[0usize] = a0lo * blo - err3;
        c = self.splitter * cex;
        abig = c - cex;
        a1hi = c - abig;
        a1lo = cex - a1hi;
        _j = cex * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = cextail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = cex * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        cxay[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        cxay[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        cxay[5usize] = around + bround;
        cxay7 = _m + _k;
        bvirt = cxay7 - _m;
        avirt = cxay7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        cxay[6usize] = around + bround;
        cxay[7usize] = cxay7;
        aclen = fast_expansion_sum_zeroelim(
                    8i32,
                    axcy.as_ptr(),
                    8i32,
                    cxay.as_ptr(),
                    ac.as_mut_ptr()
                );
        c = self.splitter * bextail;
        abig = c - bextail;
        a0hi = c - abig;
        a0lo = bextail - a0hi;
        c = self.splitter * deytail;
        abig = c - deytail;
        bhi = c - abig;
        blo = deytail - bhi;
        _i = bextail * deytail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        bxdy[0usize] = a0lo * blo - err3;
        c = self.splitter * bex;
        abig = c - bex;
        a1hi = c - abig;
        a1lo = bex - a1hi;
        _j = bex * deytail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * dey;
        abig = c - dey;
        bhi = c - abig;
        blo = dey - bhi;
        _i = bextail * dey;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxdy[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = bex * dey;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxdy[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        bxdy[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        bxdy[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        bxdy[5usize] = around + bround;
        bxdy7 = _m + _k;
        bvirt = bxdy7 - _m;
        avirt = bxdy7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        bxdy[6usize] = around + bround;
        bxdy[7usize] = bxdy7;
        negate = -bey;
        negatetail = -beytail;
        c = self.splitter * dextail;
        abig = c - dextail;
        a0hi = c - abig;
        a0lo = dextail - a0hi;
        c = self.splitter * negatetail;
        abig = c - negatetail;
        bhi = c - abig;
        blo = negatetail - bhi;
        _i = dextail * negatetail;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        dxby[0usize] = a0lo * blo - err3;
        c = self.splitter * dex;
        abig = c - dex;
        a1hi = c - abig;
        a1lo = dex - a1hi;
        _j = dex * negatetail;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _k = _i + _0;
        bvirt = _k - _i;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _1 = around + bround;
        _l = _j + _k;
        bvirt = _l - _j;
        _2 = _k - bvirt;
        c = self.splitter * negate;
        abig = c - negate;
        bhi = c - abig;
        blo = negate - bhi;
        _i = dextail * negate;
        err1 = _i - a0hi * bhi;
        err2 = err1 - a0lo * bhi;
        err3 = err2 - a0hi * blo;
        _0 = a0lo * blo - err3;
        _k = _1 + _0;
        bvirt = _k - _1;
        avirt = _k - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxby[1usize] = around + bround;
        _j = _2 + _k;
        bvirt = _j - _2;
        avirt = _j - bvirt;
        bround = _k - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _j;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _j - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _j = dex * negate;
        err1 = _j - a1hi * bhi;
        err2 = err1 - a1lo * bhi;
        err3 = err2 - a1hi * blo;
        _0 = a1lo * blo - err3;
        _n = _i + _0;
        bvirt = _n - _i;
        avirt = _n - bvirt;
        bround = _0 - bvirt;
        around = _i - avirt;
        _0 = around + bround;
        _i = _1 + _0;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxby[2usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _l = _m + _k;
        bvirt = _l - _m;
        avirt = _l - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        _2 = around + bround;
        _k = _j + _n;
        bvirt = _k - _j;
        avirt = _k - bvirt;
        bround = _n - bvirt;
        around = _j - avirt;
        _0 = around + bround;
        _j = _1 + _0;
        bvirt = _j - _1;
        avirt = _j - bvirt;
        bround = _0 - bvirt;
        around = _1 - avirt;
        dxby[3usize] = around + bround;
        _i = _2 + _j;
        bvirt = _i - _2;
        avirt = _i - bvirt;
        bround = _j - bvirt;
        around = _2 - avirt;
        _1 = around + bround;
        _m = _l + _i;
        bvirt = _m - _l;
        avirt = _m - bvirt;
        bround = _i - bvirt;
        around = _l - avirt;
        _2 = around + bround;
        _i = _1 + _k;
        bvirt = _i - _1;
        avirt = _i - bvirt;
        bround = _k - bvirt;
        around = _1 - avirt;
        dxby[4usize] = around + bround;
        _k = _2 + _i;
        bvirt = _k - _2;
        avirt = _k - bvirt;
        bround = _i - bvirt;
        around = _2 - avirt;
        dxby[5usize] = around + bround;
        dxby7 = _m + _k;
        bvirt = dxby7 - _m;
        avirt = dxby7 - bvirt;
        bround = _k - bvirt;
        around = _m - avirt;
        dxby[6usize] = around + bround;
        dxby[7usize] = dxby7;
        bdlen = fast_expansion_sum_zeroelim(
                    8i32,
                    bxdy.as_ptr(),
                    8i32,
                    dxby.as_ptr(),
                    bd.as_mut_ptr()
                );
        temp32alen = self.scale_expansion_zeroelim(
                         cdlen,
                         cd.as_ptr(),
                         -bez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         cdlen,
                         cd.as_ptr(),
                         -beztail,
                         temp32b.as_mut_ptr()
                     );
        temp64alen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64a.as_mut_ptr()
                     );
        temp32alen = self.scale_expansion_zeroelim(
                         bdlen,
                         bd.as_ptr(),
                         cez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         bdlen,
                         bd.as_ptr(),
                         ceztail,
                         temp32b.as_mut_ptr()
                     );
        temp64blen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64b.as_mut_ptr()
                     );
        temp32alen = self.scale_expansion_zeroelim(
                         bclen,
                         bc.as_ptr(),
                         -dez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         bclen,
                         bc.as_ptr(),
                         -deztail,
                         temp32b.as_mut_ptr()
                     );
        temp64clen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64c.as_mut_ptr()
                     );
        temp128len = fast_expansion_sum_zeroelim(
                         temp64alen,
                         temp64a.as_ptr(),
                         temp64blen,
                         temp64b.as_ptr(),
                         temp128.as_mut_ptr()
                     );
        temp192len = fast_expansion_sum_zeroelim(
                         temp64clen,
                         temp64c.as_ptr(),
                         temp128len,
                         temp128.as_ptr(),
                         temp192.as_mut_ptr()
                     );
        xlen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   aex,
                   detx.as_mut_ptr()
               );
        xxlen = self.scale_expansion_zeroelim(
                    xlen,
                    detx.as_ptr(),
                    aex,
                    detxx.as_mut_ptr()
                );
        xtlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    aextail,
                    detxt.as_mut_ptr()
                );
        xxtlen = self.scale_expansion_zeroelim(
                     xtlen,
                     detxt.as_ptr(),
                     aex,
                     detxxt.as_mut_ptr()
                 );
        i = 0i32;
        'loop1: loop {
            if !(i < xxtlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detxxt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        xtxtlen = self.scale_expansion_zeroelim(
                      xtlen,
                      detxt.as_ptr(),
                      aextail,
                      detxtxt.as_mut_ptr()
                  );
        x1len = fast_expansion_sum_zeroelim(
                    xxlen,
                    detxx.as_ptr(),
                    xxtlen,
                    detxxt.as_ptr(),
                    x1.as_mut_ptr()
                );
        x2len = fast_expansion_sum_zeroelim(
                    x1len,
                    x1.as_ptr(),
                    xtxtlen,
                    detxtxt.as_ptr(),
                    x2.as_mut_ptr()
                );
        ylen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   aey,
                   dety.as_mut_ptr()
               );
        yylen = self.scale_expansion_zeroelim(
                    ylen,
                    dety.as_ptr(),
                    aey,
                    detyy.as_mut_ptr()
                );
        ytlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    aeytail,
                    detyt.as_mut_ptr()
                );
        yytlen = self.scale_expansion_zeroelim(
                     ytlen,
                     detyt.as_ptr(),
                     aey,
                     detyyt.as_mut_ptr()
                 );
        i = 0i32;
        'loop3: loop {
            if !(i < yytlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detyyt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ytytlen = self.scale_expansion_zeroelim(
                      ytlen,
                      detyt.as_ptr(),
                      aeytail,
                      detytyt.as_mut_ptr()
                  );
        y1len = fast_expansion_sum_zeroelim(
                    yylen,
                    detyy.as_ptr(),
                    yytlen,
                    detyyt.as_ptr(),
                    y1.as_mut_ptr()
                );
        y2len = fast_expansion_sum_zeroelim(
                    y1len,
                    y1.as_ptr(),
                    ytytlen,
                    detytyt.as_ptr(),
                    y2.as_mut_ptr()
                );
        zlen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   aez,
                   detz.as_mut_ptr()
               );
        zzlen = self.scale_expansion_zeroelim(
                    zlen,
                    detz.as_ptr(),
                    aez,
                    detzz.as_mut_ptr()
                );
        ztlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    aeztail,
                    detzt.as_mut_ptr()
                );
        zztlen = self.scale_expansion_zeroelim(
                     ztlen,
                     detzt.as_ptr(),
                     aez,
                     detzzt.as_mut_ptr()
                 );
        i = 0i32;
        'loop5: loop {
            if !(i < zztlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detzzt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ztztlen = self.scale_expansion_zeroelim(
                      ztlen,
                      detzt.as_ptr(),
                      aeztail,
                      detztzt.as_mut_ptr()
                  );
        z1len = fast_expansion_sum_zeroelim(
                    zzlen,
                    detzz.as_ptr(),
                    zztlen,
                    detzzt.as_ptr(),
                    z1.as_mut_ptr()
                );
        z2len = fast_expansion_sum_zeroelim(
                    z1len,
                    z1.as_ptr(),
                    ztztlen,
                    detztzt.as_ptr(),
                    z2.as_mut_ptr()
                );
        xylen = fast_expansion_sum_zeroelim(
                    x2len,
                    x2.as_ptr(),
                    y2len,
                    y2.as_ptr(),
                    detxy.as_mut_ptr()
                );
        alen = fast_expansion_sum_zeroelim(
                   z2len,
                   z2.as_ptr(),
                   xylen,
                   detxy.as_ptr(),
                   adet.as_mut_ptr()
               );
        temp32alen = self.scale_expansion_zeroelim(
                         dalen,
                         da.as_ptr(),
                         cez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         dalen,
                         da.as_ptr(),
                         ceztail,
                         temp32b.as_mut_ptr()
                     );
        temp64alen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64a.as_mut_ptr()
                     );
        temp32alen = self.scale_expansion_zeroelim(
                         aclen,
                         ac.as_ptr(),
                         dez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         aclen,
                         ac.as_ptr(),
                         deztail,
                         temp32b.as_mut_ptr()
                     );
        temp64blen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64b.as_mut_ptr()
                     );
        temp32alen = self.scale_expansion_zeroelim(
                         cdlen,
                         cd.as_ptr(),
                         aez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         cdlen,
                         cd.as_ptr(),
                         aeztail,
                         temp32b.as_mut_ptr()
                     );
        temp64clen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64c.as_mut_ptr()
                     );
        temp128len = fast_expansion_sum_zeroelim(
                         temp64alen,
                         temp64a.as_ptr(),
                         temp64blen,
                         temp64b.as_ptr(),
                         temp128.as_mut_ptr()
                     );
        temp192len = fast_expansion_sum_zeroelim(
                         temp64clen,
                         temp64c.as_ptr(),
                         temp128len,
                         temp128.as_ptr(),
                         temp192.as_mut_ptr()
                     );
        xlen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   bex,
                   detx.as_mut_ptr()
               );
        xxlen = self.scale_expansion_zeroelim(
                    xlen,
                    detx.as_ptr(),
                    bex,
                    detxx.as_mut_ptr()
                );
        xtlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    bextail,
                    detxt.as_mut_ptr()
                );
        xxtlen = self.scale_expansion_zeroelim(
                     xtlen,
                     detxt.as_ptr(),
                     bex,
                     detxxt.as_mut_ptr()
                 );
        i = 0i32;
        'loop7: loop {
            if !(i < xxtlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detxxt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        xtxtlen = self.scale_expansion_zeroelim(
                      xtlen,
                      detxt.as_ptr(),
                      bextail,
                      detxtxt.as_mut_ptr()
                  );
        x1len = fast_expansion_sum_zeroelim(
                    xxlen,
                    detxx.as_ptr(),
                    xxtlen,
                    detxxt.as_ptr(),
                    x1.as_mut_ptr()
                );
        x2len = fast_expansion_sum_zeroelim(
                    x1len,
                    x1.as_ptr(),
                    xtxtlen,
                    detxtxt.as_ptr(),
                    x2.as_mut_ptr()
                );
        ylen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   bey,
                   dety.as_mut_ptr()
               );
        yylen = self.scale_expansion_zeroelim(
                    ylen,
                    dety.as_ptr(),
                    bey,
                    detyy.as_mut_ptr()
                );
        ytlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    beytail,
                    detyt.as_mut_ptr()
                );
        yytlen = self.scale_expansion_zeroelim(
                     ytlen,
                     detyt.as_ptr(),
                     bey,
                     detyyt.as_mut_ptr()
                 );
        i = 0i32;
        'loop9: loop {
            if !(i < yytlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detyyt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ytytlen = self.scale_expansion_zeroelim(
                      ytlen,
                      detyt.as_ptr(),
                      beytail,
                      detytyt.as_mut_ptr()
                  );
        y1len = fast_expansion_sum_zeroelim(
                    yylen,
                    detyy.as_ptr(),
                    yytlen,
                    detyyt.as_ptr(),
                    y1.as_mut_ptr()
                );
        y2len = fast_expansion_sum_zeroelim(
                    y1len,
                    y1.as_ptr(),
                    ytytlen,
                    detytyt.as_ptr(),
                    y2.as_mut_ptr()
                );
        zlen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   bez,
                   detz.as_mut_ptr()
               );
        zzlen = self.scale_expansion_zeroelim(
                    zlen,
                    detz.as_ptr(),
                    bez,
                    detzz.as_mut_ptr()
                );
        ztlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    beztail,
                    detzt.as_mut_ptr()
                );
        zztlen = self.scale_expansion_zeroelim(
                     ztlen,
                     detzt.as_ptr(),
                     bez,
                     detzzt.as_mut_ptr()
                 );
        i = 0i32;
        'loop11: loop {
            if !(i < zztlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detzzt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ztztlen = self.scale_expansion_zeroelim(
                      ztlen,
                      detzt.as_ptr(),
                      beztail,
                      detztzt.as_mut_ptr()
                  );
        z1len = fast_expansion_sum_zeroelim(
                    zzlen,
                    detzz.as_ptr(),
                    zztlen,
                    detzzt.as_ptr(),
                    z1.as_mut_ptr()
                );
        z2len = fast_expansion_sum_zeroelim(
                    z1len,
                    z1.as_ptr(),
                    ztztlen,
                    detztzt.as_ptr(),
                    z2.as_mut_ptr()
                );
        xylen = fast_expansion_sum_zeroelim(
                    x2len,
                    x2.as_ptr(),
                    y2len,
                    y2.as_ptr(),
                    detxy.as_mut_ptr()
                );
        blen = fast_expansion_sum_zeroelim(
                   z2len,
                   z2.as_ptr(),
                   xylen,
                   detxy.as_ptr(),
                   bdet.as_mut_ptr()
               );
        temp32alen = self.scale_expansion_zeroelim(
                         ablen,
                         ab.as_ptr(),
                         -dez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         ablen,
                         ab.as_ptr(),
                         -deztail,
                         temp32b.as_mut_ptr()
                     );
        temp64alen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64a.as_mut_ptr()
                     );
        temp32alen = self.scale_expansion_zeroelim(
                         bdlen,
                         bd.as_ptr(),
                         -aez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         bdlen,
                         bd.as_ptr(),
                         -aeztail,
                         temp32b.as_mut_ptr()
                     );
        temp64blen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64b.as_mut_ptr()
                     );
        temp32alen = self.scale_expansion_zeroelim(
                         dalen,
                         da.as_ptr(),
                         -bez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         dalen,
                         da.as_ptr(),
                         -beztail,
                         temp32b.as_mut_ptr()
                     );
        temp64clen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64c.as_mut_ptr()
                     );
        temp128len = fast_expansion_sum_zeroelim(
                         temp64alen,
                         temp64a.as_ptr(),
                         temp64blen,
                         temp64b.as_ptr(),
                         temp128.as_mut_ptr()
                     );
        temp192len = fast_expansion_sum_zeroelim(
                         temp64clen,
                         temp64c.as_ptr(),
                         temp128len,
                         temp128.as_ptr(),
                         temp192.as_mut_ptr()
                     );
        xlen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   cex,
                   detx.as_mut_ptr()
               );
        xxlen = self.scale_expansion_zeroelim(
                    xlen,
                    detx.as_ptr(),
                    cex,
                    detxx.as_mut_ptr()
                );
        xtlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    cextail,
                    detxt.as_mut_ptr()
                );
        xxtlen = self.scale_expansion_zeroelim(
                     xtlen,
                     detxt.as_ptr(),
                     cex,
                     detxxt.as_mut_ptr()
                 );
        i = 0i32;
        'loop13: loop {
            if !(i < xxtlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detxxt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        xtxtlen = self.scale_expansion_zeroelim(
                      xtlen,
                      detxt.as_ptr(),
                      cextail,
                      detxtxt.as_mut_ptr()
                  );
        x1len = fast_expansion_sum_zeroelim(
                    xxlen,
                    detxx.as_ptr(),
                    xxtlen,
                    detxxt.as_ptr(),
                    x1.as_mut_ptr()
                );
        x2len = fast_expansion_sum_zeroelim(
                    x1len,
                    x1.as_ptr(),
                    xtxtlen,
                    detxtxt.as_ptr(),
                    x2.as_mut_ptr()
                );
        ylen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   cey,
                   dety.as_mut_ptr()
               );
        yylen = self.scale_expansion_zeroelim(
                    ylen,
                    dety.as_ptr(),
                    cey,
                    detyy.as_mut_ptr()
                );
        ytlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    ceytail,
                    detyt.as_mut_ptr()
                );
        yytlen = self.scale_expansion_zeroelim(
                     ytlen,
                     detyt.as_ptr(),
                     cey,
                     detyyt.as_mut_ptr()
                 );
        i = 0i32;
        'loop15: loop {
            if !(i < yytlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detyyt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ytytlen = self.scale_expansion_zeroelim(
                      ytlen,
                      detyt.as_ptr(),
                      ceytail,
                      detytyt.as_mut_ptr()
                  );
        y1len = fast_expansion_sum_zeroelim(
                    yylen,
                    detyy.as_ptr(),
                    yytlen,
                    detyyt.as_ptr(),
                    y1.as_mut_ptr()
                );
        y2len = fast_expansion_sum_zeroelim(
                    y1len,
                    y1.as_ptr(),
                    ytytlen,
                    detytyt.as_ptr(),
                    y2.as_mut_ptr()
                );
        zlen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   cez,
                   detz.as_mut_ptr()
               );
        zzlen = self.scale_expansion_zeroelim(
                    zlen,
                    detz.as_ptr(),
                    cez,
                    detzz.as_mut_ptr()
                );
        ztlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    ceztail,
                    detzt.as_mut_ptr()
                );
        zztlen = self.scale_expansion_zeroelim(
                     ztlen,
                     detzt.as_ptr(),
                     cez,
                     detzzt.as_mut_ptr()
                 );
        i = 0i32;
        'loop17: loop {
            if !(i < zztlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detzzt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ztztlen = self.scale_expansion_zeroelim(
                      ztlen,
                      detzt.as_ptr(),
                      ceztail,
                      detztzt.as_mut_ptr()
                  );
        z1len = fast_expansion_sum_zeroelim(
                    zzlen,
                    detzz.as_ptr(),
                    zztlen,
                    detzzt.as_ptr(),
                    z1.as_mut_ptr()
                );
        z2len = fast_expansion_sum_zeroelim(
                    z1len,
                    z1.as_ptr(),
                    ztztlen,
                    detztzt.as_ptr(),
                    z2.as_mut_ptr()
                );
        xylen = fast_expansion_sum_zeroelim(
                    x2len,
                    x2.as_ptr(),
                    y2len,
                    y2.as_ptr(),
                    detxy.as_mut_ptr()
                );
        clen = fast_expansion_sum_zeroelim(
                   z2len,
                   z2.as_ptr(),
                   xylen,
                   detxy.as_ptr(),
                   cdet.as_mut_ptr()
               );
        temp32alen = self.scale_expansion_zeroelim(
                         bclen,
                         bc.as_ptr(),
                         aez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         bclen,
                         bc.as_ptr(),
                         aeztail,
                         temp32b.as_mut_ptr()
                     );
        temp64alen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64a.as_mut_ptr()
                     );
        temp32alen = self.scale_expansion_zeroelim(
                         aclen,
                         ac.as_ptr(),
                         -bez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         aclen,
                         ac.as_ptr(),
                         -beztail,
                         temp32b.as_mut_ptr()
                     );
        temp64blen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64b.as_mut_ptr()
                     );
        temp32alen = self.scale_expansion_zeroelim(
                         ablen,
                         ab.as_ptr(),
                         cez,
                         temp32a.as_mut_ptr()
                     );
        temp32blen = self.scale_expansion_zeroelim(
                         ablen,
                         ab.as_ptr(),
                         ceztail,
                         temp32b.as_mut_ptr()
                     );
        temp64clen = fast_expansion_sum_zeroelim(
                         temp32alen,
                         temp32a.as_ptr(),
                         temp32blen,
                         temp32b.as_ptr(),
                         temp64c.as_mut_ptr()
                     );
        temp128len = fast_expansion_sum_zeroelim(
                         temp64alen,
                         temp64a.as_ptr(),
                         temp64blen,
                         temp64b.as_ptr(),
                         temp128.as_mut_ptr()
                     );
        temp192len = fast_expansion_sum_zeroelim(
                         temp64clen,
                         temp64c.as_ptr(),
                         temp128len,
                         temp128.as_ptr(),
                         temp192.as_mut_ptr()
                     );
        xlen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   dex,
                   detx.as_mut_ptr()
               );
        xxlen = self.scale_expansion_zeroelim(
                    xlen,
                    detx.as_ptr(),
                    dex,
                    detxx.as_mut_ptr()
                );
        xtlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    dextail,
                    detxt.as_mut_ptr()
                );
        xxtlen = self.scale_expansion_zeroelim(
                     xtlen,
                     detxt.as_ptr(),
                     dex,
                     detxxt.as_mut_ptr()
                 );
        i = 0i32;
        'loop19: loop {
            if !(i < xxtlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detxxt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        xtxtlen = self.scale_expansion_zeroelim(
                      xtlen,
                      detxt.as_ptr(),
                      dextail,
                      detxtxt.as_mut_ptr()
                  );
        x1len = fast_expansion_sum_zeroelim(
                    xxlen,
                    detxx.as_ptr(),
                    xxtlen,
                    detxxt.as_ptr(),
                    x1.as_mut_ptr()
                );
        x2len = fast_expansion_sum_zeroelim(
                    x1len,
                    x1.as_ptr(),
                    xtxtlen,
                    detxtxt.as_ptr(),
                    x2.as_mut_ptr()
                );
        ylen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   dey,
                   dety.as_mut_ptr()
               );
        yylen = self.scale_expansion_zeroelim(
                    ylen,
                    dety.as_ptr(),
                    dey,
                    detyy.as_mut_ptr()
                );
        ytlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    deytail,
                    detyt.as_mut_ptr()
                );
        yytlen = self.scale_expansion_zeroelim(
                     ytlen,
                     detyt.as_ptr(),
                     dey,
                     detyyt.as_mut_ptr()
                 );
        i = 0i32;
        'loop21: loop {
            if !(i < yytlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detyyt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ytytlen = self.scale_expansion_zeroelim(
                      ytlen,
                      detyt.as_ptr(),
                      deytail,
                      detytyt.as_mut_ptr()
                  );
        y1len = fast_expansion_sum_zeroelim(
                    yylen,
                    detyy.as_ptr(),
                    yytlen,
                    detyyt.as_ptr(),
                    y1.as_mut_ptr()
                );
        y2len = fast_expansion_sum_zeroelim(
                    y1len,
                    y1.as_ptr(),
                    ytytlen,
                    detytyt.as_ptr(),
                    y2.as_mut_ptr()
                );
        zlen = self.scale_expansion_zeroelim(
                   temp192len,
                   temp192.as_ptr(),
                   dez,
                   detz.as_mut_ptr()
               );
        zzlen = self.scale_expansion_zeroelim(
                    zlen,
                    detz.as_ptr(),
                    dez,
                    detzz.as_mut_ptr()
                );
        ztlen = self.scale_expansion_zeroelim(
                    temp192len,
                    temp192.as_ptr(),
                    deztail,
                    detzt.as_mut_ptr()
                );
        zztlen = self.scale_expansion_zeroelim(
                     ztlen,
                     detzt.as_ptr(),
                     dez,
                     detzzt.as_mut_ptr()
                 );
        i = 0i32;
        'loop23: loop {
            if !(i < zztlen) {
                break;
            }
            let _rhs = 2.0f64;
            let _lhs = &mut detzzt[i as (usize)];
            *_lhs = *_lhs * _rhs;
            i = i + 1;
        }
        ztztlen = self.scale_expansion_zeroelim(
                      ztlen,
                      detzt.as_ptr(),
                      deztail,
                      detztzt.as_mut_ptr()
                  );
        z1len = fast_expansion_sum_zeroelim(
                    zzlen,
                    detzz.as_ptr(),
                    zztlen,
                    detzzt.as_ptr(),
                    z1.as_mut_ptr()
                );
        z2len = fast_expansion_sum_zeroelim(
                    z1len,
                    z1.as_ptr(),
                    ztztlen,
                    detztzt.as_ptr(),
                    z2.as_mut_ptr()
                );
        xylen = fast_expansion_sum_zeroelim(
                    x2len,
                    x2.as_ptr(),
                    y2len,
                    y2.as_ptr(),
                    detxy.as_mut_ptr()
                );
        dlen = fast_expansion_sum_zeroelim(
                   z2len,
                   z2.as_ptr(),
                   xylen,
                   detxy.as_ptr(),
                   ddet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        cdlen = fast_expansion_sum_zeroelim(
                    clen,
                    cdet.as_ptr(),
                    dlen,
                    ddet.as_ptr(),
                    cddet.as_mut_ptr()
                );
        deterlen = fast_expansion_sum_zeroelim(
                       ablen,
                       abdet.as_ptr(),
                       cdlen,
                       cddet.as_ptr(),
                       deter.as_mut_ptr()
                   );
        deter[(deterlen - 1i32) as (usize)]
    }

    
    pub unsafe fn insphereadapt(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64,
        mut pe : *const f64,
        mut permanent : f64
    ) -> f64 {
        let mut aex : f64;
        let mut bex : f64;
        let mut cex : f64;
        let mut dex : f64;
        let mut aey : f64;
        let mut bey : f64;
        let mut cey : f64;
        let mut dey : f64;
        let mut aez : f64;
        let mut bez : f64;
        let mut cez : f64;
        let mut dez : f64;
        let mut det : f64;
        let mut errbound : f64;
        let mut aexbey1 : f64;
        let mut bexaey1 : f64;
        let mut bexcey1 : f64;
        let mut cexbey1 : f64;
        let mut cexdey1 : f64;
        let mut dexcey1 : f64;
        let mut dexaey1 : f64;
        let mut aexdey1 : f64;
        let mut aexcey1 : f64;
        let mut cexaey1 : f64;
        let mut bexdey1 : f64;
        let mut dexbey1 : f64;
        let mut aexbey0 : f64;
        let mut bexaey0 : f64;
        let mut bexcey0 : f64;
        let mut cexbey0 : f64;
        let mut cexdey0 : f64;
        let mut dexcey0 : f64;
        let mut dexaey0 : f64;
        let mut aexdey0 : f64;
        let mut aexcey0 : f64;
        let mut cexaey0 : f64;
        let mut bexdey0 : f64;
        let mut dexbey0 : f64;
        let mut ab : [f64; 4];
        let mut bc : [f64; 4];
        let mut cd : [f64; 4];
        let mut da : [f64; 4];
        let mut ac : [f64; 4];
        let mut bd : [f64; 4];
        let mut ab3 : f64;
        let mut bc3 : f64;
        let mut cd3 : f64;
        let mut da3 : f64;
        let mut ac3 : f64;
        let mut bd3 : f64;
        let mut abeps : f64;
        let mut bceps : f64;
        let mut cdeps : f64;
        let mut daeps : f64;
        let mut aceps : f64;
        let mut bdeps : f64;
        let mut temp8a : [f64; 8];
        let mut temp8b : [f64; 8];
        let mut temp8c : [f64; 8];
        let mut temp16 : [f64; 16];
        let mut temp24 : [f64; 24];
        let mut temp48 : [f64; 48];
        let mut temp8alen : i32;
        let mut temp8blen : i32;
        let mut temp8clen : i32;
        let mut temp16len : i32;
        let mut temp24len : i32;
        let mut temp48len : i32;
        let mut xdet : [f64; 96];
        let mut ydet : [f64; 96];
        let mut zdet : [f64; 96];
        let mut xydet : [f64; 192];
        let mut xlen : i32;
        let mut ylen : i32;
        let mut zlen : i32;
        let mut xylen : i32;
        let mut adet : [f64; 288];
        let mut bdet : [f64; 288];
        let mut cdet : [f64; 288];
        let mut ddet : [f64; 288];
        let mut alen : i32;
        let mut blen : i32;
        let mut clen : i32;
        let mut dlen : i32;
        let mut abdet : [f64; 576];
        let mut cddet : [f64; 576];
        let mut ablen : i32;
        let mut cdlen : i32;
        let mut fin1 : [f64; 1152];
        let mut finlength : i32;
        let mut aextail : f64;
        let mut bextail : f64;
        let mut cextail : f64;
        let mut dextail : f64;
        let mut aeytail : f64;
        let mut beytail : f64;
        let mut ceytail : f64;
        let mut deytail : f64;
        let mut aeztail : f64;
        let mut beztail : f64;
        let mut ceztail : f64;
        let mut deztail : f64;
        let mut bvirt : f64;
        let mut avirt : f64;
        let mut bround : f64;
        let mut around : f64;
        let mut c : f64;
        let mut abig : f64;
        let mut ahi : f64;
        let mut alo : f64;
        let mut bhi : f64;
        let mut blo : f64;
        let mut err1 : f64;
        let mut err2 : f64;
        let mut err3 : f64;
        let mut _i : f64;
        let mut _j : f64;
        let mut _0 : f64;

        aex = uninitialized();
        bex = uninitialized();
        cex = uninitialized();
        dex = uninitialized();
        aey = uninitialized();
        bey = uninitialized();
        cey = uninitialized();
        dey = uninitialized();
        aez = uninitialized();
        bez = uninitialized();
        cez = uninitialized();
        dez = uninitialized();
        det = uninitialized();
        errbound = uninitialized();
        aexbey1 = uninitialized();
        bexaey1 = uninitialized();
        bexcey1 = uninitialized();
        cexbey1 = uninitialized();
        cexdey1 = uninitialized();
        dexcey1 = uninitialized();
        dexaey1 = uninitialized();
        aexdey1 = uninitialized();
        aexcey1 = uninitialized();
        cexaey1 = uninitialized();
        bexdey1 = uninitialized();
        dexbey1 = uninitialized();
        aexbey0 = uninitialized();
        bexaey0 = uninitialized();
        bexcey0 = uninitialized();
        cexbey0 = uninitialized();
        cexdey0 = uninitialized();
        dexcey0 = uninitialized();
        dexaey0 = uninitialized();
        aexdey0 = uninitialized();
        aexcey0 = uninitialized();
        cexaey0 = uninitialized();
        bexdey0 = uninitialized();
        dexbey0 = uninitialized();
        ab = uninitialized();
        bc = uninitialized();
        cd = uninitialized();
        da = uninitialized();
        ac = uninitialized();
        bd = uninitialized();
        ab3 = uninitialized();
        bc3 = uninitialized();
        cd3 = uninitialized();
        da3 = uninitialized();
        ac3 = uninitialized();
        bd3 = uninitialized();
        abeps = uninitialized();
        bceps = uninitialized();
        cdeps = uninitialized();
        daeps = uninitialized();
        aceps = uninitialized();
        bdeps = uninitialized();
        temp8a = uninitialized();
        temp8b = uninitialized();
        temp8c = uninitialized();
        temp16 = uninitialized();
        temp24 = uninitialized();
        temp48 = uninitialized();
        temp8alen = uninitialized();
        temp8blen = uninitialized();
        temp8clen = uninitialized();
        temp16len = uninitialized();
        temp24len = uninitialized();
        temp48len = uninitialized();
        xdet = uninitialized();
        ydet = uninitialized();
        zdet = uninitialized();
        xydet = uninitialized();
        xlen = uninitialized();
        ylen = uninitialized();
        zlen = uninitialized();
        xylen = uninitialized();
        adet = uninitialized();
        bdet = uninitialized();
        cdet = uninitialized();
        ddet = uninitialized();
        alen = uninitialized();
        blen = uninitialized();
        clen = uninitialized();
        dlen = uninitialized();
        abdet = uninitialized();
        cddet = uninitialized();
        ablen = uninitialized();
        cdlen = uninitialized();
        fin1 = uninitialized();
        finlength = uninitialized();
        aextail = uninitialized();
        bextail = uninitialized();
        cextail = uninitialized();
        dextail = uninitialized();
        aeytail = uninitialized();
        beytail = uninitialized();
        ceytail = uninitialized();
        deytail = uninitialized();
        aeztail = uninitialized();
        beztail = uninitialized();
        ceztail = uninitialized();
        deztail = uninitialized();
        bvirt = uninitialized();
        avirt = uninitialized();
        bround = uninitialized();
        around = uninitialized();
        c = uninitialized();
        abig = uninitialized();
        ahi = uninitialized();
        alo = uninitialized();
        bhi = uninitialized();
        blo = uninitialized();
        err1 = uninitialized();
        err2 = uninitialized();
        err3 = uninitialized();
        _i = uninitialized();
        _j = uninitialized();
        _0 = uninitialized();
        aex = *pa.offset(0isize) - *pe.offset(0isize);
        bex = *pb.offset(0isize) - *pe.offset(0isize);
        cex = *pc.offset(0isize) - *pe.offset(0isize);
        dex = *pd.offset(0isize) - *pe.offset(0isize);
        aey = *pa.offset(1isize) - *pe.offset(1isize);
        bey = *pb.offset(1isize) - *pe.offset(1isize);
        cey = *pc.offset(1isize) - *pe.offset(1isize);
        dey = *pd.offset(1isize) - *pe.offset(1isize);
        aez = *pa.offset(2isize) - *pe.offset(2isize);
        bez = *pb.offset(2isize) - *pe.offset(2isize);
        cez = *pc.offset(2isize) - *pe.offset(2isize);
        dez = *pd.offset(2isize) - *pe.offset(2isize);
        aexbey1 = aex * bey;
        c = self.splitter * aex;
        abig = c - aex;
        ahi = c - abig;
        alo = aex - ahi;
        c = self.splitter * bey;
        abig = c - bey;
        bhi = c - abig;
        blo = bey - bhi;
        err1 = aexbey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        aexbey0 = alo * blo - err3;
        bexaey1 = bex * aey;
        c = self.splitter * bex;
        abig = c - bex;
        ahi = c - abig;
        alo = bex - ahi;
        c = self.splitter * aey;
        abig = c - aey;
        bhi = c - abig;
        blo = aey - bhi;
        err1 = bexaey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bexaey0 = alo * blo - err3;
        _i = aexbey0 - bexaey0;
        bvirt = aexbey0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bexaey0;
        around = aexbey0 - avirt;
        ab[0usize] = around + bround;
        _j = aexbey1 + _i;
        bvirt = _j - aexbey1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = aexbey1 - avirt;
        _0 = around + bround;
        _i = _0 - bexaey1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - bexaey1;
        around = _0 - avirt;
        ab[1usize] = around + bround;
        ab3 = _j + _i;
        bvirt = ab3 - _j;
        avirt = ab3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ab[2usize] = around + bround;
        ab[3usize] = ab3;
        bexcey1 = bex * cey;
        c = self.splitter * bex;
        abig = c - bex;
        ahi = c - abig;
        alo = bex - ahi;
        c = self.splitter * cey;
        abig = c - cey;
        bhi = c - abig;
        blo = cey - bhi;
        err1 = bexcey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bexcey0 = alo * blo - err3;
        cexbey1 = cex * bey;
        c = self.splitter * cex;
        abig = c - cex;
        ahi = c - abig;
        alo = cex - ahi;
        c = self.splitter * bey;
        abig = c - bey;
        bhi = c - abig;
        blo = bey - bhi;
        err1 = cexbey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cexbey0 = alo * blo - err3;
        _i = bexcey0 - cexbey0;
        bvirt = bexcey0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cexbey0;
        around = bexcey0 - avirt;
        bc[0usize] = around + bround;
        _j = bexcey1 + _i;
        bvirt = _j - bexcey1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bexcey1 - avirt;
        _0 = around + bround;
        _i = _0 - cexbey1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cexbey1;
        around = _0 - avirt;
        bc[1usize] = around + bround;
        bc3 = _j + _i;
        bvirt = bc3 - _j;
        avirt = bc3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bc[2usize] = around + bround;
        bc[3usize] = bc3;
        cexdey1 = cex * dey;
        c = self.splitter * cex;
        abig = c - cex;
        ahi = c - abig;
        alo = cex - ahi;
        c = self.splitter * dey;
        abig = c - dey;
        bhi = c - abig;
        blo = dey - bhi;
        err1 = cexdey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cexdey0 = alo * blo - err3;
        dexcey1 = dex * cey;
        c = self.splitter * dex;
        abig = c - dex;
        ahi = c - abig;
        alo = dex - ahi;
        c = self.splitter * cey;
        abig = c - cey;
        bhi = c - abig;
        blo = cey - bhi;
        err1 = dexcey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dexcey0 = alo * blo - err3;
        _i = cexdey0 - dexcey0;
        bvirt = cexdey0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dexcey0;
        around = cexdey0 - avirt;
        cd[0usize] = around + bround;
        _j = cexdey1 + _i;
        bvirt = _j - cexdey1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = cexdey1 - avirt;
        _0 = around + bround;
        _i = _0 - dexcey1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dexcey1;
        around = _0 - avirt;
        cd[1usize] = around + bround;
        cd3 = _j + _i;
        bvirt = cd3 - _j;
        avirt = cd3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        cd[2usize] = around + bround;
        cd[3usize] = cd3;
        dexaey1 = dex * aey;
        c = self.splitter * dex;
        abig = c - dex;
        ahi = c - abig;
        alo = dex - ahi;
        c = self.splitter * aey;
        abig = c - aey;
        bhi = c - abig;
        blo = aey - bhi;
        err1 = dexaey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dexaey0 = alo * blo - err3;
        aexdey1 = aex * dey;
        c = self.splitter * aex;
        abig = c - aex;
        ahi = c - abig;
        alo = aex - ahi;
        c = self.splitter * dey;
        abig = c - dey;
        bhi = c - abig;
        blo = dey - bhi;
        err1 = aexdey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        aexdey0 = alo * blo - err3;
        _i = dexaey0 - aexdey0;
        bvirt = dexaey0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - aexdey0;
        around = dexaey0 - avirt;
        da[0usize] = around + bround;
        _j = dexaey1 + _i;
        bvirt = _j - dexaey1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = dexaey1 - avirt;
        _0 = around + bround;
        _i = _0 - aexdey1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - aexdey1;
        around = _0 - avirt;
        da[1usize] = around + bround;
        da3 = _j + _i;
        bvirt = da3 - _j;
        avirt = da3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        da[2usize] = around + bround;
        da[3usize] = da3;
        aexcey1 = aex * cey;
        c = self.splitter * aex;
        abig = c - aex;
        ahi = c - abig;
        alo = aex - ahi;
        c = self.splitter * cey;
        abig = c - cey;
        bhi = c - abig;
        blo = cey - bhi;
        err1 = aexcey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        aexcey0 = alo * blo - err3;
        cexaey1 = cex * aey;
        c = self.splitter * cex;
        abig = c - cex;
        ahi = c - abig;
        alo = cex - ahi;
        c = self.splitter * aey;
        abig = c - aey;
        bhi = c - abig;
        blo = aey - bhi;
        err1 = cexaey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        cexaey0 = alo * blo - err3;
        _i = aexcey0 - cexaey0;
        bvirt = aexcey0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cexaey0;
        around = aexcey0 - avirt;
        ac[0usize] = around + bround;
        _j = aexcey1 + _i;
        bvirt = _j - aexcey1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = aexcey1 - avirt;
        _0 = around + bround;
        _i = _0 - cexaey1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - cexaey1;
        around = _0 - avirt;
        ac[1usize] = around + bround;
        ac3 = _j + _i;
        bvirt = ac3 - _j;
        avirt = ac3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        ac[2usize] = around + bround;
        ac[3usize] = ac3;
        bexdey1 = bex * dey;
        c = self.splitter * bex;
        abig = c - bex;
        ahi = c - abig;
        alo = bex - ahi;
        c = self.splitter * dey;
        abig = c - dey;
        bhi = c - abig;
        blo = dey - bhi;
        err1 = bexdey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        bexdey0 = alo * blo - err3;
        dexbey1 = dex * bey;
        c = self.splitter * dex;
        abig = c - dex;
        ahi = c - abig;
        alo = dex - ahi;
        c = self.splitter * bey;
        abig = c - bey;
        bhi = c - abig;
        blo = bey - bhi;
        err1 = dexbey1 - ahi * bhi;
        err2 = err1 - alo * bhi;
        err3 = err2 - ahi * blo;
        dexbey0 = alo * blo - err3;
        _i = bexdey0 - dexbey0;
        bvirt = bexdey0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dexbey0;
        around = bexdey0 - avirt;
        bd[0usize] = around + bround;
        _j = bexdey1 + _i;
        bvirt = _j - bexdey1;
        avirt = _j - bvirt;
        bround = _i - bvirt;
        around = bexdey1 - avirt;
        _0 = around + bround;
        _i = _0 - dexbey1;
        bvirt = _0 - _i;
        avirt = _i + bvirt;
        bround = bvirt - dexbey1;
        around = _0 - avirt;
        bd[1usize] = around + bround;
        bd3 = _j + _i;
        bvirt = bd3 - _j;
        avirt = bd3 - bvirt;
        bround = _i - bvirt;
        around = _j - avirt;
        bd[2usize] = around + bround;
        bd[3usize] = bd3;
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        cd.as_ptr(),
                        bez,
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        bd.as_ptr(),
                        -cez,
                        temp8b.as_mut_ptr()
                    );
        temp8clen = self.scale_expansion_zeroelim(
                        4i32,
                        bc.as_ptr(),
                        dez,
                        temp8c.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp24len = fast_expansion_sum_zeroelim(
                        temp8clen,
                        temp8c.as_ptr(),
                        temp16len,
                        temp16.as_ptr(),
                        temp24.as_mut_ptr()
                    );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        aex,
                        temp48.as_mut_ptr()
                    );
        xlen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   -aex,
                   xdet.as_mut_ptr()
               );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        aey,
                        temp48.as_mut_ptr()
                    );
        ylen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   -aey,
                   ydet.as_mut_ptr()
               );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        aez,
                        temp48.as_mut_ptr()
                    );
        zlen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   -aez,
                   zdet.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    xdet.as_ptr(),
                    ylen,
                    ydet.as_ptr(),
                    xydet.as_mut_ptr()
                );
        alen = fast_expansion_sum_zeroelim(
                   xylen,
                   xydet.as_ptr(),
                   zlen,
                   zdet.as_ptr(),
                   adet.as_mut_ptr()
               );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        da.as_ptr(),
                        cez,
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        ac.as_ptr(),
                        dez,
                        temp8b.as_mut_ptr()
                    );
        temp8clen = self.scale_expansion_zeroelim(
                        4i32,
                        cd.as_ptr(),
                        aez,
                        temp8c.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp24len = fast_expansion_sum_zeroelim(
                        temp8clen,
                        temp8c.as_ptr(),
                        temp16len,
                        temp16.as_ptr(),
                        temp24.as_mut_ptr()
                    );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        bex,
                        temp48.as_mut_ptr()
                    );
        xlen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   bex,
                   xdet.as_mut_ptr()
               );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        bey,
                        temp48.as_mut_ptr()
                    );
        ylen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   bey,
                   ydet.as_mut_ptr()
               );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        bez,
                        temp48.as_mut_ptr()
                    );
        zlen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   bez,
                   zdet.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    xdet.as_ptr(),
                    ylen,
                    ydet.as_ptr(),
                    xydet.as_mut_ptr()
                );
        blen = fast_expansion_sum_zeroelim(
                   xylen,
                   xydet.as_ptr(),
                   zlen,
                   zdet.as_ptr(),
                   bdet.as_mut_ptr()
               );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        ab.as_ptr(),
                        dez,
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        bd.as_ptr(),
                        aez,
                        temp8b.as_mut_ptr()
                    );
        temp8clen = self.scale_expansion_zeroelim(
                        4i32,
                        da.as_ptr(),
                        bez,
                        temp8c.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp24len = fast_expansion_sum_zeroelim(
                        temp8clen,
                        temp8c.as_ptr(),
                        temp16len,
                        temp16.as_ptr(),
                        temp24.as_mut_ptr()
                    );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        cex,
                        temp48.as_mut_ptr()
                    );
        xlen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   -cex,
                   xdet.as_mut_ptr()
               );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        cey,
                        temp48.as_mut_ptr()
                    );
        ylen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   -cey,
                   ydet.as_mut_ptr()
               );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        cez,
                        temp48.as_mut_ptr()
                    );
        zlen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   -cez,
                   zdet.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    xdet.as_ptr(),
                    ylen,
                    ydet.as_ptr(),
                    xydet.as_mut_ptr()
                );
        clen = fast_expansion_sum_zeroelim(
                   xylen,
                   xydet.as_ptr(),
                   zlen,
                   zdet.as_ptr(),
                   cdet.as_mut_ptr()
               );
        temp8alen = self.scale_expansion_zeroelim(
                        4i32,
                        bc.as_ptr(),
                        aez,
                        temp8a.as_mut_ptr()
                    );
        temp8blen = self.scale_expansion_zeroelim(
                        4i32,
                        ac.as_ptr(),
                        -bez,
                        temp8b.as_mut_ptr()
                    );
        temp8clen = self.scale_expansion_zeroelim(
                        4i32,
                        ab.as_ptr(),
                        cez,
                        temp8c.as_mut_ptr()
                    );
        temp16len = fast_expansion_sum_zeroelim(
                        temp8alen,
                        temp8a.as_ptr(),
                        temp8blen,
                        temp8b.as_ptr(),
                        temp16.as_mut_ptr()
                    );
        temp24len = fast_expansion_sum_zeroelim(
                        temp8clen,
                        temp8c.as_ptr(),
                        temp16len,
                        temp16.as_ptr(),
                        temp24.as_mut_ptr()
                    );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        dex,
                        temp48.as_mut_ptr()
                    );
        xlen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   dex,
                   xdet.as_mut_ptr()
               );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        dey,
                        temp48.as_mut_ptr()
                    );
        ylen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   dey,
                   ydet.as_mut_ptr()
               );
        temp48len = self.scale_expansion_zeroelim(
                        temp24len,
                        temp24.as_ptr(),
                        dez,
                        temp48.as_mut_ptr()
                    );
        zlen = self.scale_expansion_zeroelim(
                   temp48len,
                   temp48.as_ptr(),
                   dez,
                   zdet.as_mut_ptr()
               );
        xylen = fast_expansion_sum_zeroelim(
                    xlen,
                    xdet.as_ptr(),
                    ylen,
                    ydet.as_ptr(),
                    xydet.as_mut_ptr()
                );
        dlen = fast_expansion_sum_zeroelim(
                   xylen,
                   xydet.as_ptr(),
                   zlen,
                   zdet.as_ptr(),
                   ddet.as_mut_ptr()
               );
        ablen = fast_expansion_sum_zeroelim(
                    alen,
                    adet.as_ptr(),
                    blen,
                    bdet.as_ptr(),
                    abdet.as_mut_ptr()
                );
        cdlen = fast_expansion_sum_zeroelim(
                    clen,
                    cdet.as_ptr(),
                    dlen,
                    ddet.as_ptr(),
                    cddet.as_mut_ptr()
                );
        finlength = fast_expansion_sum_zeroelim(
                        ablen,
                        abdet.as_ptr(),
                        cdlen,
                        cddet.as_ptr(),
                        fin1.as_mut_ptr()
                    );
        det = estimate(finlength,fin1.as_ptr());
        errbound = self.isperrboundB * permanent;
        if det >= errbound || -det >= errbound {
            det
        } else {
            bvirt = *pa.offset(0isize) - aex;
            avirt = aex + bvirt;
            bround = bvirt - *pe.offset(0isize);
            around = *pa.offset(0isize) - avirt;
            aextail = around + bround;
            bvirt = *pa.offset(1isize) - aey;
            avirt = aey + bvirt;
            bround = bvirt - *pe.offset(1isize);
            around = *pa.offset(1isize) - avirt;
            aeytail = around + bround;
            bvirt = *pa.offset(2isize) - aez;
            avirt = aez + bvirt;
            bround = bvirt - *pe.offset(2isize);
            around = *pa.offset(2isize) - avirt;
            aeztail = around + bround;
            bvirt = *pb.offset(0isize) - bex;
            avirt = bex + bvirt;
            bround = bvirt - *pe.offset(0isize);
            around = *pb.offset(0isize) - avirt;
            bextail = around + bround;
            bvirt = *pb.offset(1isize) - bey;
            avirt = bey + bvirt;
            bround = bvirt - *pe.offset(1isize);
            around = *pb.offset(1isize) - avirt;
            beytail = around + bround;
            bvirt = *pb.offset(2isize) - bez;
            avirt = bez + bvirt;
            bround = bvirt - *pe.offset(2isize);
            around = *pb.offset(2isize) - avirt;
            beztail = around + bround;
            bvirt = *pc.offset(0isize) - cex;
            avirt = cex + bvirt;
            bround = bvirt - *pe.offset(0isize);
            around = *pc.offset(0isize) - avirt;
            cextail = around + bround;
            bvirt = *pc.offset(1isize) - cey;
            avirt = cey + bvirt;
            bround = bvirt - *pe.offset(1isize);
            around = *pc.offset(1isize) - avirt;
            ceytail = around + bround;
            bvirt = *pc.offset(2isize) - cez;
            avirt = cez + bvirt;
            bround = bvirt - *pe.offset(2isize);
            around = *pc.offset(2isize) - avirt;
            ceztail = around + bround;
            bvirt = *pd.offset(0isize) - dex;
            avirt = dex + bvirt;
            bround = bvirt - *pe.offset(0isize);
            around = *pd.offset(0isize) - avirt;
            dextail = around + bround;
            bvirt = *pd.offset(1isize) - dey;
            avirt = dey + bvirt;
            bround = bvirt - *pe.offset(1isize);
            around = *pd.offset(1isize) - avirt;
            deytail = around + bround;
            bvirt = *pd.offset(2isize) - dez;
            avirt = dez + bvirt;
            bround = bvirt - *pe.offset(2isize);
            around = *pd.offset(2isize) - avirt;
            deztail = around + bround;
            (if aextail == 0.0f64 && (aeytail == 0.0f64) && (aeztail == 0.0f64) && (bextail == 0.0f64) && (beytail == 0.0f64) && (beztail == 0.0f64) && (cextail == 0.0f64) && (ceytail == 0.0f64) && (ceztail == 0.0f64) && (dextail == 0.0f64) && (deytail == 0.0f64) && (deztail == 0.0f64) {
                 det
             } else {
                 errbound = self.isperrboundC * permanent + self.resulterrbound * Absolute(
                                                                            det
                                                                        );
                 abeps = aex * beytail + bey * aextail - (aey * bextail + bex * aeytail);
                 bceps = bex * ceytail + cey * bextail - (bey * cextail + cex * beytail);
                 cdeps = cex * deytail + dey * cextail - (cey * dextail + dex * ceytail);
                 daeps = dex * aeytail + aey * dextail - (dey * aextail + aex * deytail);
                 aceps = aex * ceytail + cey * aextail - (aey * cextail + cex * aeytail);
                 bdeps = bex * deytail + dey * bextail - (bey * dextail + dex * beytail);
                 det = det + ((bex * bex + bey * bey + bez * bez) * (cez * daeps + dez * aceps + aez * cdeps + (ceztail * da3 + deztail * ac3 + aeztail * cd3)) + (dex * dex + dey * dey + dez * dez) * (aez * bceps - bez * aceps + cez * abeps + (aeztail * bc3 - beztail * ac3 + ceztail * ab3)) - ((aex * aex + aey * aey + aez * aez) * (bez * cdeps - cez * bdeps + dez * bceps + (beztail * cd3 - ceztail * bd3 + deztail * bc3)) + (cex * cex + cey * cey + cez * cez) * (dez * abeps + aez * bdeps + bez * daeps + (deztail * ab3 + aeztail * bd3 + beztail * da3))) + 2.0f64 * ((bex * bextail + bey * beytail + bez * beztail) * (cez * da3 + dez * ac3 + aez * cd3) + (dex * dextail + dey * deytail + dez * deztail) * (aez * bc3 - bez * ac3 + cez * ab3) - ((aex * aextail + aey * aeytail + aez * aeztail) * (bez * cd3 - cez * bd3 + dez * bc3) + (cex * cextail + cey * ceytail + cez * ceztail) * (dez * ab3 + aez * bd3 + bez * da3))));
                 (if det >= errbound || -det >= errbound {
                      det
                  } else {
                      self.insphereexact(pa,pb,pc,pd,pe)
                  })
             })
        }
    }

    
    pub unsafe fn insphere(&self,
        mut pa : *const f64,
        mut pb : *const f64,
        mut pc : *const f64,
        mut pd : *const f64,
        mut pe : *const f64
    ) -> f64 {
        let mut aex : f64;
        let mut bex : f64;
        let mut cex : f64;
        let mut dex : f64;
        let mut aey : f64;
        let mut bey : f64;
        let mut cey : f64;
        let mut dey : f64;
        let mut aez : f64;
        let mut bez : f64;
        let mut cez : f64;
        let mut dez : f64;
        let mut aexbey : f64;
        let mut bexaey : f64;
        let mut bexcey : f64;
        let mut cexbey : f64;
        let mut cexdey : f64;
        let mut dexcey : f64;
        let mut dexaey : f64;
        let mut aexdey : f64;
        let mut aexcey : f64;
        let mut cexaey : f64;
        let mut bexdey : f64;
        let mut dexbey : f64;
        let mut alift : f64;
        let mut blift : f64;
        let mut clift : f64;
        let mut dlift : f64;
        let mut ab : f64;
        let mut bc : f64;
        let mut cd : f64;
        let mut da : f64;
        let mut ac : f64;
        let mut bd : f64;
        let mut abc : f64;
        let mut bcd : f64;
        let mut cda : f64;
        let mut dab : f64;
        let mut aezplus : f64;
        let mut bezplus : f64;
        let mut cezplus : f64;
        let mut dezplus : f64;
        let mut aexbeyplus : f64;
        let mut bexaeyplus : f64;
        let mut bexceyplus : f64;
        let mut cexbeyplus : f64;
        let mut cexdeyplus : f64;
        let mut dexceyplus : f64;
        let mut dexaeyplus : f64;
        let mut aexdeyplus : f64;
        let mut aexceyplus : f64;
        let mut cexaeyplus : f64;
        let mut bexdeyplus : f64;
        let mut dexbeyplus : f64;
        let mut det : f64;
        let mut permanent : f64;
        let mut errbound : f64;
        aex = *pa.offset(0isize) - *pe.offset(0isize);
        bex = *pb.offset(0isize) - *pe.offset(0isize);
        cex = *pc.offset(0isize) - *pe.offset(0isize);
        dex = *pd.offset(0isize) - *pe.offset(0isize);
        aey = *pa.offset(1isize) - *pe.offset(1isize);
        bey = *pb.offset(1isize) - *pe.offset(1isize);
        cey = *pc.offset(1isize) - *pe.offset(1isize);
        dey = *pd.offset(1isize) - *pe.offset(1isize);
        aez = *pa.offset(2isize) - *pe.offset(2isize);
        bez = *pb.offset(2isize) - *pe.offset(2isize);
        cez = *pc.offset(2isize) - *pe.offset(2isize);
        dez = *pd.offset(2isize) - *pe.offset(2isize);
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
        permanent = ((cexdeyplus + dexceyplus) * bezplus + (dexbeyplus + bexdeyplus) * cezplus + (bexceyplus + cexbeyplus) * dezplus) * alift + ((dexaeyplus + aexdeyplus) * cezplus + (aexceyplus + cexaeyplus) * dezplus + (cexdeyplus + dexceyplus) * aezplus) * blift + ((aexbeyplus + bexaeyplus) * dezplus + (bexdeyplus + dexbeyplus) * aezplus + (dexaeyplus + aexdeyplus) * bezplus) * clift + ((bexceyplus + cexbeyplus) * aezplus + (cexaeyplus + aexceyplus) * bezplus + (aexbeyplus + bexaeyplus) * cezplus) * dlift;
        errbound = self.isperrboundA * permanent;
        if det > errbound || -det > errbound {
            det
        } else {
            self.insphereadapt(pa,pb,pc,pd,pe,permanent)
        }
    }
}


pub unsafe fn grow_expansion(
    mut elen : i32, mut e : *const f64, mut b : f64, mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut Qnew : f64;
    let mut eindex : i32;
    let mut enow : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    Q = b;
    eindex = 0i32;
    'loop1: loop {
        if !(eindex < elen) {
            break;
        }
        enow = *e.offset(eindex as (isize));
        Qnew = Q + enow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = enow - bvirt;
        around = Q - avirt;
        *h.offset(eindex as (isize)) = around + bround;
        Q = Qnew;
        eindex = eindex + 1;
    }
    *h.offset(eindex as (isize)) = Q;
    eindex + 1i32
}


pub unsafe fn grow_expansion_zeroelim(
    mut elen : i32, mut e : *const f64, mut b : f64, mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut hh : f64;
    let mut Qnew : f64;
    let mut eindex : i32;
    let mut hindex : i32;
    let mut enow : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    hindex = 0i32;
    Q = b;
    eindex = 0i32;
    'loop1: loop {
        if !(eindex < elen) {
            break;
        }
        enow = *e.offset(eindex as (isize));
        Qnew = Q + enow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = enow - bvirt;
        around = Q - avirt;
        hh = around + bround;
        Q = Qnew;
        if hh != 0.0f64 {
            *h.offset(
                 {
                     let _old = hindex;
                     hindex = hindex + 1;
                     _old
                 } as (isize)
             ) = hh;
        }
        eindex = eindex + 1;
    }
    if Q != 0.0f64 || hindex == 0i32 {
        *h.offset(
             {
                 let _old = hindex;
                 hindex = hindex + 1;
                 _old
             } as (isize)
         ) = Q;
    }
    hindex
}


pub unsafe fn expansion_sum(
    mut elen : i32,
    mut e : *const f64,
    mut flen : i32,
    mut f : *const f64,
    mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut Qnew : f64;
    let mut findex : i32;
    let mut hindex : i32;
    let mut hlast : i32;
    let mut hnow : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    Q = *f.offset(0isize);
    hindex = 0i32;
    'loop1: loop {
        if !(hindex < elen) {
            break;
        }
        hnow = *e.offset(hindex as (isize));
        Qnew = Q + hnow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = hnow - bvirt;
        around = Q - avirt;
        *h.offset(hindex as (isize)) = around + bround;
        Q = Qnew;
        hindex = hindex + 1;
    }
    *h.offset(hindex as (isize)) = Q;
    hlast = hindex;
    findex = 1i32;
    'loop3: loop {
        if !(findex < flen) {
            break;
        }
        Q = *f.offset(findex as (isize));
        hindex = findex;
        'loop6: loop {
            if !(hindex <= hlast) {
                break;
            }
            hnow = *h.offset(hindex as (isize));
            Qnew = Q + hnow;
            bvirt = Qnew - Q;
            avirt = Qnew - bvirt;
            bround = hnow - bvirt;
            around = Q - avirt;
            *h.offset(hindex as (isize)) = around + bround;
            Q = Qnew;
            hindex = hindex + 1;
        }
        *h.offset(
             {
                 hlast = hlast + 1;
                 hlast
             } as (isize)
         ) = Q;
        findex = findex + 1;
    }
    hlast + 1i32
}


pub unsafe fn expansion_sum_zeroelim1(
    mut elen : i32,
    mut e : *const f64,
    mut flen : i32,
    mut f : *const f64,
    mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut Qnew : f64;
    let mut index : i32;
    let mut findex : i32;
    let mut hindex : i32;
    let mut hlast : i32;
    let mut hnow : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    Q = *f.offset(0isize);
    hindex = 0i32;
    'loop1: loop {
        if !(hindex < elen) {
            break;
        }
        hnow = *e.offset(hindex as (isize));
        Qnew = Q + hnow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = hnow - bvirt;
        around = Q - avirt;
        *h.offset(hindex as (isize)) = around + bround;
        Q = Qnew;
        hindex = hindex + 1;
    }
    *h.offset(hindex as (isize)) = Q;
    hlast = hindex;
    findex = 1i32;
    'loop3: loop {
        if !(findex < flen) {
            break;
        }
        Q = *f.offset(findex as (isize));
        hindex = findex;
        'loop13: loop {
            if !(hindex <= hlast) {
                break;
            }
            hnow = *h.offset(hindex as (isize));
            Qnew = Q + hnow;
            bvirt = Qnew - Q;
            avirt = Qnew - bvirt;
            bround = hnow - bvirt;
            around = Q - avirt;
            *h.offset(hindex as (isize)) = around + bround;
            Q = Qnew;
            hindex = hindex + 1;
        }
        *h.offset(
             {
                 hlast = hlast + 1;
                 hlast
             } as (isize)
         ) = Q;
        findex = findex + 1;
    }
    hindex = -1i32;
    index = 0i32;
    'loop5: loop {
        if !(index <= hlast) {
            break;
        }
        hnow = *h.offset(index as (isize));
        if hnow != 0.0f64 {
            *h.offset(
                 {
                     hindex = hindex + 1;
                     hindex
                 } as (isize)
             ) = hnow;
        }
        index = index + 1;
    }
    if hindex == -1i32 { 1i32 } else { hindex + 1i32 }
}


pub unsafe fn expansion_sum_zeroelim2(
    mut elen : i32,
    mut e : *const f64,
    mut flen : i32,
    mut f : *const f64,
    mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut hh : f64;
    let mut Qnew : f64;
    let mut eindex : i32;
    let mut findex : i32;
    let mut hindex : i32;
    let mut hlast : i32;
    let mut enow : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    hindex = 0i32;
    Q = *f.offset(0isize);
    eindex = 0i32;
    'loop1: loop {
        if !(eindex < elen) {
            break;
        }
        enow = *e.offset(eindex as (isize));
        Qnew = Q + enow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = enow - bvirt;
        around = Q - avirt;
        hh = around + bround;
        Q = Qnew;
        if hh != 0.0f64 {
            *h.offset(
                 {
                     let _old = hindex;
                     hindex = hindex + 1;
                     _old
                 } as (isize)
             ) = hh;
        }
        eindex = eindex + 1;
    }
    *h.offset(hindex as (isize)) = Q;
    hlast = hindex;
    findex = 1i32;
    'loop3: loop {
        if !(findex < flen) {
            break;
        }
        hindex = 0i32;
        Q = *f.offset(findex as (isize));
        eindex = 0i32;
        'loop6: loop {
            if !(eindex <= hlast) {
                break;
            }
            enow = *h.offset(eindex as (isize));
            Qnew = Q + enow;
            bvirt = Qnew - Q;
            avirt = Qnew - bvirt;
            bround = enow - bvirt;
            around = Q - avirt;
            hh = around + bround;
            Q = Qnew;
            if hh != 0i32 as (f64) {
                *h.offset(
                     {
                         let _old = hindex;
                         hindex = hindex + 1;
                         _old
                     } as (isize)
                 ) = hh;
            }
            eindex = eindex + 1;
        }
        *h.offset(hindex as (isize)) = Q;
        hlast = hindex;
        findex = findex + 1;
    }
    hlast + 1i32
}


pub unsafe fn fast_expansion_sum(
    mut elen : i32,
    mut e : *const f64,
    mut flen : i32,
    mut f : *const f64,
    mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut Qnew : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    let mut eindex : i32;
    let mut findex : i32;
    let mut hindex : i32;
    let mut enow : f64;
    let mut fnow : f64;
    enow = *e.offset(0isize);
    fnow = *f.offset(0isize);
    eindex = {
                 findex = 0i32;
                 findex
             };
    if (fnow > enow) as (i32) == (fnow > -enow) as (i32) {
        Q = enow;
        enow = *e.offset(
                    {
                        eindex = eindex + 1;
                        eindex
                    } as (isize)
                );
    } else {
        Q = fnow;
        fnow = *f.offset(
                    {
                        findex = findex + 1;
                        findex
                    } as (isize)
                );
    }
    hindex = 0i32;
    if eindex < elen && (findex < flen) {
        if (fnow > enow) as (i32) == (fnow > -enow) as (i32) {
            Qnew = enow + Q;
            bvirt = Qnew - enow;
            *h.offset(0isize) = Q - bvirt;
            enow = *e.offset(
                        {
                            eindex = eindex + 1;
                            eindex
                        } as (isize)
                    );
        } else {
            Qnew = fnow + Q;
            bvirt = Qnew - fnow;
            *h.offset(0isize) = Q - bvirt;
            fnow = *f.offset(
                        {
                            findex = findex + 1;
                            findex
                        } as (isize)
                    );
        }
        Q = Qnew;
        hindex = 1i32;
        'loop8: loop {
            if !(eindex < elen && (findex < flen)) {
                break;
            }
            if (fnow > enow) as (i32) == (fnow > -enow) as (i32) {
                Qnew = Q + enow;
                bvirt = Qnew - Q;
                avirt = Qnew - bvirt;
                bround = enow - bvirt;
                around = Q - avirt;
                *h.offset(hindex as (isize)) = around + bround;
                enow = *e.offset(
                            {
                                eindex = eindex + 1;
                                eindex
                            } as (isize)
                        );
            } else {
                Qnew = Q + fnow;
                bvirt = Qnew - Q;
                avirt = Qnew - bvirt;
                bround = fnow - bvirt;
                around = Q - avirt;
                *h.offset(hindex as (isize)) = around + bround;
                fnow = *f.offset(
                            {
                                findex = findex + 1;
                                findex
                            } as (isize)
                        );
            }
            Q = Qnew;
            hindex = hindex + 1;
        }
    }
    'loop9: loop {
        if !(eindex < elen) {
            break;
        }
        Qnew = Q + enow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = enow - bvirt;
        around = Q - avirt;
        *h.offset(hindex as (isize)) = around + bround;
        enow = *e.offset(
                    {
                        eindex = eindex + 1;
                        eindex
                    } as (isize)
                );
        Q = Qnew;
        hindex = hindex + 1;
    }
    'loop10: loop {
        if !(findex < flen) {
            break;
        }
        Qnew = Q + fnow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = fnow - bvirt;
        around = Q - avirt;
        *h.offset(hindex as (isize)) = around + bround;
        fnow = *f.offset(
                    {
                        findex = findex + 1;
                        findex
                    } as (isize)
                );
        Q = Qnew;
        hindex = hindex + 1;
    }
    *h.offset(hindex as (isize)) = Q;
    hindex + 1i32
}


pub unsafe fn fast_expansion_sum_zeroelim(
    mut elen : i32,
    mut e : *const f64,
    mut flen : i32,
    mut f : *const f64,
    mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut Qnew : f64;
    let mut hh : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    let mut eindex : i32;
    let mut findex : i32;
    let mut hindex : i32;
    let mut enow : f64;
    let mut fnow : f64;
    enow = *e.offset(0isize);
    fnow = *f.offset(0isize);
    eindex = {
                 findex = 0i32;
                 findex
             };
    if (fnow > enow) as (i32) == (fnow > -enow) as (i32) {
        Q = enow;
        enow = *e.offset(
                    {
                        eindex = eindex + 1;
                        eindex
                    } as (isize)
                );
    } else {
        Q = fnow;
        fnow = *f.offset(
                    {
                        findex = findex + 1;
                        findex
                    } as (isize)
                );
    }
    hindex = 0i32;
    if eindex < elen && (findex < flen) {
        if (fnow > enow) as (i32) == (fnow > -enow) as (i32) {
            Qnew = enow + Q;
            bvirt = Qnew - enow;
            hh = Q - bvirt;
            enow = *e.offset(
                        {
                            eindex = eindex + 1;
                            eindex
                        } as (isize)
                    );
        } else {
            Qnew = fnow + Q;
            bvirt = Qnew - fnow;
            hh = Q - bvirt;
            fnow = *f.offset(
                        {
                            findex = findex + 1;
                            findex
                        } as (isize)
                    );
        }
        Q = Qnew;
        if hh != 0.0f64 {
            *h.offset(
                 {
                     let _old = hindex;
                     hindex = hindex + 1;
                     _old
                 } as (isize)
             ) = hh;
        }
        'loop9: loop {
            if !(eindex < elen && (findex < flen)) {
                break;
            }
            if (fnow > enow) as (i32) == (fnow > -enow) as (i32) {
                Qnew = Q + enow;
                bvirt = Qnew - Q;
                avirt = Qnew - bvirt;
                bround = enow - bvirt;
                around = Q - avirt;
                hh = around + bround;
                enow = *e.offset(
                            {
                                eindex = eindex + 1;
                                eindex
                            } as (isize)
                        );
            } else {
                Qnew = Q + fnow;
                bvirt = Qnew - Q;
                avirt = Qnew - bvirt;
                bround = fnow - bvirt;
                around = Q - avirt;
                hh = around + bround;
                fnow = *f.offset(
                            {
                                findex = findex + 1;
                                findex
                            } as (isize)
                        );
            }
            Q = Qnew;
            if !(hh != 0.0f64) {
                continue;
            }
            *h.offset(
                 {
                     let _old = hindex;
                     hindex = hindex + 1;
                     _old
                 } as (isize)
             ) = hh;
        }
    }
    'loop10: loop {
        if !(eindex < elen) {
            break;
        }
        Qnew = Q + enow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = enow - bvirt;
        around = Q - avirt;
        hh = around + bround;
        enow = *e.offset(
                    {
                        eindex = eindex + 1;
                        eindex
                    } as (isize)
                );
        Q = Qnew;
        if !(hh != 0.0f64) {
            continue;
        }
        *h.offset(
             {
                 let _old = hindex;
                 hindex = hindex + 1;
                 _old
             } as (isize)
         ) = hh;
    }
    'loop11: loop {
        if !(findex < flen) {
            break;
        }
        Qnew = Q + fnow;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = fnow - bvirt;
        around = Q - avirt;
        hh = around + bround;
        fnow = *f.offset(
                    {
                        findex = findex + 1;
                        findex
                    } as (isize)
                );
        Q = Qnew;
        if !(hh != 0.0f64) {
            continue;
        }
        *h.offset(
             {
                 let _old = hindex;
                 hindex = hindex + 1;
                 _old
             } as (isize)
         ) = hh;
    }
    if Q != 0.0f64 || hindex == 0i32 {
        *h.offset(
             {
                 let _old = hindex;
                 hindex = hindex + 1;
                 _old
             } as (isize)
         ) = Q;
    }
    hindex
}


pub unsafe fn linear_expansion_sum(
    mut elen : i32,
    mut e : *const f64,
    mut flen : i32,
    mut f : *const f64,
    mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut q : f64;
    let mut Qnew : f64;
    let mut R : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    let mut eindex : i32;
    let mut findex : i32;
    let mut hindex : i32;
    let mut enow : f64;
    let mut fnow : f64;
    let mut g0 : f64;
    enow = *e.offset(0isize);
    fnow = *f.offset(0isize);
    eindex = {
                 findex = 0i32;
                 findex
             };
    if (fnow > enow) as (i32) == (fnow > -enow) as (i32) {
        g0 = enow;
        enow = *e.offset(
                    {
                        eindex = eindex + 1;
                        eindex
                    } as (isize)
                );
    } else {
        g0 = fnow;
        fnow = *f.offset(
                    {
                        findex = findex + 1;
                        findex
                    } as (isize)
                );
    }
    if eindex < elen && (findex >= flen || (fnow > enow) as (i32) == (fnow > -enow) as (i32)) {
        Qnew = enow + g0;
        bvirt = Qnew - enow;
        q = g0 - bvirt;
        enow = *e.offset(
                    {
                        eindex = eindex + 1;
                        eindex
                    } as (isize)
                );
    } else {
        Qnew = fnow + g0;
        bvirt = Qnew - fnow;
        q = g0 - bvirt;
        fnow = *f.offset(
                    {
                        findex = findex + 1;
                        findex
                    } as (isize)
                );
    }
    Q = Qnew;
    hindex = 0i32;
    'loop7: loop {
        if !(hindex < elen + flen - 2i32) {
            break;
        }
        if eindex < elen && (findex >= flen || (fnow > enow) as (i32) == (fnow > -enow) as (i32)) {
            R = enow + q;
            bvirt = R - enow;
            *h.offset(hindex as (isize)) = q - bvirt;
            enow = *e.offset(
                        {
                            eindex = eindex + 1;
                            eindex
                        } as (isize)
                    );
        } else {
            R = fnow + q;
            bvirt = R - fnow;
            *h.offset(hindex as (isize)) = q - bvirt;
            fnow = *f.offset(
                        {
                            findex = findex + 1;
                            findex
                        } as (isize)
                    );
        }
        Qnew = Q + R;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = R - bvirt;
        around = Q - avirt;
        q = around + bround;
        Q = Qnew;
        hindex = hindex + 1;
    }
    *h.offset(hindex as (isize)) = q;
    *h.offset((hindex + 1i32) as (isize)) = Q;
    hindex + 2i32
}


pub unsafe fn linear_expansion_sum_zeroelim(
    mut elen : i32,
    mut e : *const f64,
    mut flen : i32,
    mut f : *const f64,
    mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut q : f64;
    let mut hh : f64;
    let mut Qnew : f64;
    let mut R : f64;
    let mut bvirt : f64;
    let mut avirt : f64;
    let mut bround : f64;
    let mut around : f64;
    let mut eindex : i32;
    let mut findex : i32;
    let mut hindex : i32;
    let mut count : i32;
    let mut enow : f64;
    let mut fnow : f64;
    let mut g0 : f64;
    enow = *e.offset(0isize);
    fnow = *f.offset(0isize);
    eindex = {
                 findex = 0i32;
                 findex
             };
    hindex = 0i32;
    if (fnow > enow) as (i32) == (fnow > -enow) as (i32) {
        g0 = enow;
        enow = *e.offset(
                    {
                        eindex = eindex + 1;
                        eindex
                    } as (isize)
                );
    } else {
        g0 = fnow;
        fnow = *f.offset(
                    {
                        findex = findex + 1;
                        findex
                    } as (isize)
                );
    }
    if eindex < elen && (findex >= flen || (fnow > enow) as (i32) == (fnow > -enow) as (i32)) {
        Qnew = enow + g0;
        bvirt = Qnew - enow;
        q = g0 - bvirt;
        enow = *e.offset(
                    {
                        eindex = eindex + 1;
                        eindex
                    } as (isize)
                );
    } else {
        Qnew = fnow + g0;
        bvirt = Qnew - fnow;
        q = g0 - bvirt;
        fnow = *f.offset(
                    {
                        findex = findex + 1;
                        findex
                    } as (isize)
                );
    }
    Q = Qnew;
    count = 2i32;
    'loop7: loop {
        if !(count < elen + flen) {
            break;
        }
        if eindex < elen && (findex >= flen || (fnow > enow) as (i32) == (fnow > -enow) as (i32)) {
            R = enow + q;
            bvirt = R - enow;
            hh = q - bvirt;
            enow = *e.offset(
                        {
                            eindex = eindex + 1;
                            eindex
                        } as (isize)
                    );
        } else {
            R = fnow + q;
            bvirt = R - fnow;
            hh = q - bvirt;
            fnow = *f.offset(
                        {
                            findex = findex + 1;
                            findex
                        } as (isize)
                    );
        }
        Qnew = Q + R;
        bvirt = Qnew - Q;
        avirt = Qnew - bvirt;
        bround = R - bvirt;
        around = Q - avirt;
        q = around + bround;
        Q = Qnew;
        if hh != 0i32 as (f64) {
            *h.offset(
                 {
                     let _old = hindex;
                     hindex = hindex + 1;
                     _old
                 } as (isize)
             ) = hh;
        }
        count = count + 1;
    }
    if q != 0i32 as (f64) {
        *h.offset(
             {
                 let _old = hindex;
                 hindex = hindex + 1;
                 _old
             } as (isize)
         ) = q;
    }
    if Q != 0.0f64 || hindex == 0i32 {
        *h.offset(
             {
                 let _old = hindex;
                 hindex = hindex + 1;
                 _old
             } as (isize)
         ) = Q;
    }
    hindex
}


pub unsafe fn compress(
    mut elen : i32, mut e : *const f64, mut h : *mut f64
) -> i32 {
    let mut Q : f64;
    let mut q : f64;
    let mut Qnew : f64;
    let mut eindex : i32;
    let mut hindex : i32;
    let mut bvirt : f64;
    let mut enow : f64;
    let mut hnow : f64;
    let mut top : i32;
    let mut bottom : i32;
    bottom = elen - 1i32;
    Q = *e.offset(bottom as (isize));
    eindex = elen - 2i32;
    'loop1: loop {
        if !(eindex >= 0i32) {
            break;
        }
        enow = *e.offset(eindex as (isize));
        Qnew = Q + enow;
        bvirt = Qnew - Q;
        q = enow - bvirt;
        if q != 0i32 as (f64) {
            *h.offset(
                 {
                     let _old = bottom;
                     bottom = bottom - 1;
                     _old
                 } as (isize)
             ) = Qnew;
            Q = q;
        } else {
            Q = Qnew;
        }
        eindex = eindex - 1;
    }
    top = 0i32;
    hindex = bottom + 1i32;
    'loop3: loop {
        if !(hindex < elen) {
            break;
        }
        hnow = *h.offset(hindex as (isize));
        Qnew = hnow + Q;
        bvirt = Qnew - hnow;
        q = Q - bvirt;
        if q != 0i32 as (f64) {
            *h.offset(
                 {
                     let _old = top;
                     top = top + 1;
                     _old
                 } as (isize)
             ) = q;
        }
        Q = Qnew;
        hindex = hindex + 1;
    }
    *h.offset(top as (isize)) = Q;
    top + 1i32
}


pub unsafe fn estimate(
    mut elen : i32, mut e : *const f64
) -> f64 {
    let mut Q : f64;
    let mut eindex : i32;
    Q = *e.offset(0isize);
    eindex = 1i32;
    'loop1: loop {
        if !(eindex < elen) {
            break;
        }
        Q = Q + *e.offset(eindex as (isize));
        eindex = eindex + 1;
    }
    Q
}


pub unsafe fn orient2dfast(
    mut pa : *const f64, mut pb : *const f64, mut pc : *const f64
) -> f64 {
    let mut acx : f64;
    let mut bcx : f64;
    let mut acy : f64;
    let mut bcy : f64;
    acx = *pa.offset(0isize) - *pc.offset(0isize);
    bcx = *pb.offset(0isize) - *pc.offset(0isize);
    acy = *pa.offset(1isize) - *pc.offset(1isize);
    bcy = *pb.offset(1isize) - *pc.offset(1isize);
    acx * bcy - acy * bcx
}


pub unsafe fn orient3dfast(
    mut pa : *const f64,
    mut pb : *const f64,
    mut pc : *const f64,
    mut pd : *const f64
) -> f64 {
    let mut adx : f64;
    let mut bdx : f64;
    let mut cdx : f64;
    let mut ady : f64;
    let mut bdy : f64;
    let mut cdy : f64;
    let mut adz : f64;
    let mut bdz : f64;
    let mut cdz : f64;
    adx = *pa.offset(0isize) - *pd.offset(0isize);
    bdx = *pb.offset(0isize) - *pd.offset(0isize);
    cdx = *pc.offset(0isize) - *pd.offset(0isize);
    ady = *pa.offset(1isize) - *pd.offset(1isize);
    bdy = *pb.offset(1isize) - *pd.offset(1isize);
    cdy = *pc.offset(1isize) - *pd.offset(1isize);
    adz = *pa.offset(2isize) - *pd.offset(2isize);
    bdz = *pb.offset(2isize) - *pd.offset(2isize);
    cdz = *pc.offset(2isize) - *pd.offset(2isize);
    adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady) + cdx * (ady * bdz - adz * bdy)
}



pub unsafe fn incirclefast(
    mut pa : *const f64,
    mut pb : *const f64,
    mut pc : *const f64,
    mut pd : *const f64
) -> f64 {
    let mut adx : f64;
    let mut ady : f64;
    let mut bdx : f64;
    let mut bdy : f64;
    let mut cdx : f64;
    let mut cdy : f64;
    let mut abdet : f64;
    let mut bcdet : f64;
    let mut cadet : f64;
    let mut alift : f64;
    let mut blift : f64;
    let mut clift : f64;
    adx = *pa.offset(0isize) - *pd.offset(0isize);
    ady = *pa.offset(1isize) - *pd.offset(1isize);
    bdx = *pb.offset(0isize) - *pd.offset(0isize);
    bdy = *pb.offset(1isize) - *pd.offset(1isize);
    cdx = *pc.offset(0isize) - *pd.offset(0isize);
    cdy = *pc.offset(1isize) - *pd.offset(1isize);
    abdet = adx * bdy - bdx * ady;
    bcdet = bdx * cdy - cdx * bdy;
    cadet = cdx * ady - adx * cdy;
    alift = adx * adx + ady * ady;
    blift = bdx * bdx + bdy * bdy;
    clift = cdx * cdx + cdy * cdy;
    alift * bcdet + blift * cadet + clift * abdet
}


pub unsafe fn inspherefast(
    mut pa : *const f64,
    mut pb : *const f64,
    mut pc : *const f64,
    mut pd : *const f64,
    mut pe : *const f64
) -> f64 {
    let mut aex : f64;
    let mut bex : f64;
    let mut cex : f64;
    let mut dex : f64;
    let mut aey : f64;
    let mut bey : f64;
    let mut cey : f64;
    let mut dey : f64;
    let mut aez : f64;
    let mut bez : f64;
    let mut cez : f64;
    let mut dez : f64;
    let mut alift : f64;
    let mut blift : f64;
    let mut clift : f64;
    let mut dlift : f64;
    let mut ab : f64;
    let mut bc : f64;
    let mut cd : f64;
    let mut da : f64;
    let mut ac : f64;
    let mut bd : f64;
    let mut abc : f64;
    let mut bcd : f64;
    let mut cda : f64;
    let mut dab : f64;
    aex = *pa.offset(0isize) - *pe.offset(0isize);
    bex = *pb.offset(0isize) - *pe.offset(0isize);
    cex = *pc.offset(0isize) - *pe.offset(0isize);
    dex = *pd.offset(0isize) - *pe.offset(0isize);
    aey = *pa.offset(1isize) - *pe.offset(1isize);
    bey = *pb.offset(1isize) - *pe.offset(1isize);
    cey = *pc.offset(1isize) - *pe.offset(1isize);
    dey = *pd.offset(1isize) - *pe.offset(1isize);
    aez = *pa.offset(2isize) - *pe.offset(2isize);
    bez = *pb.offset(2isize) - *pe.offset(2isize);
    cez = *pc.offset(2isize) - *pe.offset(2isize);
    dez = *pd.offset(2isize) - *pe.offset(2isize);
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
    dlift * abc - clift * dab + (blift * cda - alift * bcd)
}

