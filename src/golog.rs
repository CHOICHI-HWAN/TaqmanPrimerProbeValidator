//! Go's math.Log, bit for bit. The Go engine's numbers are the reference, and the platform `ln` differs from Go's in
//! the last bit for about 2.5 % of inputs, which can move a printed Tm or dG by 0.01.
//!
//! Translated from Go's src/math/log.go (Copyright 2009 The Go Authors, BSD-3-Clause; the licence text is in
//! THIRD_PARTY_NOTICES.md). The Go file's algorithm and constants are from FreeBSD's e_log.c, which came with this
//! notice:
//!
//! ====================================================
//! Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
//!
//! Developed at SunPro, a Sun Microsystems, Inc. business.
//! Permission to use, copy, modify, and distribute this
//! software is freely granted, provided that this notice
//! is preserved.
//! ====================================================

/// Go's math.Log (the FreeBSD e_log.c algorithm).
pub fn go_log(x: f64) -> f64 {
    const LN2_HI: f64 = 6.93147180369123816490e-01;
    const LN2_LO: f64 = 1.90821492927058770002e-10;
    const L1: f64 = 6.666666666666735130e-01;
    const L2: f64 = 3.999999999940941908e-01;
    const L3: f64 = 2.857142874366239149e-01;
    const L4: f64 = 2.222219843214978396e-01;
    const L5: f64 = 1.818357216161805012e-01;
    const L6: f64 = 1.531383769920937332e-01;
    const L7: f64 = 1.479819860511658591e-01;
    if x.is_nan() || x == f64::INFINITY {
        return x;
    }
    if x < 0.0 {
        return f64::NAN;
    }
    if x == 0.0 {
        return f64::NEG_INFINITY;
    }
    let (mut f1, mut ki) = frexp(x);
    if f1 < std::f64::consts::SQRT_2 / 2.0 {
        f1 *= 2.0;
        ki -= 1;
    }
    let f = f1 - 1.0;
    let k = ki as f64;
    let s = f / (2.0 + f);
    let s2 = s * s;
    let s4 = s2 * s2;
    let t1 = s2 * (L1 + s4 * (L3 + s4 * (L5 + s4 * L7)));
    let t2 = s4 * (L2 + s4 * (L4 + s4 * L6));
    let r = t1 + t2;
    let hfsq = 0.5 * f * f;
    k * LN2_HI - ((hfsq - (s * (hfsq + r) + k * LN2_LO)) - f)
}

/// Go's math.Frexp for finite non-zero x: frac in [0.5, 1).
fn frexp(mut x: f64) -> (f64, i32) {
    let mut exp = 0i32;
    if x.abs() < f64::MIN_POSITIVE {
        x *= (1u64 << 52) as f64;
        exp = -52;
    }
    let mut b = x.to_bits();
    exp += ((b >> 52) & 0x7ff) as i32 - 1022;
    b &= !(0x7ffu64 << 52);
    b |= 1022u64 << 52;
    (f64::from_bits(b), exp)
}

#[cfg(test)]
mod tests {
    use super::go_log;

    #[test]
    fn known_values() {
        // Go: math.Log(0.865335013001561) = -0.14463854876317195 (libm gives ...192)
        assert_eq!(go_log(0.865335013001561), -0.14463854876317195);
        assert_eq!(go_log(1.0), 0.0);
        assert_eq!(go_log(std::f64::consts::E), 1.0);
    }
}
