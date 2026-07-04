use std::{iter::Sum, ops::{Div, Mul, Neg, Sub, SubAssign}};

use num_traits::{FromPrimitive, One, Signed, Zero};

use crate::{matn::MatN, traits::{Sqrt, Trig, Two}, vecn::VecN};

use std::fmt::Debug;

impl<T> MatN<T, 4> {
    pub fn isoclinic_decomposition(&self, eps: T) -> (MatN<T, 4>, MatN<T, 4>) where T: Neg<Output = T> + Mul<Output = T> + Div<Output = T> + Sub<Output = T> + Sum + SubAssign + PartialOrd + Signed + Sqrt + Zero + One + FromPrimitive + Copy {
        // let a1: MatN<T, 4>  = MatN::new([
        //     VecN::new([T::zero(), T::zero(), T::zero(),-T::one()]),
        //     VecN::new([T::zero(), T::zero(),-T::one(), T::zero()]),
        //     VecN::new([T::zero(), T::one(), T::zero(), T::zero()]),
        //     VecN::new([T::one(), T::zero(), T::zero(), T::zero()]),
        // ]);

        // let a2: MatN<T, 4>  = MatN::new([
        //     VecN::new([ T::zero(), T::zero(), T::one(), T::zero()]),
        //     VecN::new([ T::zero(), T::zero(), T::zero(),-T::one()]),
        //     VecN::new([-T::one(), T::zero(), T::zero(), T::zero()]),
        //     VecN::new([ T::zero(), T::one(), T::zero(), T::zero()]),
        // ]);

        // let a3: MatN<T, 4>  = MatN::new([
        //     VecN::new([T::zero(),-T::one(), T::zero(), T::zero()]),
        //     VecN::new([T::one(), T::zero(), T::zero(), T::zero()]),
        //     VecN::new([T::zero(), T::zero(), T::zero(),-T::one()]),
        //     VecN::new([T::zero(), T::zero(), T::one(), T::zero()]),
        // ]);

        // let mut s = *self - a1**self*a1 - a2**self*a2 - a3**self*a3;

        let m = MatN::new([
            VecN::new([
                self.e[0].e[0]+self.e[1].e[1]+self.e[2].e[2]+self.e[3].e[3],
                self.e[0].e[1]-self.e[2].e[3]-self.e[1].e[0]+self.e[3].e[2],
                self.e[0].e[2]+self.e[1].e[3]-self.e[2].e[0]-self.e[3].e[1],
                self.e[0].e[3]-self.e[1].e[2]-self.e[3].e[0]+self.e[2].e[1]
            ]),
            VecN::new([
                self.e[1].e[0]+self.e[2].e[3]-self.e[3].e[2]-self.e[0].e[1],
                self.e[1].e[1]+self.e[2].e[2]+self.e[3].e[3]+self.e[0].e[0],
                self.e[1].e[2]-self.e[2].e[1]+self.e[3].e[0]-self.e[0].e[3],
                self.e[1].e[3]-self.e[2].e[0]-self.e[3].e[1]+self.e[0].e[2]
            ]),
            VecN::new([
                self.e[2].e[0]-self.e[1].e[3]-self.e[0].e[2]+self.e[3].e[1],
                self.e[2].e[1]-self.e[1].e[2]+self.e[0].e[3]-self.e[3].e[0],
                self.e[2].e[2]+self.e[1].e[1]+self.e[0].e[0]+self.e[3].e[3],
                self.e[2].e[3]+self.e[1].e[0]-self.e[0].e[1]-self.e[3].e[2]
            ]),
            VecN::new([
                self.e[3].e[0]-self.e[0].e[3]+self.e[1].e[2]-self.e[2].e[1],
                self.e[3].e[1]-self.e[0].e[2]-self.e[1].e[3]+self.e[2].e[0],
                self.e[3].e[2]+self.e[0].e[1]-self.e[1].e[0]-self.e[2].e[3],
                self.e[3].e[3]+self.e[0].e[0]+self.e[1].e[1]+self.e[2].e[2]
            ]),
        ]);

        // let y = s / s.determinant(eps).sqrt().sqrt();
        let y = m / m.determinant(eps).sqrt().sqrt();

        let x = *self * y.transposed();

        (x, y)
    }

    pub fn skew_exponential_4d(&self, eps: T) -> Self where T: Neg<Output = T> + Mul<Output = T> + Sub<Output = T> + Div<Output = T> + Sum + PartialOrd + Sqrt + Zero + One + Two + Trig + Copy {
        let mut a = VecN::new([ // Left quaternion
            T::zero(),
            self.e[0].e[1] + self.e[2].e[3],
            self.e[0].e[2] - self.e[1].e[3],
            self.e[0].e[3] + self.e[1].e[2],
        ]);
        a = a / T::two();

        let len = a.length(); // Exponentiate left quaternion
        a = if len < eps {VecN::zero()} else {a * len.sin() / len};
        a.e[0] = len.cos();

        let a = MatN::new([ // Left isoclinic matrix
            a,
            VecN::new([-a.e[1], a.e[0], a.e[3],-a.e[2]]),
            VecN::new([-a.e[2],-a.e[3], a.e[0], a.e[1]]),
            VecN::new([-a.e[3], a.e[2],-a.e[1], a.e[0]]),
        ]);

        let mut b = VecN::new([ // Right quaternion
            T::zero(),
            self.e[0].e[1] - self.e[2].e[3],
            self.e[0].e[2] + self.e[1].e[3],
            self.e[0].e[3] - self.e[1].e[2],
        ]);
        b = b / T::two();

        let len = b.length(); // Exponentiate right quaternion
        b = if len < eps {VecN::zero()} else {b * len.sin() / len};
        b.e[0] = len.cos();

        let b = MatN::new([ // Right isoclinic matrix
            b,
            VecN::new([-b.e[1], b.e[0],-b.e[3], b.e[2]]),
            VecN::new([-b.e[2], b.e[3], b.e[0],-b.e[1]]),
            VecN::new([-b.e[3],-b.e[2], b.e[1], b.e[0]]),
        ]);

        a * b
    }

    // pub fn skew_exponential_4d_b(&self, eps: T) -> Self where T: Neg<Output = T> + Sub<Output = T> + Div<Output = T> + Sum + PartialOrd + Sqrt + Zero + One + Two + Trig + Copy {
    //     let mut a = VecN::new([ // Left quaternion
    //         T::zero(),
    //         self.e[0].e[1] + self.e[2].e[3],
    //         self.e[0].e[2] - self.e[1].e[3],
    //         self.e[0].e[3] + self.e[1].e[2],
    //     ]);
    //     a = a / T::two();

    //     let mut a = MatN::new([ // Left isoclinic matrix
    //         a,
    //         VecN::new([-a.e[1], T::zero(), a.e[3],-a.e[2]]),
    //         VecN::new([-a.e[2],-a.e[3], T::zero(), a.e[1]]),
    //         VecN::new([-a.e[3], a.e[2],-a.e[1], T::zero()]),
    //     ]);

    //     let lensqr = a.length_sqr() / T::two(); // Exponentiate left matrix
    //     let len = lensqr.sqrt();

    //     if len < eps {
    //         a = MatN::identity();
    //     } else {
    //         a = MatN::identity() + a * (len.sin() / len) + (a * a) * ((T::one() - len.cos()) / lensqr);
    //     }

    //     let mut b = VecN::new([ // Right quaternion
    //         T::zero(),
    //         self.e[0].e[1] - self.e[2].e[3],
    //         self.e[0].e[2] + self.e[1].e[3],
    //         self.e[0].e[3] - self.e[1].e[2],
    //     ]);
    //     b = b / T::two();

    //     let mut b = MatN::new([ // Right isoclinic matrix
    //         b,
    //         VecN::new([-b.e[1], T::zero(),-b.e[3], b.e[2]]),
    //         VecN::new([-b.e[2], b.e[3], T::zero(),-b.e[1]]),
    //         VecN::new([-b.e[3],-b.e[2], b.e[1], T::zero()]),
    //     ]);

    //     let lensqr = b.length_sqr() / T::two(); // Exponentiate right matrix
    //     let len = lensqr.sqrt();

    //     if len < eps {
    //         b = MatN::identity();
    //     } else {
    //         b = MatN::identity() + b * (len.sin() / len) + (b * b) * ((T::one() - len.cos()) / lensqr);
    //     }
        
    //     a*b
    // }

    pub fn ortho_logarithm_4d(&self, eps: T) -> Self where T: Debug + Neg<Output = T> + Mul<Output = T> + Div<Output = T> + Sub<Output = T> + Sum + SubAssign + PartialOrd + Signed + Sqrt + Zero + One + Trig + FromPrimitive + Copy {
        let (x, y) = self.isoclinic_decomposition(eps);

        let mut a = x.e[0]; // Left quaternion
        let mut b = y.e[0]; // Right quaternion

        // Logarithm of a and b
        let s = a.e[0];
        a.e[0] = T::zero();
        let l = a.length();
        a = if l < eps {VecN::zero()} else {a * s.acos() / l};

        let s = b.e[0];
        b.e[0] = T::zero();
        let l = b.length();
        b = if l < eps {VecN::zero()} else {b * s.acos() / l};

        let x = MatN::new([
            a,
            // VecN::new([T::zero(), a.e[1], a.e[2], a.e[3]]),
            VecN::new([-a.e[1], T::zero(), a.e[3],-a.e[2]]),
            VecN::new([-a.e[2],-a.e[3], T::zero(), a.e[1]]),
            VecN::new([-a.e[3], a.e[2],-a.e[1], T::zero()]),
        ]);

        let y = MatN::new([
            b,
            // VecN::new([T::zero(), b.e[1], b.e[2], b.e[3]]),
            VecN::new([-b.e[1], T::zero(),-b.e[3], b.e[2]]),
            VecN::new([-b.e[2], b.e[3], T::zero(),-b.e[1]]),
            VecN::new([-b.e[3],-b.e[2], b.e[1], T::zero()]),
        ]);

        x + y
    }
}