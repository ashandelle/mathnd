use std::{iter::Sum, ops::{Div, Mul, Sub}};

use num_traits::{FromPrimitive, One, Zero};

use crate::{matn::MatN, traits::{Sqrt, Trig, Two}};

impl<T> MatN<T, 3> {
    pub fn skew_exponential_3d(&self, eps: T) -> Self where T: Mul<Output = T> + Sub<Output = T> + Div<Output = T> + Sum + Sqrt + PartialOrd + Zero + One + Trig + Copy {
        let lensqr: T = self.e[0].e[1]*self.e[0].e[1] + self.e[0].e[2]*self.e[0].e[2] + self.e[1].e[2]*self.e[1].e[2];
        let len = lensqr.sqrt();

        if len < eps {
            return MatN::identity();
        }

        MatN::identity() + *self * (len.sin() / len) + (*self * *self) * ((T::one() - len.cos()) / lensqr)
    }

    pub fn ortho_logarithm_3d(&self, eps: T) -> Self where T: Mul<Output = T> + Sub<Output = T> + Div<Output = T> + Sum + Sqrt + PartialOrd + Zero + One + Two + Trig + FromPrimitive + Copy {
        let trace: T = self.trace();

        let theta = ((trace - T::one()) / T::two()).acos();

        if theta < eps {
            return MatN::zero();
        }

        let mut s = trace - T::one();
        s = s * s;
        s = (T::from_usize(4).unwrap() - s).sqrt();

        (*self - self.transposed()) * (theta / s)
    }
}