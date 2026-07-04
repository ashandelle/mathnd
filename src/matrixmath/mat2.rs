use std::ops::Neg;

use num_traits::Zero;

use crate::{matn::MatN, traits::Trig, vecn::VecN};

impl<T> MatN<T, 2> {
    pub fn skew_exponential_2d(&self) -> Self where T: Neg<Output = T> + Trig + Copy {
        let theta: T = self.e[0].e[1];
        let cos = theta.cos();
        let sin = theta.sin();
        MatN::new([
            VecN::new([cos, sin]),
            VecN::new([-sin, cos]),
        ])
    }

    pub fn ortho_logarithm_2d(&self) -> Self where T: Neg<Output = T> + Trig + Zero + Copy {
        let theta = self.e[0].e[1].atan2(self.e[0].e[0]);
        MatN::new([
            VecN::new([T::zero(), theta]),
            VecN::new([-theta, T::zero()]),
        ])
    }
}