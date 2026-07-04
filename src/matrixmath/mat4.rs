use std::{iter::Sum, ops::{Div, Mul, Neg, Sub, SubAssign}};

use num_traits::{FromPrimitive, One, Signed, Zero};

use crate::{matn::MatN, traits::Sqrt, vecn::VecN};

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

    // pub fn skew_exponential_4d(&self, eps: T) -> Self where T: Copy {
        
    // }

    // pub fn ortho_logarithm_4d(&self, eps: T) -> Self where T: Neg<Output = T> + Mul<Output = T> + Div<Output = T> + Sub<Output = T> + Sum + SubAssign + PartialOrd + Signed + Sqrt + Zero + One + FromPrimitive + Copy {
    //     let (x, y) = self.isoclinic_decomposition(eps);
    // }
}