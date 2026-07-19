use std::{iter::Sum, ops::{Div, Mul, Neg, Sub}};
use std::fmt::Debug;

use bit_iter::BitIter;
use num_traits::{FromPrimitive, One, Signed, Zero};

use crate::{matn::MatN, traits::{Sqrt, Two}, util::factorial};

impl<T, const N: usize> MatN<T, N> {
    pub fn orthonormalized(&self, eps: T, max: usize) -> MatN<T, N> where
        T: Sub<Output = T> + Mul<Output = T> + Div<Output = T> + PartialOrd + Sum + Sqrt + Signed + Zero + FromPrimitive + Copy {
        let mut mat = self.normalized_basis();

        let n = T::from_usize((N - 1).max(2)).unwrap();
        let nn = T::from_usize((N*N - N) / 2).unwrap();

        let mut dot: T = T::zero();
        let mut iter = 0;

        while iter == 0 || (dot > eps && iter < max) {
            let mut tmp = mat.clone();

            dot = T::zero();
            for i in 0..N {
                for j in (i+1)..N {
                    let d = mat.e[i].dot(mat.e[j]);
                    dot = dot + d.abs();
                    tmp.e[i] = tmp.e[i] - (mat.e[j] * d / n);
                    tmp.e[j] = tmp.e[j] - (mat.e[i] * d / n);
                }
            }
            dot = dot / nn;

            mat = tmp.normalized_basis();
            iter += 1;
        }

        mat
    }

    // pub fn inverse(&self, eps: T, max: usize) -> MatN<T, N> where
    //     T: Debug + Mul<Output = T> + Sub<Output = T> + Div<Output = T> + PartialOrd + Sum + Zero + One + Two + FromPrimitive + Copy {
    //     let mut inv: MatN<T, N> = MatN::identity();

    //     let mut iter = 0;
    //     let mut len = T::zero();
    //     let identity = MatN::identity();

    //     let nn = T::from_usize(N*N).unwrap();

    //     while iter == 0 || (len > eps && iter < max) {
    //         let mut id = inv * *self;
    //         inv = (inv * T::two()) - (id * inv);
    //         id = id - identity;
    //         len = id.dot(&id) / nn;
    //         iter+=1;
    //         println!("{:?}", inv);
    //     }

    //     inv
    // }

    pub fn determinant(&self, eps: T) -> T where T: Neg<Output = T> + PartialOrd + Signed + One + Copy {
        let mut determinant = T::one();
        let mut swapcount = 0;
        let mut tmp = self.clone();

        for i in 0..N {
            let mut maxrow = i;
            for k in (i+1)..N {
                if tmp.e[k].e[i].abs() > tmp.e[maxrow].e[i].abs() {
                    maxrow = k;
                }
            }

            if maxrow != i {
                swapcount += 1;
                tmp.e.swap(i, maxrow);
            }

            if tmp.e[i].e[i].abs() < eps {
                return T::zero();
            }

            determinant = determinant * tmp.e[i].e[i];

            for k in (i+1)..N {
                let factor = tmp.e[k].e[i] / tmp.e[i].e[i];
                for j in (i+1)..N {
                    tmp.e[k].e[j] = tmp.e[k].e[j] - factor * tmp.e[i].e[j];
                }
            }
        }

        if swapcount % 2 != 0 {
            determinant = -determinant;
        }

        determinant
    }

    pub fn flip_if_negative(&mut self, eps: T) where T: Neg<Output = T> + PartialOrd + Signed + One + Copy {
        let d = self.determinant(eps);

        if d < T::zero() {
            self.e[N - 1] = -self.e[N - 1];
        }
    }

    pub fn flip_if_positive(&mut self, eps: T) where T: Neg<Output = T> + PartialOrd + Signed + One + Copy {
        let d = self.determinant(eps);

        if d > T::zero() {
            self.e[N - 1] = -self.e[N - 1];
        }
    }

    pub fn exponential_taylor(&self, iter: usize) -> Self where T: Mul<Output = T> + Div<Output = T> + Sum + Zero + One + FromPrimitive + Copy {
        let mut out: MatN<T, N> = Self::identity();

        let mut pows: Vec<MatN<T, N>> = vec![*self];

        for i in 1..(iter.ilog2() as usize) {
            let m = *pows.get(i-1).unwrap();
            pows.push(m * m);
        }

        for i in 1..iter {
            let pow = BitIter::from(i).map(|j| pows[j]).product::<MatN<T, N>>();
            out = out + pow / factorial(i);
        }

        out
    }

    pub fn logarithm_taylor(&self, iter: usize) -> Self where T: Neg<Output = T> + Mul<Output = T> + Sub<Output = T> + Div<Output = T> + Sum + Zero + One + FromPrimitive + Copy {
        let mut out: MatN<T, N> = Self::zero();

        let mut pows: Vec<MatN<T, N>> = vec![Self::identity() - *self];

        for i in 1..(iter.ilog2() as usize) {
            let m = *pows.get(i-1).unwrap();
            pows.push(m * m);
        }

        for i in 1..iter {
            let pow = BitIter::from(i).map(|j| pows[j]).product::<MatN<T, N>>();
            out = out + pow / T::from_usize(i).unwrap();
        }

        -out
    }

    // pub fn skew_exponential(&self, eps: T) -> Self where Assert::<{ N > 4 }>: True, T: Zero + One + Copy {
    //     MatN::identity()
    // }

    // pub fn ortho_logarithm(&self, eps: T) -> Self where Assert::<{ N > 4 }>: True, T: Zero + One + Copy {
    //     MatN::identity()
    // }
}