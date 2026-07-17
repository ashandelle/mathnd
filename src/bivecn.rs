use std::{iter::Sum, ops::{Add, AddAssign, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign}};

use num_traits::{One, Zero};
use rand::prelude::*;
use rand_distr::StandardNormal;

use crate::{matn::MatN, traits::{Sqrt, Two}, vecn::VecN};

#[derive(Debug, Clone, Copy)]
pub struct BiVecN<T, const N: usize> {
    pub m: MatN<T, N>,
}

impl<T, const N: usize> PartialEq for BiVecN<T, N> where T: PartialEq {
    fn eq(&self, other: &Self) -> bool {
        self.m == other.m
    }
}

// Unary minus
impl<T, const N: usize> Neg for BiVecN<T, N> where T: Neg<Output = T> + Copy {
    type Output = BiVecN<T, N>;
    fn neg(self) -> BiVecN<T, N> {
        BiVecN {
            m: -self.m,
        }
    }
}

// Vector addition
impl<T, const N: usize> Add for BiVecN<T, N> where T: Add<Output = T> + Copy {
    type Output = BiVecN<T, N>;
    fn add(self, v: BiVecN<T, N>) -> BiVecN<T, N> {
        BiVecN {
            m: self.m + v.m,
        }
    }
}

impl<T, const N: usize> AddAssign for BiVecN<T, N> where T: AddAssign + Copy {
    fn add_assign(&mut self, v: BiVecN<T, N>) {
        self.m += v.m;
    }
}

// Vector subtraction
impl<T, const N: usize> Sub for BiVecN<T, N> where T: Sub<Output = T> + Copy {
    type Output = BiVecN<T, N>;
    fn sub(self, v: BiVecN<T, N>) -> BiVecN<T, N> {
        BiVecN {
            m: self.m - v.m,
        }
    }
}

impl<T, const N: usize> SubAssign for BiVecN<T, N> where T: SubAssign + Copy {
    fn sub_assign(&mut self, v: BiVecN<T, N>) {
        self.m -= v.m;
    }
}

// Scalar multiplication
impl<T, const N: usize> Mul<T> for BiVecN<T, N> where T: Mul<Output = T> + Copy {
    type Output = BiVecN<T, N>;
    fn mul(self, s: T) -> BiVecN<T, N> {
        BiVecN {
            m: self.m * s,
        }
    }
}

impl<T, const N: usize> MulAssign<T> for BiVecN<T, N> where T: MulAssign + Copy {
    fn mul_assign(&mut self, s: T) {
        self.m *= s;
    }
}

// Scalar division
impl<T, const N: usize> Div<T> for BiVecN<T, N> where T: Div<Output = T> + Copy {
    type Output = BiVecN<T, N>;
    fn div(self, s: T) -> BiVecN<T, N> {
        BiVecN {
            m: self.m / s,
        }
    }
}

impl<T, const N: usize> DivAssign<T> for BiVecN<T, N> where T: DivAssign + Copy {
    fn div_assign(&mut self, s: T) {
        self.m /= s;
    }
}

impl<T, const N: usize> Sum for BiVecN<T, N> where T: Zero + Copy {
    fn sum<I: Iterator<Item = Self>>(iter: I) -> Self {
        let mut total: BiVecN<T, N> = Self::zero();
        for item in iter {
            total = total + item;
        }
        total
    }
}

// impl fmt::Display for BiVecN {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         // write!(f, "{:?}", self.e)
//         let dim = self.m.e.len();
//         write!(f, "[");
//         for i in 0..dim {
//             for j in (i+1)..dim {
//                 write!(f, "{}", self.m.e[i].e[j]);
//                 if (i != dim-2) || (j != dim-1) {
//                     write!(f, ", ");
//                 }
//             }
//         }
//         write!(f, "]")
//     }
// }

impl<T, const N: usize> From<MatN<T, N>> for BiVecN<T, N> {
    fn from(mat: MatN<T, N>) -> Self {
        Self {
            m: mat,
        }
    }
}

impl<T, const N: usize> Into<MatN<T, N>> for BiVecN<T, N> {
    fn into(self) -> MatN<T, N> {
        self.m
    }
}

impl<T, const N: usize> BiVecN<T, N> {
    pub fn new(m: MatN<T, N>) -> Self {
        BiVecN {
            m: m,
        }
    }

    pub fn rand_normal(rng: &mut ThreadRng) -> Self where T: Neg<Output = T> + Zero + Copy, StandardNormal: Distribution<T> {
        let mut m = MatN::new(std::array::from_fn(|i| VecN::new(std::array::from_fn(|j| if i < j {rng.sample(StandardNormal)} else {T::zero()}))));
        for i in 0..N {
            for j in 0..i {
                m.e[i].e[j] = -m.e[j].e[i];
            }
        }
        Self::new(m)
    }

    // Dot product
    pub fn dot(&self, b: BiVecN<T, N>) -> T where T: Mul<Output = T> + Div<Output = T> + Sum + Two + Copy {
        self.m.dot(&b.m) / T::two()
    }

    // Length
    pub fn length(&self) -> T where T: Mul<Output = T> + Div<Output = T> + Sqrt + Sum + Two + Copy {
        (self.m.length_sqr() / T::two()).sqrt()
    }

    // Length squared
    pub fn length_sqr(&self) -> T where T: Mul<Output = T> + Div<Output = T> + Sum + Two + Copy {
        self.m.length_sqr() / T::two()
    }

    // Skew
    pub fn skew(&self) -> BiVecN<T, N> where T: Sub<Output = T> + Div<Output = T> + Two + Copy {
        BiVecN {
            m: (self.m - self.m.transposed()) / T::two(),
        }
    }

    // pub fn exponential(&self) -> MatN<T, N> where T: Copy {
    // }

    pub fn get(&self, i: usize, j: usize) -> T where T: Copy {
        self.m.get(i, j)
    }

    // To MatN
    pub fn to_matn(self) -> MatN<T, N> {
        self.m
    }

    // Zero
    pub fn zero() -> Self where T: Zero + Copy {
        Self {
            m: MatN::zero(),
        }
    }

    // Basis element
    pub fn basis(i: usize, j: usize) -> Self where T: Neg<Output = T> + Zero + One + Copy {
        let mut mat = MatN::zero();
        if i != j {
            mat.e[i].e[j] = T::one();
            mat.e[j].e[i] = -T::one();
        }
        Self {
            m: mat,
        }
    }
}