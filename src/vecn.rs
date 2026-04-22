use std::{iter::Sum, ops::{Add, AddAssign, BitXor, Div, DivAssign, Mul, MulAssign, Neg, Sub, SubAssign}};

use num_traits::{One, Signed, Zero};

use crate::{bivecn::BiVecN, matn::MatN, traits::{CosSin, Sqrt, Two}};

#[derive(Debug, Clone, Copy)]
pub struct VecN<T, const N: usize> {
    pub e: [T; N]
}

impl<T, const N: usize> PartialEq for VecN<T, N> where T: PartialEq {
    fn eq(&self, other: &Self) -> bool {
        self.e == other.e
    }
}

impl<T, const N: usize> Neg for VecN<T, N> where T: Neg<Output = T> + Copy {
    type Output = VecN<T, N>;
    fn neg(self) -> VecN<T, N> {
        VecN::new(std::array::from_fn(|i| -self.e[i]))
    }
}

impl<T, const N: usize> Add for VecN<T, N> where T: Add<Output = T> + Copy {
    type Output = VecN<T, N>;
    fn add(self, v: VecN<T, N>) -> VecN<T, N> {
        VecN::new(std::array::from_fn(|i| self.e[i] + v.e[i]))
    }
}

impl<T, const N: usize> AddAssign for VecN<T, N> where T: AddAssign + Copy {
    fn add_assign(&mut self, v: VecN<T, N>) {
        for (i, val) in self.e.iter_mut().enumerate() {
            *val += v.e[i];
        }
    }
}

impl<T, const N: usize> Sub for VecN<T, N> where T: Sub<Output = T> + Copy {
    type Output = VecN<T, N>;
    fn sub(self, v: VecN<T, N>) -> VecN<T, N> {
        VecN::new(std::array::from_fn(|i| self.e[i] - v.e[i]))
    }
}

impl<T, const N: usize> SubAssign for VecN<T, N> where T: SubAssign + Copy {
    fn sub_assign(&mut self, v: VecN<T, N>) {
        for (i, val) in self.e.iter_mut().enumerate() {
            *val -= v.e[i];
        }
    }
}

// impl<T, const N: usize> Mul<VecN<T, N>> for T where T: Mul<Output = T> + Copy {
//     type Output = VecN<T, N>;
//     fn mul(self, v: VecN<T, N>) -> VecN<T, N> {
//         VecN::new(std::array::from_fn(|i| self * v.e[i]))
//     }
// }
impl<T, const N: usize> Mul<T> for VecN<T, N> where T: Mul<Output = T> + Copy {
    type Output = VecN<T, N>;
    fn mul(self, s: T) -> VecN<T, N> {
        VecN::new(std::array::from_fn(|i| self.e[i] * s))
    }
}

impl<T, const N: usize> MulAssign<T> for VecN<T, N> where T: MulAssign + Copy {
    fn mul_assign(&mut self, s: T) {
        for val in self.e.iter_mut() {
            *val *= s;
        }
    }
}

impl<T, const N: usize> Div<T> for VecN<T, N> where T: Div<Output = T> + Copy {
    type Output = VecN<T, N>;
    fn div(self, s: T) -> VecN<T, N> {
        VecN::new(std::array::from_fn(|i| self.e[i] / s))
    }
}

impl<T, const N: usize> DivAssign<T> for VecN<T, N> where T: DivAssign + Copy {
    fn div_assign(&mut self, s: T) {
        for val in self.e.iter_mut() {
            *val /= s;
        }
    }
}

impl<T, const N: usize> BitXor for VecN<T, N> where T: Sub<Output = T> + Mul<Output = T> + Div<Output = T> + Two + Copy {
    type Output = BiVecN<T, N>;
    fn bitxor(self, v: VecN<T, N>) -> BiVecN<T, N> {
        BiVecN {
            m: MatN::mult_transpose_vecn(self, v) - MatN::mult_transpose_vecn(v, self),
        }.skew()
    }
}

// impl fmt::Display for VecN {
//     fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
//         write!(f, "{:?}", self.e)
//     }
// }

impl<T, const N: usize> VecN<T, N> {
    pub fn new(e: [T; N]) -> Self {
        VecN {
            e: e,
        }
    }

    // Dot product
    pub fn dot(&self, v: VecN<T, N>) -> T where T: Mul<Output = T> + Sum + Copy {
        (self.e).iter()
                .zip((v.e).iter())
                .map(|(&x, &y)| x * y)
                .sum()
    }

    // Left contraction
    pub fn left_contract(&self, v: BiVecN<T, N>) -> VecN<T, N> where T: Neg<Output = T> + Mul<Output = T> + Sum + Copy {
        -(v.m * *self)
    }

    // Normalize
    pub fn normalize(&mut self) where T: Mul<Output = T> + DivAssign + Sqrt + Sum + Copy {
        *self /= self.length();
    }

    pub fn normalized(&self) -> VecN<T, N> where T: Mul<Output = T> + Div<Output = T> + Sqrt + Sum + Copy {
        *self / self.length()
    }

    // Length
    pub fn length(&self) -> T where T: Mul<Output = T> + Sqrt + Sum + Copy {
        (self.e).iter()
                .map(|&x| x*x)
                .sum::<T>().sqrt()
    }

    // Length squared
    pub fn length_sqr(&self) -> T where T: Mul<Output = T> + Sum + Copy {
        (self.e).iter()
                .map(|&x| x*x)
                .sum::<T>()
    }

    pub fn orthonormal_basis(&self) -> [VecN<T, N>; N-1] where T: Neg<Output = T> + Sub<Output = T> + Mul<Output = T> + Div<Output = T> + AddAssign + SubAssign + PartialOrd + Sum + Sqrt + Signed + Zero + One + Copy {
        let normal: VecN<T, N> = self.normalized();
    
        let mut vecs: [VecN<T, N>; N-1] = [VecN::zero(); N-1];
        let mut maxdot: T = T::zero();
        let mut maxi: usize = 0;
    
        for i in 0..N {
            let d = normal.e[i].abs();
            if d > maxdot {
                maxdot = d;
                maxi = i;
            }
        }
        
        let mut l = 0;
        for i in 0..N {
            if i != maxi {
                let mut v: VecN<T, N> = normal * -normal.e[i];
                v.e[i] += T::one();
                vecs[l] = v;
                l+=1;
            }
        }
    
        for j in 0..(N-1) {
            let vec = vecs[j].normalized();
            for k in j+1..(N-1) {
                vecs[k] -= vec * vec.dot(vecs[k]);
            }
        }
    
        vecs
    }

    pub fn orthogonal_product(vecs: &[VecN<T, N>; N-1], eps: T) -> VecN<T, N> where T: Neg<Output = T> + SubAssign + PartialOrd + Signed + Zero + One + Copy {
        let mut mat = MatN::new(std::array::from_fn(|i| if i < N-1 {vecs[i]} else {VecN::zero()}));

        VecN::new(std::array::from_fn(|i| {
            mat.e[N-1] = VecN::basis(i);
            mat.determinant(eps)
        }))
    }

    pub fn rotate(&self, i: usize, j: usize, angle: T) -> VecN<T, N> where T: Mul<Output = T> + Add<Output = T> + Sub<Output = T> + CosSin + Copy {
        let mut r = self.clone();

        let a = r.e[i];
        let b = r.e[j];

        r.e[i] = a * angle.cos() + b * angle.sin();
        r.e[j] = b * angle.cos() - a * angle.sin();

        r
    }

    // Zero
    pub fn zero() -> Self where T: Zero + Copy {
        Self::new([T::zero(); N])
    }

    pub fn default() -> Self where T: Default + Copy {
        Self::new([T::default(); N])
    }

    // Basis element
    pub fn basis(element: usize) -> Self where T: Zero + One + Copy {
        let mut vec = Self::zero();
        vec.e[element] = T::one();
        vec
    }
}