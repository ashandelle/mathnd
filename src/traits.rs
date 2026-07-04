use std::ops::Add;

use num_traits::{One};

pub trait Sqrt{
    fn sqrt(self) -> Self;
}

macro_rules! impl_sqrt {
    ( $($ty:ty),* ) => {
        $(
            impl Sqrt for $ty {
                fn sqrt(self) -> Self {
                    self.sqrt()
                }
            }
        )*
    };
}

impl_sqrt!(f32, f64);

pub trait Trig{
    fn cos(self) -> Self;
    fn sin(self) -> Self;
    fn tan(self) -> Self;

    fn acos(self) -> Self;
    fn asin(self) -> Self;
    fn atan(self) -> Self;

    fn cosh(self) -> Self;
    fn sinh(self) -> Self;
    fn tanh(self) -> Self;

    fn acosh(self) -> Self;
    fn asinh(self) -> Self;
    fn atanh(self) -> Self;

    fn atan2(self, other: Self) -> Self;
}

macro_rules! impl_trig {
    ( $($ty:ty),* ) => {
        $(
            impl Trig for $ty {
                fn cos(self) -> Self {
                    self.cos()
                }
                fn sin(self) -> Self {
                    self.sin()
                }
                fn tan(self) -> Self {
                    self.tan()
                }

                fn acos(self) -> Self {
                    self.acos()
                }
                fn asin(self) -> Self {
                    self.asin()
                }
                fn atan(self) -> Self {
                    self.atan()
                }

                fn cosh(self) -> Self {
                    self.cosh()
                }
                fn sinh(self) -> Self {
                    self.sinh()
                }
                fn tanh(self) -> Self {
                    self.tanh()
                }

                fn acosh(self) -> Self {
                    self.acosh()
                }
                fn asinh(self) -> Self {
                    self.asinh()
                }
                fn atanh(self) -> Self {
                    self.atanh()
                }

                fn atan2(self, other: Self) -> Self {
                    self.atan2(other)
                }
            }
        )*
    };
}

impl_trig!(
    f32, f64
);

pub trait MinMax{
    fn min(self, other: Self) -> Self;
    fn max(self, other: Self) -> Self;
}

macro_rules! impl_minmax {
    ( $($ty:ty),* ) => {
        $(
            impl MinMax for $ty {
                fn min(self, other: Self) -> Self {
                    std::cmp::min(self, other)
                }

                fn max(self, other: Self) -> Self {
                    std::cmp::max(self, other)
                }
            }
        )*
    };
}

macro_rules! impl_minmaxf {
    ( $($ty:ty),* ) => {
        $(
            impl MinMax for $ty {
                fn min(self, other: Self) -> Self {
                    self.min(other)
                }

                fn max(self, other: Self) -> Self {
                    self.max(other)
                }
            }
        )*
    };
}

impl_minmax!(
    isize, usize,
    i8, i16, i32, i64,
    u8, u16, u32, u64
);

impl_minmaxf!(f32, f64);

pub trait Two{
    fn two() -> Self;
}

impl<T: Add<Output = T> + One> Two for T {
    fn two() -> Self {
        T::one() + T::one()
    }
}

// pub trait MinMaxValue{
//     fn minimum() -> Self;
//     fn maximum() -> Self;
// }

// macro_rules! impl_minmaxvalue {
//     ( $( $ty:ty ),* ) => {
//         $(
//             impl MinMaxValue for $ty {
//                 fn minimum() -> Self {
//                     Self::MIN
//                 }

//                 fn maximum() -> Self {
//                     Self::MAX
//                 }
//             }
//         )*
//     };
// }

// impl_minmaxvalue!(
//     isize, usize,
//     i8, i16, i32, i64,
//     u8, u16, u32, u64,
//     f32, f64
// );