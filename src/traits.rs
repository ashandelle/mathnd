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

pub trait CosSin{
    fn cos(self) -> Self;
    fn sin(self) -> Self;
}

macro_rules! impl_cossin {
    ( $($ty:ty),* ) => {
        $(
            impl CosSin for $ty {
                fn cos(self) -> Self {
                    self.cos()
                }
                fn sin(self) -> Self {
                    self.sin()
                }
            }
        )*
    };
}

impl_cossin!(
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