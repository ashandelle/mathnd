use std::ops::Mul;

use num_traits::{FromPrimitive};

pub fn factorial<T>(i: usize) -> T where T: Mul<Output = T> + FromPrimitive {
    (1..=i).map(|j| T::from_usize(j).unwrap()).reduce(|acc: T, j: T| {
        acc * j
    }).unwrap()
}

pub trait True {}
pub trait False {}

pub struct Assert<const COND: bool> {}
impl True for Assert<true> {}
impl False for Assert<false> {}