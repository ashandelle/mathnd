#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

// pub fn add(left: u64, right: u64) -> u64 {
//     left + right
// }

pub mod vecn;
pub mod matn;
pub mod bivecn;
pub mod traits;

#[cfg(test)]
mod tests {
    use super::*;

    // #[test]
    // fn it_works() {
    //     let result = add(2, 2);
    //     assert_eq!(result, 4);
    // }
}
