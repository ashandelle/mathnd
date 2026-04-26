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
    use crate::{matn::MatN, vecn::VecN};

    use rand::prelude::*;
    use rand_distr::StandardNormal;

    #[test]
    fn determinant() {
        let mut rng = rand::rng();

        let identity: MatN<f64, 4> = MatN::identity();

        assert_eq!(identity.determinant(1e-12), 1.0);

        for _n in 0..100 {
            let rand1: MatN<f64, 4> = MatN {
                    e: std::array::from_fn(|_i|
                        VecN::new(std::array::from_fn(|_j| rng.sample(StandardNormal)))
                    ),
                };
            
            let rand2: MatN<f64, 4> = MatN {
                    e: std::array::from_fn(|_i|
                        VecN::new(std::array::from_fn(|_j| rng.sample(StandardNormal)))
                    ),
                };

            let trans = rand1.transposed();

            assert!((rand1.determinant(1e-12) - trans.determinant(1e-12)).abs() < 1e-12);

            assert!(((rand1.determinant(1e-12) * rand2.determinant(1e-12)) - (rand1 * rand2).determinant(1e-12)).abs() < 1e-12);

            let diag: MatN<f64, 4> = MatN {
                    e: std::array::from_fn(|i|
                        VecN::new(std::array::from_fn(|j| if i == j {rng.sample(StandardNormal)} else {0.0}))
                    ),
                };

            let d: f64 = (0..4).into_iter().map(|i| diag.e[i].e[i]).product();
            assert!((diag.determinant(1e-12) - d).abs() < 1e-12);

            let tria: MatN<f64, 4> = MatN {
                    e: std::array::from_fn(|i|
                        VecN::new(std::array::from_fn(|j| if i <= j {rng.sample(StandardNormal)} else {0.0}))
                    ),
                };

            let d: f64 = (0..4).into_iter().map(|i| tria.e[i].e[i]).product();
            assert!((tria.determinant(1e-12) - d).abs() < 1e-12);
        }
    }

    #[test]
    fn orthogonal_product() {
        let mut arr: [VecN<f64, 4>; 3] = [
                VecN::basis(0),
                VecN::basis(1),
                VecN::basis(2),
            ];
        
        assert_eq!(VecN::orthogonal_product(&arr, 1e-12), VecN::basis(3));

        arr.swap(0, 1);
        
        assert_eq!(VecN::orthogonal_product(&arr, 1e-12), -VecN::basis(3));
    }
}
