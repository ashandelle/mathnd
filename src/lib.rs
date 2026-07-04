#![feature(generic_const_exprs)]
#![allow(incomplete_features)]

// pub fn add(left: u64, right: u64) -> u64 {
//     left + right
// }

pub mod util;
pub mod vecn;
pub mod matn;
pub mod bivecn;
pub mod traits;

#[cfg(test)]
mod tests {
    use crate::{bivecn::BiVecN, matn::MatN, vecn::VecN};

    use rand::prelude::*;
    use rand_distr::StandardNormal;

    #[test]
    fn determinant() {
        let mut rng = rand::rng();

        let identity: MatN<f64, 4> = MatN::identity();

        assert_eq!(identity.determinant(1e-12), 1.0);

        for _n in 0..1000 {
            let rand1: MatN<f64, 4> = MatN::rand_normal(&mut rng);
            
            let rand2: MatN<f64, 4> = MatN::rand_normal(&mut rng);

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

    #[test]
    fn bivec_reflection() {
        let mut rng = rand::rng();

        for _n in 0..1000 {
            let rand1: VecN<f64, 4> = VecN::rand_normal(&mut rng);
            let rand2: VecN<f64, 4> = VecN::rand_normal(&mut rng);

            let rand3: VecN<f64, 4> = VecN::rand_normal(&mut rng);

            let reflect1 = rand3.reflect(rand1);
            let reflect2 = rand3.reflect(rand2);

            let mut bivec1 = rand1 ^ rand2;
            let bivec2 = reflect1 ^ reflect2;

            bivec1 = rand3.reflect_bivec(bivec1);

            let dif = bivec1 - bivec2;
            let dot = dif.dot(dif);
            assert!(dot < 1e-12);
        }
    }

    #[test]
    fn reflection() {
        let mut rng = rand::rng();

        for _n in 0..1000 {
            let rand1: VecN<f64, 4> = VecN::rand_normal(&mut rng);
            let rand2: VecN<f64, 4> = VecN::rand_normal(&mut rng);

            let mat = rand1.reflect_mat(MatN::identity());

            let reflect1 = rand1.reflect(rand2);
            let reflect2 = mat * rand2;

            assert!((reflect1 - reflect2).length_sqr() < 1e-12);
        }
    }

    #[test]
    fn explog_taylor() {
        let mut rng = rand::rng();

        for _n in 0..1000 {
            let mat: MatN<f64, 4> = MatN::rand_normal(&mut rng) * 0.1;

            let exp = mat.exponential_taylor(32);

            let log = exp.logarithm_taylor(128);

            assert!((mat - log).length_sqr() < 1e-8);
        }
    }

    #[test]
    fn skewortho_explog() {
        let mut rng = rand::rng();

        for n in 0..1000 {
            let mat: MatN<f64, 2> = if n == 0 {MatN::zero()} else {BiVecN::rand_normal(&mut rng).to_matn()};

            let expt = mat.exponential_taylor(32);
            let exp1 = mat.skew_exponential(1e-8);

            assert!((exp1 - expt).length_sqr() < 1e-8);

            let log = exp1.ortho_logarithm(1e-8);
            let exp2 = log.skew_exponential(1e-8);

            assert!((exp1 - exp2).length_sqr() < 1e-8);
        }

        for n in 0..1000 {
            let mat: MatN<f64, 3> = if n == 0 {MatN::zero()} else {BiVecN::rand_normal(&mut rng).to_matn()};

            let expt = mat.exponential_taylor(64);
            let exp1 = mat.skew_exponential(1e-8);

            assert!((exp1 - expt).length_sqr() < 1e-8);

            let log = exp1.ortho_logarithm(1e-8);
            let exp2 = log.skew_exponential(1e-8);

            assert!((exp1 - exp2).length_sqr() < 1e-8);
        }

        for n in 0..2 {
            // let mat: MatN<f64, 4> = if n == 0 {MatN::zero()} else {BiVecN::rand_normal(&mut rng).to_matn()};

            // let expt = mat.exponential_taylor(64);
            // let exp1 = mat.skew_exponential(1e-8);

            // assert!((exp1 - expt).length_sqr() < 1e-8);

            // let log = exp1.ortho_logarithm(1e-8);
            // let exp2 = log.skew_exponential(1e-8);

            // assert!((exp1 - exp2).length_sqr() < 1e-8);
        }
    }

    #[test]
    fn isoclinic_decomposition() {
        let mut rng = rand::rng();

        for n in 0..1000 {
            let mat: MatN<f64, 4> = if n == 0 {MatN::identity()} else {MatN::rand_normal(&mut rng).orthonormalized(1e-8, 32)};

            let (x, y) = mat.isoclinic_decomposition(1e-8);

            assert!((mat - x*y).length_sqr() < 1e-8);
        }
    }
}
