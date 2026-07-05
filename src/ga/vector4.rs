use std::ops::Mul;

use crate::ga::{Bivector4, Reflector4, Rotor4, Trivector4};

pub struct Vector4 {
    pub e1: f64,
    pub e2: f64,
    pub e3: f64,
    pub e4: f64,
}

impl Mul<Vector4> for Vector4 {
    type Output = (f64, Bivector4);
    fn mul(self, v: Vector4) -> (f64, Bivector4) {
        
    }
}

impl Mul<Bivector4> for Vector4 {
    type Output = (Vector4, Trivector4);
    fn mul(self, b: Bivector4) -> (Vector4, Trivector4) {
        
    }
}

impl Mul<Trivector4> for Vector4 {
    type Output = (Bivector4, f64);
    fn mul(self, t: Trivector4) -> (Bivector4, f64) {
        
    }
}

impl Mul<Rotor4> for Vector4 {
    type Output = Reflector4;
    fn mul(self, r: Rotor4) -> Reflector4 {
        
    }
}

impl Mul<Reflector4> for Vector4 {
    type Output = Rotor4;
    fn mul(self, r: Reflector4) -> Rotor4 {
        
    }
}