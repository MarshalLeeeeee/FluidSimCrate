/// Module for integration
use ndarray as nd;
use crate::grid::sample_with_spatial_index;

pub trait Integrate {

    /// Get the incremental part for integration
    ///
    /// Params:
    /// - ```d```: data to be integrated.
    /// - ```dd```: the first derivative of ```d```.
    /// - ```step```: discrete step for integration.
    fn get_integration_step(&self, d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64>;
    
    /// Get the new data after integration
    ///
    /// Params:
    /// - ```d```: data to be integrated.
    /// - ```dd```: the first derivative of ```d```.
    /// - ```step```: discrete step for integration.
    fn integrate(&self, d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64> {
        d + self.get_integration_step(&d, &dd, step)
    }
}

pub struct ForwardEuler;
impl Integrate for ForwardEuler {
    fn get_integration_step(&self, _d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64> {
        dd * step
    }
}

pub struct RK2;
impl Integrate for RK2 {
    fn get_integration_step(&self, _d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64> {
        dd * step
    }
}

/// Get the incremental part for integration
///
/// Params:
/// - ```integ```: specific integration type.
/// - ```d```: data to be integrated.
/// - ```dd```: the first derivative of ```d```.
/// - ```step```: discrete step for integration.
pub fn get_integration_step(integ: Box<dyn Integrate>, d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64> {
    integ.get_integration_step(&d, &dd, step)
}

/// Get the new data after integration
///
/// Params:
/// - ```d```: data to be integrated.
/// - ```dd```: the first derivative of ```d```.
/// - ```step```: discrete step for integration.
pub fn integrate(integ: Box<dyn Integrate>, d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64> {
    integ.integrate(&d, &dd, step)
}