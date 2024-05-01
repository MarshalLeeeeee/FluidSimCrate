/// Module for integration
use ndarray as nd;

pub trait Integrate {
    fn get_integration_step(&self, d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64>;
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

pub fn get_integration_step(integ: Box<dyn Integrate>, d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64> {
    integ.get_integration_step(&d, &dd, step)
}

pub fn integrate(integ: Box<dyn Integrate>, d: &nd::Array2::<f64>, dd: &nd::Array2::<f64>, step: f64) -> nd::Array2::<f64> {
    integ.integrate(&d, &dd, step)
}