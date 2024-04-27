/// Module for boundary conditions

use ndarray as nd;

pub trait Boundary {
    /// Return the conditioned value according to the boundary value 
    fn boundary(&self, d: &nd::Array1::<f64>) -> nd::Array1::<f64>;
}

/// Dirichlet boundary condition
///
/// Take one argument as the boundary value
pub struct DirichletBoundary(pub f64);
impl Boundary for DirichletBoundary {
    fn boundary(&self, d: &nd::Array1::<f64>) -> nd::Array1::<f64> {
        nd::Array1::<f64>::from_elem(d.dim(), self.0)
    }
}

/// Neumann boundary condition
///
/// Implement zero 1st derivative neumann boundary condition, i.e., the same as the boundary value
pub struct NeumannBoundary;
impl Boundary for NeumannBoundary {
    fn boundary(&self, d: &nd::Array1::<f64>) -> nd::Array1::<f64> {
        d.clone()
    }
}
