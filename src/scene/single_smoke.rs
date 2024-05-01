/// Module for single smoke simulation
///
/// The scene has rigid boundary without any other abstacles
use ndarray as nd;
use serde_json::{Value, Map, Number};
use crate::canvas;
use crate::grid::*;
use crate::parser;
use crate::boundary;
use crate::integration;
use crate::laplacian_solver;

/// Smoke simulation performed with Euler method
///
/// Examples can be found in ```src/examples/```
pub struct SingleSmokeGridScene {
    /// discrete of time
    dt: f64,
    /// discrete of space
    ds: f64,
    /// gravity
    g: nd::Array1::<f64>,
    /// density of air
    rho_air: f64,
    /// density of smoke
    rho_smoke: f64,
    /// color of air for visualization
    c_air: canvas::RGBAColor,
    /// color of smoke for visualization
    c_smoke: canvas::RGBAColor,
    /// portion of density of smoke, range from [0, 1], scalar field
    sf_d: nd::Array2::<f64>,
    /// velocity of axis x, staggered grid, vector field
    vf_vx: nd::Array2::<f64>,
    /// velocity of axis y, staggered grid, vector field
    vf_vy: nd::Array2::<f64>,
    /// pressure solver
    pressure_solver: laplacian_solver::LaplacianSolver,
}
impl SingleSmokeGridScene {
    
    /// Get the config table
    pub fn get_config_map() -> Map<String, Value> {
        let mut m = Map::new();
        m.insert(String::from("dt"), Value::Number(Number::from_f64(0.5_f64).unwrap()));
        m.insert(String::from("ds"), Value::Number(Number::from_f64(1.0_f64).unwrap()));
        m.insert(String::from("gx"), Value::Number(Number::from_f64(0_f64).unwrap()));
        m.insert(String::from("gy"), Value::Number(Number::from_f64(9.8_f64).unwrap()));
        m.insert(String::from("rho_air"), Value::Number(Number::from_f64(1_f64).unwrap()));
        m.insert(String::from("rho_smoke"), Value::Number(Number::from_f64(0.5_f64).unwrap()));
        m.insert(String::from("c_air_r"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_air_g"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_air_b"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_smoke_r"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_smoke_g"), Value::Number(Number::from(0_u8)));
        m.insert(String::from("c_smoke_b"), Value::Number(Number::from(0_u8)));
        m
    }

    /// Create instance from parser
    pub fn new_by_parser(parser: &Value) -> Self {
        let dt = parser::get_from_parser_f64(parser, "dt");
        let ds = parser::get_from_parser_f64(parser, "ds");
        let gx = parser::get_from_parser_f64(parser, "gx");
        let gy = parser::get_from_parser_f64(parser, "gy");
        let rho_air = parser::get_from_parser_f64(parser, "rho_air");
        let rho_smoke = parser::get_from_parser_f64(parser, "rho_smoke");
        let c_air_r = parser::get_from_parser_u8(parser, "c_air_r");
        let c_air_g = parser::get_from_parser_u8(parser, "c_air_g");
        let c_air_b = parser::get_from_parser_u8(parser, "c_air_b");
        let c_smoke_r = parser::get_from_parser_u8(parser, "c_smoke_r");
        let c_smoke_g = parser::get_from_parser_u8(parser, "c_smoke_g");
        let c_smoke_b = parser::get_from_parser_u8(parser, "c_smoke_b");
        let w = parser::get_from_parser_usize(parser, "width");
        let h = parser::get_from_parser_usize(parser, "height");
        Self {
            dt: dt,
            ds: ds,
            g: nd::Array1::<f64>::from_vec(vec![gx, gy]),
            rho_air: rho_air,
            rho_smoke: rho_smoke,
            c_air: canvas::RGBAColor::new(c_air_r, c_air_g, c_air_b, 255_u8),
            c_smoke: canvas::RGBAColor::new(c_smoke_r, c_smoke_g, c_smoke_b, 255_u8),
            sf_d: nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{
                if 16 <= i && i <= 48 && 1 <= j && j <= 32 { 1_f64 }
                else { 0_f64 }
            }),
            vf_vx: nd::Array2::<f64>::from_elem((w+1, h), 0_f64),
            vf_vy: nd::Array2::<f64>::from_elem((w, h+1), 0_f64),
            pressure_solver: laplacian_solver::LaplacianSolver::new(w, h),
        }
    }

    fn _solve_rho(&self) -> nd::Array2::<f64> {
        self.sf_d.clone() * (self.rho_smoke - self.rho_air) + self.rho_air
    }
    
    fn _solve_force(&self, rho: &nd::Array2::<f64>) -> (nd::Array2::<f64>, nd::Array2::<f64>) {
        let rho_x = to_staggered_x_grid(&rho, Box::new(boundary::NeumannBoundary));
        let fx = (self.rho_air - rho_x.clone()) / rho_x * self.g[0];
        
        let rho_y = to_staggered_y_grid(&rho, Box::new(boundary::NeumannBoundary));
        let fy = (self.rho_air - rho_y.clone()) / rho_y * self.g[1];
    
        (fx, fy)
    }
    
    fn _velocity_advection(&mut self, force_x: nd::Array2::<f64>, force_y: nd::Array2::<f64>) {
        let vx_g = staggered_x_grid_gradient_as_staggered_x(&self.vf_vx, self.ds);
        let vx_pre = self.vf_vx.clone() + force_x * self.dt - vx_g * self.dt;
        let vx_pre_g = staggered_x_grid_gradient_as_staggered_x(&vx_pre, self.ds);
        let vx_new = (self.vf_vx.clone() + vx_pre) * 0.5 - vx_pre_g * 0.5 * self.dt;
        self.vf_vx.assign(&vx_new);
        
        let vy_g = staggered_y_grid_gradient_as_staggered_y(&self.vf_vy, self.ds);
        let vy_pre = self.vf_vy.clone() + force_y * self.dt - vy_g * self.dt;
        let vy_pre_g = staggered_y_grid_gradient_as_staggered_y(&vy_pre, self.ds);
        let vy_new = (self.vf_vy.clone() + vy_pre) * 0.5 - vy_pre_g * 0.5 * self.dt;
        self.vf_vy.assign(&vy_new);
    }
    
    fn _velocity_boundary_adjust(&mut self) {
        let (w, h) = self.vf_vx.dim();
        let w = w - 1;
        for j in 0..h {
            self.vf_vx[[0, j]] = 0_f64;
            self.vf_vx[[w, j]] = 0_f64;
        }
        for i in 0..w {
            self.vf_vy[[i, 0]] = 0_f64;
            self.vf_vy[[i, h]] = 0_f64;
        }
    }
    
    fn _solve_velocity_divergence(&self) -> nd::Array2::<f64> {
        staggered_x_grid_gradient(&self.vf_vx, self.ds) + staggered_y_grid_gradient(&self.vf_vy, self.ds)
    }
    
    fn _pressure_projection(&mut self, rho: &nd::Array2::<f64>) {
        let (w, h) = rho.dim();
        let sz = w * h;

        // pressure solver
        let divergence = self._solve_velocity_divergence();
        let b = divergence * self.ds * self.ds;
        let b = b.into_shape(sz).unwrap();
        let p = self.pressure_solver.solve(b).into_shape((w, h)).unwrap();
    
        // pressure projection
        let px = to_staggered_x_grid(&p, Box::new(boundary::NeumannBoundary));
        let pxg = staggered_x_grid_gradient_as_staggered_x(&px, self.ds);
        self.vf_vx = self.vf_vx.clone() - pxg;
        
        let py = to_staggered_y_grid(&p, Box::new(boundary::NeumannBoundary));
        let pyg = staggered_y_grid_gradient_as_staggered_y(&py, self.ds);
        self.vf_vy = self.vf_vy.clone() - pyg;
    }
    
    fn _density_advection(&mut self) {
        let (w, h) = self.sf_d.dim();
        
        let (vx_0, vx_1) = split_staggered_x_grid(&self.vf_vx);
        let vx_mid = (vx_0 + vx_1) * 0.5_f64;
        let xs = nd::Array2::<f64>::from_shape_fn((w, h), |(i, _j)|{ i as f64 });
        let xs = integration::integrate(Box::new(integration::ForwardEuler), &xs, &vx_mid, -self.dt);
        
        let (vy_0, vy_1) = split_staggered_y_grid(&self.vf_vy);
        let vy_mid = (vy_0 + vy_1) * 0.5_f64;
        let ys = nd::Array2::<f64>::from_shape_fn((w, h), |(_i, j)|{ j as f64 });
        let ys = integration::integrate(Box::new(integration::ForwardEuler), &ys, &vy_mid, -self.dt);

        self.sf_d.assign(&sample_with_spatial_index(&self.sf_d, xs, ys));
    }

    /// Do simulation in the temporal step
    pub fn sim(&mut self) {
        let rho = self._solve_rho();
        let (fx, fy) = self._solve_force(&rho);
        self._velocity_advection(fx, fy);
        self._velocity_boundary_adjust();
        self._pressure_projection(&rho);
        self._velocity_boundary_adjust();
        self._density_advection();
    }

    /// Visualization of smoke density
    pub fn visualize_density(&self) -> Vec<canvas::RGBAColor> {
        let (w, h) = self.sf_d.dim();
        let mut buffer = Vec::with_capacity(w * h);
        for j in 0..h {
            for i in 0..w {
                buffer.push(canvas::RGBAColor::mix(&self.c_smoke, &self.c_air, self.sf_d[[i, j]]));
            }
        }
        buffer.reverse();
        buffer
    }
}
