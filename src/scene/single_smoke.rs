/// Module for single smoke simulation
///
/// The scene has rigid boundary without any other abstacles
use ndarray as nd;
use serde_json::{Value, Map, Number};
use crate::canvas;
use crate::grid::*;
use crate::parser;
use crate::render;
use crate::boundary;
use crate::laplacian_solver;

/// Smoke simulation performed with Euler method, using symplectic Euler integration
///
/// Examples can be found in ```src/examples/```
pub struct SingleSmokeGridScene<C: canvas::Color> {
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
    c_air: C,
    /// color of smoke for visualization
    c_smoke: C,
    /// portion of density of smoke, range from [0, 1], scalar field
    sf_d: nd::Array2::<f64>,
    /// initial smoke pos x min
    init_d_x_min: usize,
    /// initial smoke pos x max
    init_d_x_max: usize,
    /// initial smoke pos y min
    init_d_y_min: usize,
    /// initial smoke pos y max
    init_d_y_max: usize,
    /// velocity of axis x, staggered grid, vector field
    vf_vx: nd::Array2::<f64>,
    /// velocity of axis y, staggered grid, vector field
    vf_vy: nd::Array2::<f64>,
    /// pressure
    sf_p: nd::Array2::<f64>,
    /// pressure solver
    pressure_solver: laplacian_solver::LaplacianSolver,
    /// if has infinity smoke from source
    is_infinite_smoke: bool,
}
impl<C: canvas::Color + 'static> SingleSmokeGridScene<C> {
    
    /// Get the config table
    pub fn get_config_map() -> Map<String, Value> {
        let mut m = Map::new();
        m.insert(String::from("dt"), Value::Number(Number::from_f64(0.5_f64).unwrap()));
        m.insert(String::from("ds"), Value::Number(Number::from_f64(1.0_f64).unwrap()));
        m.insert(String::from("gx"), Value::Number(Number::from_f64(0_f64).unwrap()));
        m.insert(String::from("gy"), Value::Number(Number::from_f64(9.8_f64).unwrap()));
        m.insert(String::from("rho_air"), Value::Number(Number::from_f64(1_f64).unwrap()));
        m.insert(String::from("rho_smoke"), Value::Number(Number::from_f64(0.5_f64).unwrap()));
        m.insert(String::from("scene_width"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("scene_height"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_air_r"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_air_g"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_air_b"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_smoke_r"), Value::Number(Number::from(255_u8)));
        m.insert(String::from("c_smoke_g"), Value::Number(Number::from(0_u8)));
        m.insert(String::from("c_smoke_b"), Value::Number(Number::from(0_u8)));
        m.insert(String::from("init_d_x_min"), Value::Number(Number::from(16_usize)));
        m.insert(String::from("init_d_x_max"), Value::Number(Number::from(48_usize)));
        m.insert(String::from("init_d_y_min"), Value::Number(Number::from(2_usize)));
        m.insert(String::from("init_d_y_max"), Value::Number(Number::from(32_usize)));
        m.insert(String::from("is_infinite_smoke"), Value::Number(Number::from(0_u8)));
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
        let w = parser::get_from_parser_usize(parser, "scene_width");
        let h = parser::get_from_parser_usize(parser, "scene_height");
        let init_d_x_min = parser::get_from_parser_usize(parser, "init_d_x_min");
        let init_d_x_max = parser::get_from_parser_usize(parser, "init_d_x_max");
        let init_d_y_min = parser::get_from_parser_usize(parser, "init_d_y_min");
        let init_d_y_max = parser::get_from_parser_usize(parser, "init_d_y_max");
        let is_infinite_smoke = parser::get_from_parser_u8(parser, "is_infinite_smoke");
        Self {
            dt: dt,
            ds: ds,
            g: nd::Array1::<f64>::from_vec(vec![gx, gy]),
            rho_air: rho_air,
            rho_smoke: rho_smoke,
            c_air: C::new(c_air_r, c_air_g, c_air_b),
            c_smoke: C::new(c_smoke_r, c_smoke_g, c_smoke_b),
            sf_d: nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{
                if init_d_x_min <= i && i <= init_d_x_max && init_d_y_min <= j && j <= init_d_y_max { 1_f64 }
                else { 0_f64 }
            }),
            init_d_x_min: init_d_x_min,
            init_d_x_max: init_d_x_max,
            init_d_y_min: init_d_y_min,
            init_d_y_max: init_d_y_max,
            vf_vx: nd::Array2::<f64>::from_elem((w+1, h), 0_f64),
            vf_vy: nd::Array2::<f64>::from_elem((w, h+1), 0_f64),
            sf_p: nd::Array2::<f64>::zeros((w, h)),
            pressure_solver: laplacian_solver::LaplacianSolver::new(w, h),
            is_infinite_smoke: is_infinite_smoke != 0,
        }
    }

    fn _apply_buoyancy(&mut self) {
        let rho = self.sf_d.clone() * (self.rho_smoke - self.rho_air) + self.rho_air;

        let rho_x = to_staggered_x_grid(&rho, Box::new(boundary::NeumannBoundary));
        let fx = (self.rho_air - rho_x.clone()) / rho_x * self.g[0];
        
        let rho_y = to_staggered_y_grid(&rho, Box::new(boundary::NeumannBoundary));
        let fy = (self.rho_air - rho_y.clone()) / rho_y * self.g[1];
    
        self.vf_vx = self.vf_vx.clone() + fx * self.dt;
        self.vf_vy = self.vf_vy.clone() + fy * self.dt;
        self._velocity_boundary_adjust();
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
    
    fn _pressure_projection(&mut self) {
        let (w, h) = self.sf_d.dim();
        let sz = w * h;

        // pressure solver
        let divergence = self._solve_velocity_divergence();
        let b = divergence * self.ds * self.ds;
        let b = b.into_shape(sz).unwrap();
        self.sf_p = self.pressure_solver.solve(b).into_shape((w, h)).unwrap();
    
        // pressure projection
        let px = to_staggered_x_grid(&self.sf_p, Box::new(boundary::NeumannBoundary));
        let pxg = staggered_x_grid_gradient_as_staggered_x(&px, self.ds);
        self.vf_vx = self.vf_vx.clone() - pxg;
        
        let py = to_staggered_y_grid(&self.sf_p, Box::new(boundary::NeumannBoundary));
        let pyg = staggered_y_grid_gradient_as_staggered_y(&py, self.ds);
        self.vf_vy = self.vf_vy.clone() - pyg;

        self._velocity_boundary_adjust();
    }

    fn _velocity_advection(&mut self) {
        let (w, h) = self.sf_d.dim();
    
        let xs = nd::Array2::<f64>::from_shape_fn((w+1, h), |(i, _j)|{ i as f64 }) - self.vf_vx.clone() * self.dt / self.ds;
        let ys = nd::Array2::<f64>::from_shape_fn((w+1, h), |(_i, j)|{ j as f64 });
        self.vf_vx.assign(&sample_with_spatial_index(&self.vf_vx, xs, ys));
    
        let xs = nd::Array2::<f64>::from_shape_fn((w, h+1), |(i, _j)|{ i as f64 });
        let ys = nd::Array2::<f64>::from_shape_fn((w, h+1), |(_i, j)|{ j as f64 }) - self.vf_vy.clone() * self.dt / self.ds;
        self.vf_vy.assign(&sample_with_spatial_index(&self.vf_vy, xs, ys));

        self._velocity_boundary_adjust();
    }
    
    fn _density_advection(&mut self) {
        let (w, h) = self.sf_d.dim();
        
        let (vx_0, vx_1) = split_staggered_x_grid(&self.vf_vx);
        let vx_mid = (vx_0 + vx_1) * 0.5_f64;
        let xs = nd::Array2::<f64>::from_shape_fn((w, h), |(i, _j)|{ i as f64 }) - vx_mid * self.dt / self.ds;
        
        let (vy_0, vy_1) = split_staggered_y_grid(&self.vf_vy);
        let vy_mid = (vy_0 + vy_1) * 0.5_f64;
        let ys = nd::Array2::<f64>::from_shape_fn((w, h), |(_i, j)|{ j as f64 }) - vy_mid * self.dt / self.ds;

        self.sf_d = sample_with_spatial_index(&self.sf_d, xs, ys);

        if self.is_infinite_smoke {
            let sf_d_source = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{
                if self.init_d_x_min <= i && i <= self.init_d_x_max && self.init_d_y_min <= j && j <= self.init_d_y_max { 1_f64 }
                else { 0_f64 }
            });
            self.sf_d.zip_mut_with(&sf_d_source, |a, b| {
                *a = a.max(*b);
            });
        }
    }

    /// Do simulation in the temporal step
    pub fn sim(&mut self) {
        self._apply_buoyancy();
        self._pressure_projection();
        self._velocity_advection();
        self._density_advection();
    }

    /// Visualization of smoke density
    pub fn visualize_density(&self) -> Vec<Box<dyn render::Shape::<C>>> {
        let (w, h) = self.sf_d.dim();
        let mut shapes: Vec<Box<dyn render::Shape::<C>>> = Vec::new();
        for i in 0..w {
            for j in 0..h {
                let mut cs = self.c_smoke.clone();
                cs.set_opacity(self.sf_d[[i, j]]);
                let mut ca = self.c_air.clone();
                ca.set_opacity(1.0 - self.sf_d[[i, j]]);
                shapes.push(Box::new(render::RectShape::<C>::new(
                    cs.mix(ca),
                    (w-1-i) as f64 + 0.5,
                    (h-1-j) as f64 + 0.5,
                    1.0,
                    1.0,
                )));
            }
        }
        shapes
    }
}

/// Smoke simulation performed with Lagrangian method
///
/// Examples can be found in ```src/examples/```
pub struct SingleSmokeParticleScene<C: canvas::Color> {
    /// discrete of time
    dt: f64,
    /// discrete of space
    ds: f64,
    /// color of air for visualization
    c_air: C,
    /// color of smoke for visualization
    c_smoke: C,
    /// width, index coordinate
    w: usize,
    /// height, index coordinate
    h: usize,
    /// particle count per cell dimension
    p_cnt_per_cell_dimension: usize,
    /// particle total count
    p_cnt: usize,
    /// particle pos x, world coordinate
    p_px: nd::Array1::<f64>,
    /// particle pos y, world coordinate
    p_py: nd::Array1::<f64>,
    /// particle velocity x
    p_vx: nd::Array1::<f64>,
    /// particle velocity y
    p_vy: nd::Array1::<f64>,
    /// particle acceleration x
    p_ax: nd::Array1::<f64>,
    /// particle acceleration y
    p_ay: nd::Array1::<f64>,
    /// particle is smoke or not
    p_is_smoke: nd::Array1::<bool>,
    /// particle buoyancy force x, world coordinate
    p_buo_fx: nd::Array1::<f64>,
    /// particle buoyancy force y, world coordinate
    p_buo_fy: nd::Array1::<f64>,
    /// particle pressure
    p_pres: nd::Array1::<f64>,
    /// particle volume
    p_volume: f64,
    /// particle spatial index
    p_grid_index: nd::Array2::<Vec<usize>>,
}
impl<C: canvas::Color + 'static> SingleSmokeParticleScene<C> {
     
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
        m.insert(String::from("init_d_x_min"), Value::Number(Number::from(16_usize)));
        m.insert(String::from("init_d_x_max"), Value::Number(Number::from(48_usize)));
        m.insert(String::from("init_d_y_min"), Value::Number(Number::from(8_usize)));
        m.insert(String::from("init_d_y_max"), Value::Number(Number::from(16_usize)));
        m.insert(String::from("particle_cnt_per_cell"), Value::Number(Number::from(1_u8)));
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
        let init_d_x_min = parser::get_from_parser_usize(parser, "init_d_x_min");
        let init_d_x_max = parser::get_from_parser_usize(parser, "init_d_x_max");
        let init_d_y_min = parser::get_from_parser_usize(parser, "init_d_y_min");
        let init_d_y_max = parser::get_from_parser_usize(parser, "init_d_y_max");
        let particle_cnt_per_cell = parser::get_from_parser_usize(parser, "particle_cnt_per_cell"); // number of particle in one cell dimension
        let pds = ds / particle_cnt_per_cell as f64;
        let pds_half = pds * 0.5;
        let particle_cnt_total = particle_cnt_per_cell * particle_cnt_per_cell * w * h;
        let mut v_px: Vec<f64> = Vec::with_capacity(particle_cnt_total);
        let mut v_py: Vec<f64> = Vec::with_capacity(particle_cnt_total);
        let mut v_is_smoke: Vec<bool> = Vec::with_capacity(particle_cnt_total);
        let mut v_buo_fx: Vec<f64> = Vec::with_capacity(particle_cnt_total);
        let mut v_buo_fy: Vec<f64> = Vec::with_capacity(particle_cnt_total);
        let mut v_grid_index: Vec<Vec<usize>> = Vec::with_capacity(w * h);
        let mut cnt = 0;
        for i in 0..w {
            for j in 0..h {
                let is_smoke =  {
                    if init_d_x_min <= i && i <= init_d_x_max && init_d_y_min <= j && j <= init_d_y_max { true }
                    else { false }
                };
                let rho = {
                    if is_smoke { rho_smoke }
                    else { rho_air }
                };
                let rho_eff = (rho_air - rho) / rho;
                let mut ids: Vec<usize> = Vec::new();
                for ii in 0..particle_cnt_per_cell {
                    for jj in 0..particle_cnt_per_cell {
                        v_px.push(ds * i as f64 + pds_half + pds * ii as f64);
                        v_py.push(ds * j as f64 + pds_half + pds * jj as f64);
                        v_is_smoke.push(is_smoke);
                        v_buo_fx.push(rho_eff * gx);
                        v_buo_fy.push(rho_eff * gy);
                        ids.push(cnt);
                        cnt = cnt + 1;
                    }
                }
                v_grid_index.push(ids);
            }
        }
        Self {
            dt: dt,
            ds: ds,
            c_air: C::new(c_air_r, c_air_g, c_air_b),
            c_smoke: C::new(c_smoke_r, c_smoke_g, c_smoke_b),
            w: w,
            h: h,
            p_cnt_per_cell_dimension: particle_cnt_per_cell,
            p_cnt: particle_cnt_total,
            p_px: nd::Array1::<f64>::from_vec(v_px),
            p_py: nd::Array1::<f64>::from_vec(v_py),
            p_vx: nd::Array1::<f64>::zeros(particle_cnt_total),
            p_vy: nd::Array1::<f64>::zeros(particle_cnt_total),
            p_ax: nd::Array1::<f64>::zeros(particle_cnt_total),
            p_ay: nd::Array1::<f64>::zeros(particle_cnt_total),
            p_is_smoke: nd::Array1::<bool>::from_vec(v_is_smoke),
            p_buo_fx: nd::Array1::<f64>::from_vec(v_buo_fx),
            p_buo_fy: nd::Array1::<f64>::from_vec(v_buo_fy),
            p_pres: nd::Array1::<f64>::zeros(particle_cnt_total),
            p_volume: ds * ds / particle_cnt_per_cell as f64 / particle_cnt_per_cell as f64,
            p_grid_index: nd::Array2::<Vec<usize>>::from_shape_vec((w, h), v_grid_index).unwrap(),
        }
    }

    fn _update_grid_index(&mut self) {
        let mut v_grid_index: Vec<Vec<usize>> = vec![Vec::<usize>::new(); self.w * self.h];
        for pi in 0..self.p_cnt {
            let px = self.p_px[pi];
            let py = self.p_py[pi];
            let wi = (px / self.ds).floor() as usize;
            let hi = (py / self.ds).floor() as usize;
            v_grid_index[wi*self.h+hi].push(pi);
        }
        self.p_grid_index = nd::Array2::<Vec<usize>>::from_shape_vec((self.w, self.h), v_grid_index).unwrap();
    }

    fn _kernel(&self, dist: f64) -> f64 {
        if dist > self.ds + 0.0001_f64 { 0_f64 }
        else { 1_f64 / (3.1415926_f64 * self.ds.powi(2)) * (dist * dist / self.ds / self.ds).exp() }
    }

    fn _boundary_kernel(&self, dx: f64, dy: f64) -> f64 {
        let pds = self.ds / self.p_cnt_per_cell_dimension as f64;
        let pds_half = pds * 0.5;
        let ps_half = self.ds * 0.5;
        let mut res = 0_f64;
        for i in 0..self.p_cnt_per_cell_dimension {
            for j in 0..self.p_cnt_per_cell_dimension {
                let dpx = pds_half + pds * i as f64 - ps_half + dx;
                let dpy = pds_half + pds * j as f64 - ps_half + dy;
                let d = (dpx * dpx + dpy * dpy).sqrt();
                res += self._kernel(d);
            }
        }
        res
    }

    fn _pressure_projection(&mut self) {
        let mut p_pressure: Vec<f64> = vec![0_f64; self.p_cnt]; // particle pressure from spatial density
        let ds_half = self.ds * 0.5_f64;
        // traverse particles in every grid
        for x in 0..self.w {
            for y in 0..self.h {
                for pi in &self.p_grid_index[[x, y]] {
                    // println!("Pi: {}", pi);
                    let mut local_volume = self.p_volume * self._kernel(0_f64);
                    let p_px = self.p_px[*pi];
                    let p_py = self.p_py[*pi];
                    // traverse particles in the neighbor grids
                    for dx in 0..3 {
                        for dy in 0..3 {
                            // println!("dx {}, dy {}", dx, dy);
                            let bdpx = {
                                if x+dx < 1 { -1 }
                                else if x+dx >= self.w+1 { 1 }
                                else { 0 }
                            };
                            let bdpy = {
                                if y+dy < 1 { -1 }
                                else if y+dy >= self.h+1 { 1 }
                                else { 0 }
                            };
                            if bdpx == 0 && bdpy == 0 { // has the corresponding adjacent cell
                                for pj in &self.p_grid_index[[x+dx-1, y+dy-1]] {
                                    if pi == pj { continue }
                                    let dpx = p_px - self.p_px[*pj];
                                    let dpy = p_py - self.p_py[*pj];
                                    let d = (dpx * dpx + dpy * dpy).sqrt();
                                    if d < self.ds {
                                        local_volume += self.p_volume * self._kernel(d);
                                    }
                                }
                            }
                            else {
                                let bddpx = {
                                    if bdpx < 0 { p_px + ds_half }
                                    else if bdpx > 0 { p_px - (self.w as f64 * self.ds + ds_half) }
                                    else { p_px - (self.ds * (x + dx - 1) as f64 + ds_half) }
                                };
                                let bddpy = {
                                    if bdpy < 0 { p_py + ds_half }
                                    else if bdpy > 0 { p_py - (self.h as f64 * self.ds + ds_half) }
                                    else { p_py - (self.ds * (y + dy - 1) as f64 + ds_half) }
                                };
                                local_volume += self.p_volume * self._boundary_kernel(bddpx, bddpy);
                            }
                        }
                    }
                    // p_pressure[*pi] = (local_volume.powi(7) - 1_f64).max(0_f64) * 0.001_f64; // no attraction force
                    // p_pressure[*pi] = (local_volume - 1_f64); // no attraction force
                    p_pressure[*pi] = ((local_volume / 1.44572775).powi(1) - 1_f64).max(0_f64) * 0.01_f64; // no attraction force
                    // p_pressure[*pi] = local_volume; // no attraction force
                }
            }
        }
        self.p_pres = nd::Array1::<f64>::from_vec(p_pressure);
        // println!("{}", self.p_pres);
    }

    fn _apply_buoyancy(&mut self) {
        let mut v_p_ax: Vec<f64> = vec![0_f64; self.p_cnt];
        let mut v_p_ay: Vec<f64> = vec![0_f64; self.p_cnt];
        // traverse particles in every grid
        for x in 0..self.w {
            for y in 0..self.h {
                for pi in &self.p_grid_index[[x, y]] {
                    let mut ax = 0_f64;
                    let mut ay = 0_f64;
                    // traverse particles in the neighbor grids
                    for dx in 0..3 {
                        if x+dx < 1 { continue }
                        if x+dx >= self.w+1 { continue }
                        for dy in 0..3 {
                            if y+dy < 1 { continue }
                            if y+dy >= self.h+1 { continue }
                            for pj in &self.p_grid_index[[x+dx-1, y+dy-1]] {
                                if pi == pj { continue }
                                let dpx = self.p_px[*pi] - self.p_px[*pj];
                                let dpy = self.p_py[*pi] - self.p_py[*pj];
                                let d = (dpx * dpx + dpy * dpy).sqrt();
                                if d < self.ds {
                                    let pp = self.p_pres[*pi] + self.p_pres[*pj] * self._kernel(d) / d;
                                    ax += dpx * pp;
                                    ay += dpy * pp;
                                }
                            }
                        }
                    }
                    v_p_ax[*pi] += self.p_buo_fx[*pi] - ax;
                    v_p_ay[*pi] += self.p_buo_fy[*pi] - ay;
                }
            }
        }
        self.p_ax = nd::Array1::<f64>::from_vec(v_p_ax);
        self.p_ay = nd::Array1::<f64>::from_vec(v_p_ay);
    }

    fn _advection(&mut self) {
        for pi in 0..self.p_cnt {
            self.p_vx[pi] += self.p_ax[pi] * self.dt;
            self.p_vy[pi] += self.p_ay[pi] * self.dt;
            self.p_px[pi] += self.p_vx[pi] * self.dt;
            self.p_py[pi] += self.p_vy[pi] * self.dt;
            if self.p_px[pi] < 0_f64 {
                self.p_px[pi] = self.ds * 0.1_f64;
                self.p_vx[pi] = 0_f64;
                // self.p_px[pi] = -self.p_px[pi];
                // self.p_vx[pi] = -self.p_vx[pi];
            }
            if self.p_py[pi] < 0_f64 {
                self.p_py[pi] = self.ds * 0.1_f64;
                self.p_vy[pi] = 0_f64;
                // self.p_py[pi] = -self.p_py[pi];
                // self.p_vy[pi] = -self.p_vy[pi];
            }
            if self.p_px[pi] > self.w as f64 {
                self.p_px[pi] = self.w as f64 - self.ds * 0.1_f64;
                self.p_vx[pi] = 0_f64;
                // self.p_px[pi] = 2_f64 * self.w as f64 - self.p_px[pi];
                // self.p_vx[pi] = -self.p_vx[pi];
            }
            if self.p_py[pi] > self.h as f64 {
                self.p_py[pi] = self.h as f64 - self.ds * 0.1_f64;
                self.p_vy[pi] = 0_f64;
                // self.p_py[pi] = 2_f64 * self.h as f64 - self.p_py[pi];
                // self.p_vy[pi] = -self.p_vy[pi];
            }
        }
    }

    // fn _clamp(d: &mut nd::Array1::<f64>, low: f64, high: f64) {
    //     let d_new = d.map(|o| {
    //         if *o < low { 2.0 * low - *o}
    //         else if *o > high { 2.0 * high - *o }
    //         else { *o }
    //     });
    //     d.assign(&d_new);
    // }

    // fn _collide(&mut self) {
    //     let (w, h) = self.p_grid_index.dim();
    //     self._clamp(&mut self.p_px, 0.0, w as f64);
    //     self._clamp(&mut self.p_py, 0.0, h as f64);
    // }

    // fn _position_advection(&mut self) {
    //     self.p_px = self.p_px + self.p_vx.clone() * self.dt / self.ds;
    //     self.p_py = self.p_py + self.p_vy.clone() * self.dt / self.ds;
    //     self._collide()
    // }

    /// Do simulation in the temporal step
    pub fn sim(&mut self) {
        self._update_grid_index();
        self._pressure_projection();
        self._apply_buoyancy();
        self._advection();
    }

    /// Visualization of smoke density
    pub fn visualize_density(&self) -> Vec<C> {
        let (w, h) = self.p_grid_index.dim();
        let mut buffer = Vec::with_capacity(w* h);
        for j in 0..h {
            for i in 0..w {
                let mut smoke_particle_cnt = 0;
                let mut total_particle_cnt = 0;
                for pi in &self.p_grid_index[[i, j]] {
                    if self.p_is_smoke[*pi] {
                        smoke_particle_cnt += 1;
                    }
                    total_particle_cnt += 1;
                }
                let mut smoke_portion = 0_f64;
                if total_particle_cnt > 0 {
                    smoke_portion = smoke_particle_cnt as f64 / total_particle_cnt as f64;
                }
                let mut cs = self.c_smoke.clone();
                cs.set_opacity(smoke_portion);
                let mut ca = self.c_air.clone();
                ca.set_opacity(1.0 - smoke_portion);
                buffer.push(cs.mix(ca));
            }
        }
        buffer.reverse();
        buffer
    }
}
