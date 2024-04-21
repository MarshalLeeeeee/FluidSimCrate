use ndarray as nd;
use sprs::{CsMat, CsVec, TriMat};
use sprs_ldl::*;
use fluid_sim::parser;
use fluid_sim::canvas;

fn main() {
    let raw_parser_maps = vec!(
        canvas::Canvas::get_config_map(),
    );
    let parser_map = parser::register(raw_parser_maps);
    let parser = parser::parse(parser_map);

    let mut canvas = canvas::Canvas::new_by_parser(&parser);
    let width = canvas.get_width();
    let height = canvas.get_height();

    let dt = 0.01_f64;
    let g = nd::Array1::<f64>::from_vec(vec![0_f64, 1_f64]);
    let mut d = nd::Array2::<f64>::from_shape_fn((width, height), |(i, j)|{
        if 5 <= i && i <= 10 && 2 <= j && j <= 4 { 1_f64 }
        else { 0_f64 }
    });
    let mut vx = nd::Array2::<f64>::from_elem((width+1, height), 0_f64);
    let mut vy = nd::Array2::<f64>::from_elem((width, height+1), 0_f64);
    let c_air = canvas::RGBAColor::new(0, 0, 0, 255);
    let c_fluid = canvas::RGBAColor::new(0, 255, 0, 255);
    let rho_air = 1_f64;
    let rho_fluid = 0.5_f64;

    while canvas.is_valid() {
        let rho = nd_solve_rho(&d, rho_air, rho_fluid);
        let (fx, fy) = nd_solve_force(&rho, &g, rho_air);
        nd_velocity_advection(&mut vx, &mut vy, fx, fy, dt);
        nd_fix_velocit_boundary(&mut vx, &mut vy);
        nd_pressure_projection(&mut vx, &mut vy, &rho);
        nd_fix_velocit_boundary(&mut vx, &mut vy);
        nd_density_advection(&mut d, &vx, &vy, dt);
        let buffer = na_to_buffer(&d, &c_air, &c_fluid);
        canvas.refresh(&buffer);
    }
}

fn jagged_x_split(d: &nd::Array2::<f64>) -> (nd::Array2::<f64>, nd::Array2::<f64>) {
    let d0 = d.slice(nd::s![..-1, ..]);
    let d1 = d.slice(nd::s![1.., ..]);
    (d0.to_owned(), d1.to_owned())
}

fn jagged_x_gradient(d: &nd::Array2::<f64>) -> nd::Array2::<f64> {
    let (d0, d1) = jagged_x_split(d);
    d1 - d0
}

fn jagged_x_gradient_jagged(d: &nd::Array2::<f64>) -> nd::Array2::<f64> {
    let (d0, d1) = jagged_x_split(d);
    let dg = d1 - d0;
    let dg_x0 = dg.slice(nd::s![0..1, ..]).to_owned();
    let mut dgg = nd::Array2::zeros(d.dim());
    dgg.slice_mut(nd::s![0..1, ..]).assign(&dg_x0);
    dgg.slice_mut(nd::s![1.., ..]).assign(&dg);
    dgg
}

fn to_jagged_x(d: &nd::Array2::<f64>) -> nd::Array2::<f64> {
    let (w, h) = d.dim();
    let (d0, d1) = jagged_x_split(d);
    let d_mid = (d0 + d1) * 0.5;
    let mut d_new = nd::Array2::<f64>::from_elem((w+1,h), 0_f64);
    d_new.slice_mut(nd::s![..1, ..]).assign(&d.slice(nd::s![..1, ..]).to_owned());
    d_new.slice_mut(nd::s![1..-1, ..]).assign(&d_mid);
    d_new.slice_mut(nd::s![-1.., ..]).assign(&d.slice(nd::s![-1.., ..]).to_owned());
    d_new
}

fn jagged_y_split(d: &nd::Array2::<f64>) -> (nd::Array2::<f64>, nd::Array2::<f64>) {
    let d0 = d.slice(nd::s![.., ..-1]);
    let d1 = d.slice(nd::s![.., 1..]);
    (d0.to_owned(), d1.to_owned())
}

fn jagged_y_gradient(d: &nd::Array2::<f64>) -> nd::Array2::<f64> {
    let (d0, d1) = jagged_y_split(d);
    d1 - d0
}

fn jagged_y_gradient_jagged(d: &nd::Array2::<f64>) -> nd::Array2::<f64> {
    let (d0, d1) = jagged_y_split(d);
    let dg = d1 - d0;
    let dg_y0 = dg.slice(nd::s![.., ..1]).to_owned();
    let mut dgg = nd::Array2::zeros(d.dim());
    dgg.slice_mut(nd::s![.., ..1]).assign(&dg_y0);
    dgg.slice_mut(nd::s![.., 1..]).assign(&dg);
    dgg
}

fn to_jagged_y(d: &nd::Array2::<f64>) -> nd::Array2::<f64> {
    let (w, h) = d.dim();
    let (d0, d1) = jagged_y_split(d);
    let d_mid = (d0 + d1) * 0.5;
    let mut d_new = nd::Array2::<f64>::from_elem((w,h+1), 0_f64);
    d_new.slice_mut(nd::s![.., ..1]).assign(&d.slice(nd::s![.., ..1]).to_owned());
    d_new.slice_mut(nd::s![.., 1..-1]).assign(&d_mid);
    d_new.slice_mut(nd::s![.., -1..]).assign(&d.slice(nd::s![.., -1..]).to_owned());
    d_new
}

fn nd_solve_rho(d: &nd::Array2::<f64>, rho_air: f64, rho_fluid: f64) -> nd::Array2::<f64> {
    d * (rho_fluid - rho_air) + rho_air
}

fn nd_solve_force(rho: &nd::Array2::<f64>, g: &nd::Array1::<f64>, rho_air: f64) -> (nd::Array2::<f64>, nd::Array2::<f64>) {
    let rho_x = to_jagged_x(&rho);
    let fx = (rho_air - rho_x.clone()) / rho_x * g[0];
    
    let rho_y = to_jagged_y(&rho);
    let fy = (rho_air - rho_y.clone()) / rho_y * g[1];

    (fx, fy)
}

fn nd_velocity_advection(vx: &mut nd::Array2::<f64>, vy: &mut nd::Array2::<f64>, force_x: nd::Array2::<f64>, force_y: nd::Array2::<f64>, dt: f64) {
    let vx_g = jagged_x_gradient_jagged(&vx);
    let vx_pre = vx.clone() + force_x * dt - vx_g * dt;
    let vx_pre_g = jagged_x_gradient_jagged(&vx_pre);
    let vx_new = (vx.clone() + vx_pre) * 0.5 - vx_pre_g * 0.5 * dt;
    vx.assign(&vx_new);
    
    let vy_g = jagged_y_gradient_jagged(&vy);
    let vy_pre = vy.clone() + force_y * dt - vy_g * dt;
    let vy_pre_g = jagged_y_gradient_jagged(&vy_pre);
    let vy_new = (vy.clone() + vy_pre) * 0.5 - vy_pre_g * 0.5 * dt;
    vy.assign(&vy_new);
}

fn nd_fix_velocit_boundary(vx: &mut nd::Array2::<f64>, vy: &mut nd::Array2::<f64>) {
    let (w, h) = vx.dim();
    let w = w - 1;
    for j in 0..h {
        if vx[[0, j]] < 0_f64 { vx[[0, j]] = 0_f64; }
        if vx[[w-1, j]] > 0_f64 { vx[[w-1, j]] = 0_f64; }
    }
    for i in 0..w {
        if vy[[i, 0]] < 0_f64 { vy[[i, 0]] = 0_f64; }
        if vy[[i, h-1]] > 0_f64 { vy[[i, h-1]] = 0_f64; }
    }
}

fn nd_solve_velocity_divengence(vx: &nd::Array2::<f64>, vy: &nd::Array2::<f64>) -> nd::Array2::<f64> {
    jagged_x_gradient(&vx) + jagged_y_gradient(&vy)
}

fn nd_pressure_projection(vx: &mut nd::Array2::<f64>, vy: &mut nd::Array2::<f64>, rho: &nd::Array2::<f64>) {
    let (w, h) = rho.dim();
    let sz = w * h;

    let divergence = nd_solve_velocity_divengence(&vx, &vy);
    let mut pa = TriMat::<f32>::new((sz, sz));
    let mut pbi = Vec::with_capacity(sz);
    let mut pbv: Vec<f32> = Vec::with_capacity(sz);
    for i in 0..w {
        for j in 0..h {
            let mut cnt = 0;
            if i > 0 {
                pa.add_triplet(i*h+j, (i-1)*h+j, -1_f32);
                cnt += 1;
            }
            if i < w-1 {
                pa.add_triplet(i*h+j, (i+1)*h+j, -1_f32);
                cnt += 1;
            }
            if j > 0 {
                pa.add_triplet(i*h+j, i*h+j-1, -1_f32);
                cnt += 1;
            }
            if j < h-1 {
                pa.add_triplet(i*h+j, i*h+j+1, -1_f32);
                cnt += 1;
            }
            pa.add_triplet(i*h+j, i*h+j, cnt as f32);

            pbi.push(i*h+j);
            pbv.push((divergence[[i, j]] * rho[[i, j]]) as f32);
        }
    }
    let pa: CsMat<_> = pa.to_csr();
    let pb = CsVec::new(sz, pbi, pbv);
    let ldl = Ldl::default();
    let system = ldl.numeric(pa.view()).unwrap();
    let p = system.solve(pb.to_dense());
    let p = nd::Array2::<f64>::from_shape_fn((w, h), |(i,j)|{
        p[i*h+j] as f64
    });

    // pressure projection
    let px = to_jagged_x(&p);
    let pxg = jagged_x_gradient_jagged(&px);
    *vx = vx.clone() + pxg / to_jagged_x(&rho);
    
    let py = to_jagged_y(&p);
    let pyg = jagged_y_gradient_jagged(&py);
    *vy = vy.clone() + pyg / to_jagged_y(&rho);
}

fn nd_density_advection(d: &mut nd::Array2::<f64>, vx: &nd::Array2::<f64>, vy: &nd::Array2::<f64>, dt: f64) {
    let (w, h) = d.dim();
    let xs = nd::Array2::<f64>::from_shape_fn((w, h), |(i, _j)|{ i as f64 });
    let ys = nd::Array2::<f64>::from_shape_fn((w, h), |(_i, j)|{ j as f64 });

    let (vx_0, vx_1) = jagged_x_split(&vx);
    let vx_mid = (vx_0 + vx_1) * 0.5_f64;
    let dx = vx_mid * dt;
    
    let (vy_0, vy_1) = jagged_y_split(&vy);
    let vy_mid = (vy_0 + vy_1) * 0.5_f64;
    let dy = vy_mid * dt;
    
    let dxs = (xs - dx).map(|o|{ o.clamp(0_f64, w as f64 - 1_f64) });
    let dys = (ys - dy).map(|o|{ o.clamp(0_f64, h as f64 - 1_f64) });
    
    let dxs0 = dxs.map(|o|{ o.floor() });
    let dxs1 = dxs.map(|o|{ o.ceil() });
    let dys0 = dys.map(|o|{ o.floor() });
    let dys1 = dys.map(|o|{ o.ceil() });
    
    let dx0y0 = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{ d[[dxs0[[i,j]] as usize, dys0[[i,j]] as usize]] });
    let dx0y1 = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{ d[[dxs0[[i,j]] as usize, dys1[[i,j]] as usize]] });
    let dx1y0 = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{ d[[dxs1[[i,j]] as usize, dys0[[i,j]] as usize]] });
    let dx1y1 = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{ d[[dxs1[[i,j]] as usize, dys1[[i,j]] as usize]] });
    
    let ddxs = dxs - dxs0;
    let ddys = dys - dys0;
    
    let dxy0 = (dx1y0 - dx0y0.clone()) * ddxs.clone() + dx0y0;
    let dxy1 = (dx1y1 - dx0y1.clone()) * ddxs + dx0y1;
    let dxy = (dxy1 - dxy0.clone()) * ddys + dxy0;

    d.assign(&dxy);
}

fn na_to_buffer(d: &nd::Array2::<f64>, c_air: &canvas::RGBAColor, c_fluid: &canvas::RGBAColor) -> Vec<canvas::RGBAColor> {
    let (width, height) = d.dim();
    let mut buffer = Vec::with_capacity(width * height);
    for j in 0..height {
        for i in 0..width {
            buffer.push(canvas::RGBAColor::mix(&c_fluid, &c_air, d[[i, j]]));
        }
    }
    buffer.reverse();
    buffer
}
