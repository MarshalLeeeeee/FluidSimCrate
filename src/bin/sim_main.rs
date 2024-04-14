use nalgebra as na;
use std::time;
use std::thread;
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

    let dt = 0.1_f64;
    let gravity: na::Vector2<f64> = na::Vector2::new(0_f64, 10_f64);
    let rho_air = 1.0;
    let color_air = canvas::RGBAColor::new(0, 0, 0, 255);
    let rho_smoke = 0.8;
    let color_smoke = canvas::RGBAColor::new(0, 255, 0, 255);
    let mut density = na::DMatrix::<f64>::from_fn(
        width,
        height,
        |i, j| { 
            if i > 40 && i < 60 && j > 40 && j < 60 { 1_f64 }
            else { 0_f64 }
        }
    );
    let mut velocity_x = na::DMatrix::<f64>::from_fn(width+1, height, |i, j| { 0_f64 });
    let mut velocity_y = na::DMatrix::<f64>::from_fn(width, height+1, |i, j| { 0_f64 });

    while canvas.is_valid() {
        (velocity_x, velocity_y) = advection(&density, rho_air, rho_smoke, velocity_x, velocity_y, &gravity, dt);
        density = update_density(density, &velocity_x, &velocity_y, dt);
        let buffer = update_buffer(&density, &color_air, &color_smoke);
        canvas.refresh(&buffer);
    }
}

fn velocity_x_split(velocity_x: &na::DMatrix::<f64>) -> (na::DMatrix::<f64>, na::DMatrix::<f64>) {
    let (velocity_x_width, velocity_x_height) = velocity_x.shape();
    let velocity_x_h = velocity_x.view((1,0), (velocity_x_width-1, velocity_x_height));
    let velocity_x_l = velocity_x.view((0,0), (velocity_x_width-1, velocity_x_height));
    (velocity_x_h.into(), velocity_x_l.into())
}

fn velocity_x_gradient(velocity_x: &na::DMatrix::<f64>) -> na::DMatrix::<f64> {
    let (velocity_x_h, velocity_x_l) = velocity_x_split(&velocity_x);
    let mut velocity_x_g = (velocity_x_h - velocity_x_l).insert_row(0, 0_f64);
    let velocity_x_g_r1 = velocity_x_g.row(1).into_owned();
    velocity_x_g.set_row(0, &velocity_x_g_r1);
    velocity_x_g
}

fn velocity_y_split(velocity_y: &na::DMatrix::<f64>) -> (na::DMatrix::<f64>, na::DMatrix::<f64>) {
    let (velocity_y_width, velocity_y_height) = velocity_y.shape();
    let velocity_y_h = velocity_y.view((0,1), (velocity_y_width, velocity_y_height-1));
    let velocity_y_l = velocity_y.view((0,0), (velocity_y_width, velocity_y_height-1));
    (velocity_y_h.into(), velocity_y_l.into())
}

fn velocity_y_gradient(velocity_y: &na::DMatrix::<f64>) -> na::DMatrix::<f64> {
    let (velocity_y_h, velocity_y_l) = velocity_y_split(&velocity_y);
    let mut velocity_y_g = (velocity_y_h - velocity_y_l).insert_column(0, 0_f64);
    let velocity_y_g_c1 = velocity_y_g.column(1).into_owned();
    velocity_y_g.set_column(0, &velocity_y_g_c1);
    velocity_y_g
}

/// convert (w,h) to (w+1, h)
fn to_edged_x(data: &na::DMatrix::<f64>) -> na::DMatrix::<f64> {
    let (data_width, data_height) = data.shape();
    let data_h = data.view((1,0), (data_width-1, data_height));
    let data_l = data.view((0,0), (data_width-1, data_height));
    let mut data_mid = (data_h + data_l) * 0.5;
    data_mid = data_mid.insert_row(0, 0_f64);
    data_mid.set_row(0, &data.row(0).into_owned());
    data_mid = data_mid.insert_row(data_width, 0_f64);
    data_mid.set_row(data_width, &data.row(data_width-1).into_owned());
    data_mid
}

/// convert (w,h) to (w, h+1)
fn to_edged_y(data: &na::DMatrix::<f64>) -> na::DMatrix::<f64> {
    let (data_width, data_height) = data.shape();
    let data_h = data.view((0,1), (data_width, data_height-1));
    let data_l = data.view((0,0), (data_width, data_height-1));
    let mut data_mid = (data_h + data_l) * 0.5;
    data_mid = data_mid.insert_column(0, 0_f64);
    data_mid.set_column(0, &data.column(0).into_owned());
    data_mid = data_mid.insert_column(data_height, 0_f64);
    data_mid.set_column(data_height, &data.column(data_height-1).into_owned());
    data_mid
}

fn advection(density: &na::DMatrix::<f64>, rho_air: f64, rho_smoke: f64, velocity_x: na::DMatrix::<f64>, velocity_y: na::DMatrix::<f64>, gravity: &na::Vector2<f64>, dt: f64) -> (na::DMatrix::<f64>, na::DMatrix::<f64>) {
    let density_x = to_edged_x(density);
    let density_smoke = density_x.clone() * rho_smoke;
    let density_air = (density_x.clone()*(-1_f64)).add_scalar(1_f64) * rho_air;
    let density_mix = density_smoke + density_air;
    let force_x = density_mix.add_scalar(-rho_air).component_div(&density_mix) * gravity.y;
    
    let (vx_w, vx_h) = velocity_x.shape();
    let velocity_x_g = velocity_x_gradient(&velocity_x);
    let mut velocity_x_p = na::DMatrix::<f64>::zeros(vx_w, vx_h);
    let euler_advection_x = force_x * dt - velocity_x_g * dt;
    velocity_x.add_to(&euler_advection_x, &mut velocity_x_p);
    let velocity_x_p_g = velocity_x_gradient(&velocity_x_p);
    velocity_x_p = (velocity_x + velocity_x_p) * 0.5 - velocity_x_p_g * 0.5 * dt;
    
    let density_y = to_edged_y(density);
    let density_smoke = density_y.clone() * rho_smoke;
    let density_air = (density_y.clone()*(-1_f64)).add_scalar(1_f64) * rho_air;
    let density_mix = density_smoke + density_air;
    let force_y = density_mix.add_scalar(-rho_air).component_div(&density_mix) * gravity.x;
    
    let (vy_w, vy_h) = velocity_y.shape();
    let velocity_y_g = velocity_y_gradient(&velocity_y);
    let mut velocity_y_p = na::DMatrix::<f64>::zeros(vy_w, vy_h);
    let euler_advection_y = force_y * dt - velocity_y_g * dt;
    velocity_y.add_to(&euler_advection_y, &mut velocity_y_p);
    let velocity_y_p_g = velocity_y_gradient(&velocity_y_p);
    velocity_y_p = (velocity_y + velocity_y_p) * 0.5 - velocity_y_p_g * 0.5 * dt;

    (velocity_x_p, velocity_y_p)
}

fn pressure_solution(density: &na::DMatrix::<f64>, velocity: &na::DMatrix::<na::Vector2<f64>>, force: &na::Vector2<f64>) {
    // TODO: A linear system solver
    // let sz = 1000;
    // let r = na::DMatrix::<f64>::from_fn(sz, sz, |i, j| {
    //     if rand::random::<f64>() > 0.9_f64 { rand::random::<f64>() }
    //     else { 0_f64 }
    // });
    // let mut b = na::DMatrix::<f64>::zeros(sz, 1);
    // b = b.add_scalar(1_f64);

    // let lu = na::LU::new(r);
    // if let Some(x) = lu.solve(&b) {
    //     println!("ok!");
    //     // println!("{}", x);
    // }

    // let qr = na::QR::new(r);
    // if let Some(x) = qr.solve(&b) {
    //     println!("ok!");
    //     // println!("{}", x);
    // }

    // let fplu = na::FullPivLU::new(r);
    // if let Some(x) = fplu.solve(&b) {
    //     println!("ok!");
    //     // println!("{}", x);
    // }


    // use sprs::{CsMat, TriMat, SpSolver};
    // use sprs::SprsShape;
    // // 创建一个高维稀疏矩阵和一个向量b
    // let matrix_a = CsMat::new((3, 3), vec![0, 2, 3], vec![0, 1, 2], vec![4.0, 5.0, 6.0]);
    // let vector_b = vec![1.0, 2.0, 3.0];

    // // 解线性系统 Ax=b
    // let solver = SpSolver::factorize(&matrix_a).unwrap();
    // let solution = solver.solve(&vector_b).unwrap();
    // println!("Solution: {:?}", solution);
}

fn projection(density: &na::DMatrix::<f64>, velocity: &na::DMatrix::<na::Vector2<f64>>, force: &na::Vector2<f64>) {
    // TODO: apply projection
    // TODO: process boundary velocity
}

fn update_density(density: na::DMatrix::<f64>, velocity_x: &na::DMatrix::<f64>, velocity_y: &na::DMatrix::<f64>, dt: f64) -> na::DMatrix::<f64> {
    let (velocity_x_h, velocity_x_l) = velocity_x_split(&velocity_x);
    let velocity_x_mid = (velocity_x_h + velocity_x_l) * 0.5;
    let dx = velocity_x_mid * dt;
    
    let (velocity_y_h, velocity_y_l) = velocity_y_split(&velocity_y);
    let velocity_y_mid = (velocity_y_h + velocity_y_l) * 0.5;
    let dy = velocity_y_mid * dt;

    let (density_width, density_height) = density.shape();
    let mut new_density = na::DMatrix::zeros(density_width, density_height);
    for i in 0..density_width {
        for j in 0..density_height {
            let index = i*density_width+j;

            let ii = i as f64 - dx.get(index).unwrap();
            let jj = j as f64 - dy.get(index).unwrap();
            let ii0 = ii.floor() as usize;
            let ii1 = ii.ceil() as usize;
            let jj0 = jj.floor() as usize;
            let jj1 = jj.ceil() as usize;

            let d00 = get_density_by_index(&density, ii0, jj0);
            let d01 = get_density_by_index(&density, ii0, jj1);
            let d10 = get_density_by_index(&density, ii1, jj0);
            let d11 = get_density_by_index(&density, ii1, jj1);
            let d0 = d00 + (d01 - d00) * (ii - ii0 as f64);
            let d1 = d10 + (d11 - d10) * (ii - ii0 as f64);
            let d = d0 + (d1 - d0) * (jj - jj0 as f64);

            if let Some(v) = new_density.get_mut(index) {
                *v = d.max(0_f64).min(1_f64);
            }
        }
    }

    new_density
}

fn get_density_by_index(density: &na::DMatrix::<f64>, i: usize, j: usize) -> f64 {
    let (density_width, density_height) = density.shape();
    if i >= density_width { 0_f64 }
    else if j >= density_height { 0_f64 }
    else {
        match density.get(i * density_width + j) {
            Some(v) => *v,
            None => 0_f64,
        }
    }
}

fn update_buffer(density: &na::DMatrix::<f64>, color_air: &canvas::RGBAColor, color_smoke: &canvas::RGBAColor) -> Vec<canvas::RGBAColor> {
    let (density_width, density_height) = density.shape();
    let mut buffer = Vec::with_capacity(density_width * density_height);
    for d in density.iter() {
        buffer.push(color_smoke.mix(*d, color_air));
    }
    buffer.reverse();
    buffer
}
