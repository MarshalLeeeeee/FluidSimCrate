/// Module for grid (Euler) related functions
use ndarray as nd;
use crate::boundary;

pub fn split_staggered_x_grid(d: &nd::Array2::<f64>) -> (nd::Array2::<f64>, nd::Array2::<f64>) {
    let d0 = d.slice(nd::s![..-1, ..]);
    let d1 = d.slice(nd::s![1.., ..]);
    (d0.to_owned(), d1.to_owned())
}

pub fn staggered_x_grid_gradient(d: &nd::Array2::<f64>, ds: f64) -> nd::Array2::<f64> {
    let (d0, d1) = split_staggered_x_grid(d);
    (d1 - d0) / ds
}

pub fn staggered_x_grid_gradient_as_staggered_x(d: &nd::Array2::<f64>, ds: f64) -> nd::Array2::<f64> {
    let g = staggered_x_grid_gradient(&d, ds);
    to_staggered_x_grid(&g, Box::new(boundary::DirichletBoundary(0_f64)))
}

pub fn to_staggered_x_grid(d: &nd::Array2::<f64>, boundary: Box<dyn boundary::Boundary>) -> nd::Array2::<f64> {
    let (w, h) = d.dim();
    let (d0, d1) = split_staggered_x_grid(d);
    let d_mid = (d0 + d1) * 0.5;
    let mut d_new = nd::Array2::<f64>::from_elem((w+1,h), 0_f64);
    d_new.slice_mut(nd::s![..1, ..]).assign(&boundary.boundary(&d.row(0).to_owned()).insert_axis(nd::Axis(0)));
    d_new.slice_mut(nd::s![1..-1, ..]).assign(&d_mid);
    d_new.slice_mut(nd::s![-1.., ..]).assign(&boundary.boundary(&d.row(w-1).to_owned()).insert_axis(nd::Axis(0)));
    d_new
}

pub fn split_staggered_y_grid(d: &nd::Array2::<f64>) -> (nd::Array2::<f64>, nd::Array2::<f64>) {
    let d0 = d.slice(nd::s![.., ..-1]);
    let d1 = d.slice(nd::s![.., 1..]);
    (d0.to_owned(), d1.to_owned())
}

pub fn staggered_y_grid_gradient(d: &nd::Array2::<f64>, ds: f64) -> nd::Array2::<f64> {
    let (d0, d1) = split_staggered_y_grid(d);
    (d1 - d0) / ds
}

pub fn staggered_y_grid_gradient_as_staggered_y(d: &nd::Array2::<f64>, ds: f64) -> nd::Array2::<f64> {
    let g = staggered_y_grid_gradient(&d, ds);
    to_staggered_y_grid(&g, Box::new(boundary::DirichletBoundary(0_f64)))
}

pub fn to_staggered_y_grid(d: &nd::Array2::<f64>, boundary: Box<dyn boundary::Boundary>) -> nd::Array2::<f64> {
    let (w, h) = d.dim();
    let (d0, d1) = split_staggered_y_grid(d);
    let d_mid = (d0 + d1) * 0.5;
    let mut d_new = nd::Array2::<f64>::from_elem((w,h+1), 0_f64);
    d_new.slice_mut(nd::s![.., ..1]).assign(&boundary.boundary(&d.column(0).to_owned()).insert_axis(nd::Axis(1)));
    d_new.slice_mut(nd::s![.., 1..-1]).assign(&d_mid);
    d_new.slice_mut(nd::s![.., -1..]).assign(&boundary.boundary(&d.column(h-1).to_owned()).insert_axis(nd::Axis(1)));
    d_new
}

pub fn sample_with_offset(d: &nd::Array2::<f64>, dx: nd::Array2::<f64>, dy: nd::Array2::<f64>) -> nd::Array2::<f64> {
    let (w, h) = d.dim();
    let xs = nd::Array2::<f64>::from_shape_fn((w, h), |(i, _j)|{ i as f64 });
    let ys = nd::Array2::<f64>::from_shape_fn((w, h), |(_i, j)|{ j as f64 });
    
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
    dxy
}
