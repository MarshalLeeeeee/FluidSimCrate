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

pub fn sample_with_spatial_index(d: &nd::Array2::<f64>, xs: nd::Array2::<f64>, ys: nd::Array2::<f64>) -> nd::Array2::<f64> {
    let (w, h) = d.dim();
    let xs = xs.map(|o|{ o.clamp(0_f64, w as f64 - 1_f64) });
    let ys = ys.map(|o|{ o.clamp(0_f64, h as f64 - 1_f64) });
    
    let xs0 = xs.map(|o|{ o.floor() });
    let xs1 = xs.map(|o|{ o.ceil() });
    let ys0 = ys.map(|o|{ o.floor() });
    let ys1 = ys.map(|o|{ o.ceil() });
    
    let dx0y0 = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{ d[[xs0[[i,j]] as usize, ys0[[i,j]] as usize]] });
    let dx0y1 = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{ d[[xs0[[i,j]] as usize, ys1[[i,j]] as usize]] });
    let dx1y0 = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{ d[[xs1[[i,j]] as usize, ys0[[i,j]] as usize]] });
    let dx1y1 = nd::Array2::<f64>::from_shape_fn((w, h), |(i, j)|{ d[[xs1[[i,j]] as usize, ys1[[i,j]] as usize]] });
    
    let kxs0 = xs - xs0;
    let kys0 = ys - ys0;
    
    let dxy0 = (dx1y0 - dx0y0.clone()) * kxs0.clone() + dx0y0;
    let dxy1 = (dx1y1 - dx0y1.clone()) * kxs0 + dx0y1;
    let dxy = (dxy1 - dxy0.clone()) * kys0 + dxy0;
    dxy
}
