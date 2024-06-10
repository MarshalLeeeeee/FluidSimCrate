use std::cmp::{min, max};
use ndarray as nd;
use crate::canvas;

struct Vec2<T>(T, T);
impl<T> Vec2<T> {
    fn new(x: T, y: T) -> Self {
        Self(x, y)
    }
}

type UsizeVec2 = Vec2<usize>;
type F64Vec2 = Vec2<f64>;
impl F64Vec2 {
    fn ref_mul(&self, rhs: usize) -> F64Vec2 {
        F64Vec2::new(
            self.0 * rhs as f64,
            self.1 * rhs as f64,
        )
    }
}

pub struct Entry<C: canvas::Color> {
    color: C,
    pos: UsizeVec2, // pos in canvas coordinate
}

pub trait Shape<C: canvas::Color> {
    fn rasterize(&self, width: usize, height: usize, ratio: usize) -> Vec<Entry<C>>;
}

pub struct RectShape<C: canvas::Color> {
    color: C,
    pos: F64Vec2, // pos in scene coordinate
    size: F64Vec2, // size in scene coordinate
}
impl<C: canvas::Color> Shape<C> for RectShape<C> {
    fn rasterize(&self, width: usize, height: usize, ratio: usize) -> Vec<Entry<C>> {
        let canvas_pos_x = self.pos.0 * ratio as f64;
        let canvas_pos_y = self.pos.1 * ratio as f64;
        let canvas_size_w = self.size.0 * ratio as f64;
        let canvas_size_w_half = canvas_size_w * 0.5;
        let canvas_size_h = self.size.1 * ratio as f64;
        let canvas_size_h_half = canvas_size_h * 0.5;
        let canvas_shape_x_min = canvas_pos_x - canvas_size_w_half;
        let canvas_shape_x_max = canvas_pos_x + canvas_size_w_half;
        let canvas_shape_y_min = canvas_pos_y - canvas_size_h_half;
        let canvas_shape_y_max = canvas_pos_y + canvas_size_h_half;
        let canvas_shape_x_min_usize = max(canvas_shape_x_min as usize, 0);
        let canvas_shape_x_max_usize = min(canvas_shape_x_max as usize, width);
        let canvas_shape_y_min_usize = max(canvas_shape_y_min as usize, 0);
        let canvas_shape_y_max_usize = min(canvas_shape_y_max as usize, height);
        let mut xr: Vec<f64> = Vec::new();
        let mut yr: Vec<f64> = Vec::new();
        for xi in canvas_shape_x_min_usize..(canvas_shape_x_max_usize + 1) {
            let x = xi as f64;
            if x < canvas_pos_x {
                if x >= canvas_shape_x_min {
                    xr.push(1.0);
                }
                else if x+1.0 < canvas_shape_x_min {
                    xr.push(0.0);
                }
                else {
                    xr.push(x+1.0-canvas_shape_x_min);
                }
            }
            else {
                if x <= canvas_shape_x_max {
                    xr.push(1.0);
                }
                else if x-1.0 > canvas_shape_x_max {
                    xr.push(0.0);
                }
                else {
                    xr.push(canvas_shape_x_max-x+1.0);
                }
            }
        }
        for yi in canvas_shape_y_min_usize..(canvas_shape_y_max_usize + 1) {
            let y = yi as f64;
            if y < canvas_pos_y {
                if y >= canvas_shape_y_min {
                    yr.push(1.0);
                }
                else if y+1.0 < canvas_shape_y_min {
                    yr.push(0.0);
                }
                else {
                    yr.push(y+1.0-canvas_shape_y_min);
                }
            }
            else {
                if y <= canvas_shape_y_max {
                    yr.push(1.0);
                }
                else if y-1.0 > canvas_shape_y_max {
                    yr.push(0.0);
                }
                else {
                    yr.push(canvas_shape_y_max-y+1.0);
                }
            }
        }
        let mut entries: Vec<Entry<C>> = Vec::new();
        for xi in canvas_shape_x_min_usize..canvas_shape_x_max_usize {
            for yi in canvas_shape_y_min_usize..canvas_shape_y_max_usize {
                let mut cc = self.color.clone();
                cc.set_opacity(xr[xi-canvas_shape_x_min_usize].min(xr[xi+1-canvas_shape_x_min_usize]) * yr[yi-canvas_shape_y_min_usize].min(yr[yi+1-canvas_shape_y_min_usize]));
                entries.push(Entry::<C>{
                    color: cc,
                    pos: UsizeVec2::new(xi, yi),
                })
            }
        }
        entries
    }
}
impl<C: canvas::Color> RectShape<C> {
    pub fn new(c: C, x: f64, y: f64, w: f64, h: f64) -> Self {
        Self {
            color: c,
            pos: F64Vec2::new(x, y),
            size: F64Vec2::new(w, h),
        }
    }
}

pub struct CircleShape<C: canvas::Color> {
    color: C,
    pos: F64Vec2, // pos in scene coordinate
    radius: f64, // radius in scene coordinate
}
impl<C: canvas::Color> Shape<C> for CircleShape<C> {
    fn rasterize(&self, width: usize, height: usize, ratio: usize) -> Vec<Entry<C>> {
        let canvas_pos_x = self.pos.0 * ratio as f64;
        let canvas_pos_y = self.pos.1 * ratio as f64;
        let canvas_radius = self.radius * ratio as f64;
        Vec::new() // TODO
    }
}

pub struct Renderer<C: canvas::Color> {
    width: usize, // canvas_width
    height: usize, // canvas height
    ratio: usize, // canvas ratio
    bg_color: C,
}
impl<C: canvas::Color> Renderer<C> {
    pub fn new(width: usize, height: usize, ratio: usize) -> Self {
        Self {
            width,
            height,
            ratio,
            bg_color: C::default(),
        }
    }

    pub fn draw(&self, shapes: &Vec<Box<dyn Shape<C>>>) -> Vec<C> {
        let mut color_buff = nd::Array2::<Vec<Box<C>>>::from_elem((self.width, self.height), Vec::new());
        for shape in shapes {
            let entries = shape.rasterize(self.width, self.height, self.ratio);
            for entry in entries {
                let i = entry.pos.0;
                let j = entry.pos.1;
                color_buff[[i, j]].push(Box::new(entry.color));
            }
        }
        let mut res: Vec<C> = Vec::new();
        for c in color_buff.map(|colors|{
            *(C::stack_colors_with_default(colors.to_vec(), Box::new(self.bg_color.clone())))
        }).t().to_owned().iter() {
            res.push(c.clone());
        }
        res
    }
}