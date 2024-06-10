/// Module for rendering utilization

use std::cmp::{min, max};
use ndarray as nd;
use crate::canvas;

/// Simple struction holding two data with the same type
struct Vec2<T>(T, T);
impl<T> Vec2<T> {
    fn new(x: T, y: T) -> Self {
        Self(x, y)
    }
}
/// Simple struction holding two number as usize
type UsizeVec2 = Vec2<usize>;
/// Simple struction holding two number as f64
type F64Vec2 = Vec2<f64>;

/// Data structure to locate color with position in grid coordinate
pub struct Entry<C: canvas::Color> {
    /// color
    color: C,
    /// pos in canvas coordinate
    pos: UsizeVec2,
}

/// Any customized shape type should implement this trait so as to be rendered on the canvas
///
/// Used to bridge scene visualization with canvas rendering
pub trait Shape<C: canvas::Color> {
    /// Rasterize shape to list of entries
    ///
    /// Parameter
    /// - width: canvas width
    /// - height: canvas height
    /// - ratio: ratio between canvas and scene
    fn rasterize(&self, width: usize, height: usize, ratio: usize) -> Vec<Entry<C>>;
}

/// Rectangle shape to be rendered given with
/// - color: color for rendering
/// - pos: center of the rectangle in scene coordinate
/// - size: width and height of the rectangle in scene coordinate
pub struct RectShape<C: canvas::Color> {
    color: C,
    pos: F64Vec2,
    size: F64Vec2,
}
impl<C: canvas::Color> Shape<C> for RectShape<C> {
    /// Rasterize rectangle to list of entries
    ///
    /// Parameter
    /// - width: canvas width
    /// - height: canvas height
    /// - ratio: ratio between canvas and scene
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
    /// Construct rectangle shape
    pub fn new(c: C, x: f64, y: f64, w: f64, h: f64) -> Self {
        Self {
            color: c,
            pos: F64Vec2::new(x, y),
            size: F64Vec2::new(w, h),
        }
    }
}

/// Circle shape to be rendered given with
/// - color: color for rendering
/// - pos: center of the circle in scene coordinate
/// - radius: radius of the circle in scene coordinate
pub struct CircleShape<C: canvas::Color> {
    color: C,
    pos: F64Vec2, // pos in scene coordinate
    radius: f64, // radius in scene coordinate
}
impl<C: canvas::Color> Shape<C> for CircleShape<C> {
    /// Rasterize circle to list of entries
    ///
    /// Parameter
    /// - width: canvas width
    /// - height: canvas height
    /// - ratio: ratio between canvas and scene
    fn rasterize(&self, width: usize, height: usize, ratio: usize) -> Vec<Entry<C>> {
        let canvas_pos_x = self.pos.0 * ratio as f64;
        let canvas_pos_y = self.pos.1 * ratio as f64;
        let canvas_radius = self.radius * ratio as f64;
        Vec::new() // TODO
    }
}
impl<C: canvas::Color> CircleShape<C> {
    /// Construct circle shape
    pub fn new(c: C, x: f64, y: f64, r: f64) -> Self {
        Self {
            color: c,
            pos: F64Vec2::new(x, y),
            radius: r,
        }
    }
}

/// Simple renderer to convert shapes to canvas coordinate
///
/// Included by canvas::Canvas 
pub struct Renderer<C: canvas::Color> {
    width: usize, // canvas_width
    height: usize, // canvas height
    ratio: usize, // canvas ratio
    bg_color: C,
}
impl<C: canvas::Color> Renderer<C> {
    /// Construct renderer with
    /// - width: canvas width
    /// - height: canvas height
    /// - ratio: ratio between canvas and scene
    pub fn new(width: usize, height: usize, ratio: usize) -> Self {
        Self {
            width,
            height,
            ratio,
            bg_color: C::default(),
        }
    }

    /// Convert list of shapes to canvas buffer
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