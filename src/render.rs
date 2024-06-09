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

pub struct Entry<C: canvas::Color + std::clone::Clone> {
    color: C,
    pos: UsizeVec2, // pos in canvas coordinate
}
impl<C: canvas::Color + std::clone::Clone> Entry<C> {
    fn print_self(&self) {
        self.color.print_color();
        println!("Canvas pos: {} {}", self.pos.0, self.pos.1);
    }
}

pub trait Shape<C: canvas::Color + std::clone::Clone> {
    fn rasterize(&self, width: usize, height: usize, ratio: usize) -> Vec<Entry<C>>;
}

pub struct RectShape<C: canvas::Color + std::clone::Clone> {
    color: C,
    pos: F64Vec2, // pos in scene coordinate
    size: F64Vec2, // size in scene coordinate
}
impl<C: canvas::Color + std::clone::Clone> Shape<C> for RectShape<C> {
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
        let mut xr: Vec<f64> = Vec::new();
        let mut yr: Vec<f64> = Vec::new();
        let mut xr_min_i: i32 = -1;
        let mut xr_max_i: i32 = -1;
        let mut yr_min_i: i32 = -1;
        let mut yr_max_i: i32 = -1;
        for xi in 0..(width+1) {
            let x = xi as f64;
            if x < canvas_pos_x {
                if x >= canvas_shape_x_min {
                    xr.push(1.0);
                    if xr_min_i == -1 {
                        xr_min_i = xi as i32;
                    }
                    if xr_min_i != -1 {
                        xr_max_i = xi as i32;
                    }
                }
                else if x+1.0 < canvas_shape_x_min {
                    xr.push(0.0);
                }
                else {
                    xr.push(x+1.0-canvas_shape_x_min);
                    if xr_min_i == -1 {
                        xr_min_i = xi as i32;
                    }
                    if xr_min_i != -1 {
                        xr_max_i = xi as i32;
                    }
                }
            }
            else {
                if x <= canvas_shape_x_max {
                    xr.push(1.0);
                    if xr_min_i == -1 {
                        xr_min_i = xi as i32;
                    }
                    if xr_min_i != -1 {
                        xr_max_i = xi as i32;
                    }
                }
                else if x-1.0 > canvas_shape_x_max {
                    xr.push(0.0);
                }
                else {
                    xr.push(canvas_shape_x_max-x+1.0);
                    if xr_min_i == -1 {
                        xr_min_i = xi as i32;
                    }
                    if xr_min_i != -1 {
                        xr_max_i = xi as i32;
                    }
                }
            }
        }
        for yi in 0..(height+1) {
            let y = yi as f64;
            if y < canvas_pos_y {
                if y >= canvas_shape_y_min {
                    yr.push(1.0);
                    if yr_min_i == -1 {
                        yr_min_i = yi as i32;
                    }
                    if yr_min_i != -1 {
                        yr_max_i = yi as i32;
                    }
                }
                else if y+1.0 < canvas_shape_y_min {
                    yr.push(0.0);
                }
                else {
                    yr.push(y+1.0-canvas_shape_y_min);
                    if yr_min_i == -1 {
                        yr_min_i = yi as i32;
                    }
                    if yr_min_i != -1 {
                        yr_max_i = yi as i32;
                    }
                }
            }
            else {
                if y <= canvas_shape_y_max {
                    yr.push(1.0);
                    if yr_min_i == -1 {
                        yr_min_i = yi as i32;
                    }
                    if yr_min_i != -1 {
                        yr_max_i = yi as i32;
                    }
                }
                else if y-1.0 > canvas_shape_y_max {
                    yr.push(0.0);
                }
                else {
                    yr.push(canvas_shape_y_max-y+1.0);
                    if yr_min_i == -1 {
                        yr_min_i = yi as i32;
                    }
                    if yr_min_i != -1 {
                        yr_max_i = yi as i32;
                    }
                }
            }
        }
        let mut entries: Vec<Entry<C>> = Vec::new();
        for xi in (xr_min_i as usize)..(xr_max_i as usize - 1) {
            for yi in (yr_min_i as usize)..(yr_max_i as usize - 1) {
                let mut cc = self.color.clone();
                cc.set_opacity(xr[xi].min(xr[xi+1]) * yr[yi].min(yr[yi+1]));
                entries.push(Entry::<C>{
                    color: cc,
                    pos: UsizeVec2::new(xi, yi),
                })
            }
        }
        entries
    }
}
impl<C: canvas::Color + std::clone::Clone> RectShape<C> {
    pub fn new(c: C, x: f64, y: f64, w: f64, h: f64) -> Self {
        Self {
            color: c,
            pos: F64Vec2::new(x, y),
            size: F64Vec2::new(w, h),
        }
    }
}

pub struct CircleShape<C: canvas::Color + std::clone::Clone> {
    color: C,
    pos: F64Vec2, // pos in scene coordinate
    radius: f64, // radius in scene coordinate
}

pub struct Renderer<C: canvas::Color + std::clone::Clone> {
    width: usize, // canvas_width
    height: usize, // canvas height
    ratio: usize, // canvas ratio
    bg_color: C,
}
impl<C: canvas::Color + std::clone::Clone> Renderer<C> {
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