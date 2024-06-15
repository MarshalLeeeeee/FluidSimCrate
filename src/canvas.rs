/// Module for canvas rendering in a separate window

use std::ops::{Add, Mul, Div};
use std::time;
use serde_json::{Value, Map, Number};
use minifb::{Window, WindowOptions, Scale, ScaleMode, Key};
use crate::parser;
use crate::render::{self, Shape};

/// Any customized temporary color type should implement this trait
pub trait TempColor {
    /// The corresponding color type fed to renderer
    type TypeColor;
    /// Convert temporary color to color for utilization
    fn to_c(self) -> Self::TypeColor where Self::TypeColor:Color;
}

/// Any customized color type should implement this trait
///
/// The color type should implement trait Clone
pub trait Color : std::clone::Clone {
    /// The corresponding temporary color type for caculation purpose (to prevent data overflow)
    type TypeTempColor;
    /// Convert the customized color struct representation to u32
    fn to_u32(&self) -> u32;
    /// Convert the customized color to corresponding temporary color
    fn to_tc(&self) -> Self::TypeTempColor where Self::TypeTempColor:TempColor+Add+Mul<u32>+Div<u32>;
    /// Weighted color mixture
    fn stack_colors(colors: Vec<Box<Self>>) -> Option<Box<Self>>;
    /// Weighted color mixture with default color
    ///
    /// Equivelant to stack_colors(colors).unwrap_or(default_color)
    fn stack_colors_with_default(colors: Vec<Box<Self>>, default_color: Box<Self>) -> Box<Self>;
    /// Mix with another color
    fn mix(self, other: Self) -> Self;
    /// New color
    fn new(x: u8, y: u8, z: u8) -> Self;
    /// Default color of the type
    fn default() -> Self;
    /// Copy color
    fn copy_color(&self) -> Self;
    /// Print color
    fn print_color(&self);
    /// Remove opacity to be no transparent as all
    fn remove_opacity(self) -> Self;
    /// Adjust opacity based on the current opacity
    fn set_opacity(&mut self, r: f64);
}

/// Temporary RGBA color for calculation
pub struct RGBATempColor(u32, u32, u32, u32);
impl TempColor for RGBATempColor {
    /// The corresponding color type fed to renderer
    type TypeColor = RGBAColor;
    /// Convert temporary color to color for utilization
    fn to_c(self) -> Self::TypeColor where Self::TypeColor:Color {
        RGBAColor::new_with_a(
            self.0 as u8,
            self.1 as u8,
            self.2 as u8,
            self.3 as u8,
        )
    }
}
impl Add for RGBATempColor {
    type Output = RGBATempColor;
    fn add(self, rhs: RGBATempColor) -> RGBATempColor {
        RGBATempColor(
            self.0 + rhs.0,
            self.1 + rhs.1,
            self.2 + rhs.2,
            self.3 + rhs.3,
        )
    }
}
impl Mul<u32> for RGBATempColor {
    type Output = RGBATempColor;
    fn mul(self, rhs: u32) -> RGBATempColor {
        RGBATempColor(
            self.0 * rhs,
            self.1 * rhs,
            self.2 * rhs,
            self.3 * rhs,
        )
    }
}
impl Div<u32> for RGBATempColor {
    type Output = RGBATempColor;
    fn div(self, rhs: u32) -> RGBATempColor {
        RGBATempColor(
            self.0 / rhs,
            self.1 / rhs,
            self.2 / rhs,
            self.3 / rhs,
        )
    }
}

/// RGBA color
#[derive(std::clone::Clone)]
pub struct RGBAColor(u8, u8, u8, u8); // r, g, b, a
impl Color for RGBAColor {
    /// The corresponding temporary color type for caculation purpose (to prevent data overflow)
    type TypeTempColor = RGBATempColor;
    /// Convert the customized color struct representation to u32
    fn to_u32(&self) -> u32 {
        ((self.3 as u32) << 24) | ((self.0 as u32) << 16) | ((self.1 as u32) << 8) | (self.2 as u32)
    }
    /// Convert the customized color to corresponding temporary color
    fn to_tc(&self) -> Self::TypeTempColor where Self::TypeTempColor:TempColor+Add+Mul<u32>+Div<u32> {
        RGBATempColor(self.0 as u32, self.1 as u32, self.2 as u32, self.3 as u32)
    }
    /// Weighted color mixture
    fn stack_colors(colors: Vec<Box<Self>>) -> Option<Box<Self>> {
        if colors.is_empty() {
            None
        }
        else {
            let mut tc = Self(0, 0, 0, 0).to_tc();
            let mut a_sum: u32 = 0;
            for color in colors {
                let mut tcc = color.to_tc();
                let a = tcc.3;
                tcc = tcc * a / 255;
                tc = tc + tcc;
                a_sum += a;
            }
            if a_sum > 0 {
                tc = tc * 255 / a_sum;
                Some(Box::new(tc.to_c().remove_opacity()))
            }
            else {
                None
            }
        }
    }
    /// Weighted color mixture with default color
    ///
    /// Equivelant to stack_colors(colors).unwrap_or(default_color)
    fn stack_colors_with_default(colors: Vec<Box<Self>>, default_color: Box<Self>) -> Box<Self> {
        Self::stack_colors(colors).unwrap_or(default_color)
    }
    /// Mix with another color
    fn mix(self, other: Self) -> Self {
        let self_box = Box::new(self);
        let other_box = Box::new(other);
        *(Self::stack_colors(vec![self_box, other_box]).unwrap())
    }
    /// New color
    fn new(x: u8, y: u8, z: u8) -> Self {
        Self(x, y, z, 255)
    }
    /// Default color of the type
    fn default() -> Self {
        Self(0,0,0,255)
    }
    /// Copy color
    fn copy_color(&self) -> Self {
        Self(
            self.0,
            self.1,
            self.2,
            self.3,
        )
    }
    /// Print color
    fn print_color(&self) {
        println!("RGBAColor: {} {} {} {}", self.0, self.1, self.2, self.3);
    }
    /// Remove opacity to be no transparent as all
    fn remove_opacity(self) -> Self {
        Self(
            self.0,
            self.1,
            self.2,
            255,
        )
    }
    /// Adjust opacity based on the current opacity
    fn set_opacity(&mut self, r: f64) {
        let a = self.3 as f64;
        let a = a * r;
        self.3 = a as u8;
    }
}
impl RGBAColor {
    /// Construct RGBAColor with all dimension given
    pub fn new_with_a(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self(r, g, b, a)
    }
}

/// The window for simple render of color
/// 
/// The following params are configurable
/// - width: the width of the scene
/// - height: the height of the scene
/// - tick_dt: the minimal update interval in millis
pub struct Canvas<C: Color> {
    width: usize,
    height: usize,
    buffer: Vec<u32>,
    window: Window,
    tick_ts: time::SystemTime,
    tick_dt: u64,
    renderer: render::Renderer<C>,
}
impl<C: Color> Canvas<C> {
    /// Construct a canvas from json configs
    pub fn new_by_parser(parser: &Value) -> Self {
        let ratio = parser::get_from_parser_usize(parser, "canvas_ratio");
        let w = parser::get_from_parser_usize(parser, "scene_width");
        let h = parser::get_from_parser_usize(parser, "scene_height");
        let tick_dt = parser::get_from_parser_u64(parser, "tick_dt");
        Canvas::new(w, h, tick_dt, ratio)
    }

    /// Provide default configs for Canvas
    pub fn get_config_map() -> Map<String, Value> {
        let mut m = Map::new();
        m.insert(String::from("canvas_ratio"), Value::Number(Number::from(1_usize)));
        m.insert(String::from("tick_dt"), Value::Number(Number::from(5_u64)));
        m
    }

    /// Construct a canvas with configs
    pub fn new (width: usize, height: usize, tick_dt: u64, ratio: usize) -> Self {
        let width = width * ratio;
        let height = height * ratio;
        let buffer = vec![0; width*height];
        let window = Window::new(
            "Rust Example - Display Color Data",
            width,
            height,
            WindowOptions {
                resize: true,
                scale: Scale::X2,
                scale_mode: ScaleMode::AspectRatioStretch,
                ..WindowOptions::default()
            },
        ).unwrap_or_else(|e| {
            panic!("Create canvas error: {}", e);
        });
        Self {
            width,
            height,
            buffer,
            window,
            tick_ts: time::SystemTime::now(),
            tick_dt: tick_dt,
            renderer: render::Renderer::new(width, height, ratio),
        }
    }

    /// Check if the canvas is still available, can be used as the tick condition in the main thread
    ///
    /// ```
    /// while canvas.is_valid() {...}
    /// ```
    pub fn is_valid(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(Key::Escape)
    }

    /// Refresh the buffer and update the window
    ///
    /// Refresh rate is controlled by the tick interval (as the minimal refresh rate)
    ///
    /// ```
    /// let buffer: Vec<Color> = Vec::new();
    /// ...
    /// canvas.refresh(buffer);
    /// ```
    pub fn refresh(&mut self, buffer: &Vec<C>) {
        while (std::time::SystemTime::now().duration_since(self.tick_ts).expect("Expect delta time").as_millis() as u64) < self.tick_dt {}
        self.buffer.iter_mut().zip(buffer.iter()).for_each(|(target, &ref source)| {
            *target = source.to_u32();
        });
        self.window.update_with_buffer(&self.buffer, self.width, self.height).unwrap_or_else(|e| {
            panic!("Tick canvas error: {}", e);
        });
        self.tick_ts = std::time::SystemTime::now();
    }

    /// Refresh the buffer and update the window with list of render shapes
    pub fn draw(&mut self, shapes: &Vec<Box<dyn Shape<C>>>) {
        self.refresh(&self.renderer.draw(&shapes));
    }

    /// Get the width of the canvas
    pub fn get_width(&self) -> usize {
        self.width
    }

    /// Get the height of the canvas
    pub fn get_height(&self) -> usize {
        self.height
    }
}

