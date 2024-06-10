/// Module for simple render

use std::ops::{Add, Mul, Div};
use std::{thread, time};
use serde_json::{Value, Map, Number};
use minifb::{Window, WindowOptions, Scale, ScaleMode, Key};
use crate::parser;
use crate::render::{self, Shape};

pub trait TempColor {
    type TypeColor;
    fn to_c(self) -> Self::TypeColor where Self::TypeColor:Color;
}

/// Any customized color type should implement this trait
pub trait Color : std::clone::Clone {
    type TypeTempColor;
    /// Convert the customized color struct representation to u32
    fn to_u32(&self) -> u32;
    fn to_tc(&self) -> Self::TypeTempColor where Self::TypeTempColor:TempColor+Add+Mul<u32>+Div<u32>;
    fn remove_opacity(self) -> Self;
    fn set_opacity(&mut self, r: f64);
    fn _do_stack_color(colors: Vec<Box<Self>>) -> Option<Box<Self>>;
    fn stack_colors(colors: Vec<Box<Self>>) -> Option<Box<Self>>;
    fn stack_colors_with_default(colors: Vec<Box<Self>>, default_color: Box<Self>) -> Box<Self>;
    fn default() -> Self;
    fn copy_color(&self) -> Self;
    fn print_color(&self);
    fn new(x: u8, y: u8, z: u8) -> Self;
    fn mix(self, other: Self) -> Self;
}

/// RGBA color for calculation
pub struct RGBATempColor(u32, u32, u32, u32);
impl TempColor for RGBATempColor {
    type TypeColor = RGBAColor;
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
    type TypeTempColor = RGBATempColor;
    fn to_u32(&self) -> u32 {
        ((self.3 as u32) << 24) | ((self.0 as u32) << 16) | ((self.1 as u32) << 8) | (self.2 as u32)
    }
    fn remove_opacity(self) -> Self {
        Self(
            self.0,
            self.1,
            self.2,
            255,
        )
    }
    fn set_opacity(&mut self, r: f64) {
        let a = self.3 as f64;
        let a = a * r;
        self.3 = a as u8;
    }
    fn to_tc(&self) -> Self::TypeTempColor where Self::TypeTempColor:TempColor+Add+Mul<u32>+Div<u32> {
        RGBATempColor(self.0 as u32, self.1 as u32, self.2 as u32, self.3 as u32)
    }
    fn _do_stack_color(colors: Vec<Box<Self>>) -> Option<Box<Self>> {
        assert!(!colors.is_empty(), "Colors should not be empty");
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
    fn stack_colors(colors: Vec<Box<Self>>) -> Option<Box<Self>> {
        if colors.is_empty() {
            None
        }
        else {
            Self::_do_stack_color(colors)
        }
    }
    fn stack_colors_with_default(colors: Vec<Box<Self>>, default_color: Box<Self>) -> Box<Self> {
        if colors.is_empty() {
            default_color
        }
        else {
            match Self::_do_stack_color(colors) {
                Some(c) => c,
                None => default_color,
            }
        }
    }
    fn default() -> Self {
        Self(0,0,0,255)
    }
    fn copy_color(&self) -> Self {
        Self(
            self.0,
            self.1,
            self.2,
            self.3,
        )
    }
    fn print_color(&self) {
        println!("RGBAColor: {} {} {} {}", self.0, self.1, self.2, self.3);
    }
    fn new(x: u8, y: u8, z: u8) -> Self {
        Self(x, y, z, 255)
    }
    fn mix(self, other: Self) -> Self {
        let self_box = Box::new(self);
        let other_box = Box::new(other);
        *(Self::_do_stack_color(vec![self_box, other_box]).unwrap())
    }
}
impl RGBAColor {
    pub fn new_with_a(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self(r, g, b, a)
    }
}

/// The window for simple render of color
/// 
/// The following params are configurable
/// - width: the width of the window
/// - height: the height of the window
/// - tick_dt: the minimal update interval in millis
pub struct Canvas<C: Color> {
    ratio: usize,
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
            ratio,
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
        let dt: u64 = std::time::SystemTime::now().duration_since(self.tick_ts).expect("Expect delta time").as_millis() as u64;
        if dt < self.tick_dt {
            thread::sleep(time::Duration::from_millis(self.tick_dt - dt));
        }
        self.buffer.iter_mut().zip(buffer.iter()).for_each(|(target, &ref source)| {
            *target = source.to_u32();
        });
        self.window.update_with_buffer(&self.buffer, self.width, self.height).unwrap_or_else(|e| {
            panic!("Tick canvas error: {}", e);
        });
        self.tick_ts = std::time::SystemTime::now();
    }

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

