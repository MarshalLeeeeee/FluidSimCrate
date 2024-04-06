use std::{thread, time};
use minifb::{Window, WindowOptions, Scale, Key};


pub trait Color {
    fn to_u32(&self) -> u32;
}

pub struct RGBAColor(u8, u8, u8, u8); // r, g, b, a
impl RGBAColor {
    pub fn new(r: u8, g: u8, b: u8, a: u8) -> Self {
        Self(r, g, b, a)
    }
}
impl Color for RGBAColor {
    fn to_u32(&self) -> u32 {
        ((self.3 as u32) << 24) | ((self.0 as u32) << 16) | ((self.1 as u32) << 8) | (self.2 as u32)
    }
}

pub struct Canvas {
    width: usize,
    height: usize,
    buffer: Vec<u32>,
    window: Window,
    tick_ts: f64,
    tick_dt: u32,
}

impl Canvas {
    pub fn new (width: usize, height: usize) -> Self {
        let buffer = vec![0; width*height];
        let window = Window::new(
            "Rust Example - Display Color Data",
            width,
            height,
            WindowOptions {
                resize: false,
                scale: Scale::X1,
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
            tick_ts: 0.0,
            tick_dt: 100,
        }
    }
    pub fn is_valid(&self) -> bool {
        self.window.is_open() && !self.window.is_key_down(Key::Escape)
    }
    pub fn refresh<T>(&mut self, buffer: &Vec<T>) where T:Color {
        thread::sleep(time::Duration::from_millis(100)); // TODO: const
        // self.buffer.clone_from_slice(buffer);
        self.buffer.iter_mut().zip(buffer.iter()).for_each(|(target, &ref source)| {
            *target = source.to_u32();
        });
        self.window.update_with_buffer(&self.buffer, self.width, self.height).unwrap_or_else(|e| {
            panic!("Tick canvas error: {}", e);
        });
    }
}
