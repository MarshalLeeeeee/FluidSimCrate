use std::{thread, time};
use minifb::{Window, WindowOptions, Scale, Key};

/// Any customized color type should implement this trait
pub trait Color {
    /// Convert the customized color struct representation to u32
    fn to_u32(&self) -> u32;
}

/// RBGA color
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

/// The window for simple render of color
/// 
/// The following params are configurable
/// - width: the width of the window
/// - height: the height of the window
/// - tick_dt: the minimal update interval in millis
pub struct Canvas {
    width: usize,
    height: usize,
    buffer: Vec<u32>,
    window: Window,
    tick_ts: time::SystemTime,
    tick_dt: u64,
}
impl Canvas {
    /// Construct a canvas with configs
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
            tick_ts: time::SystemTime::now(),
            tick_dt: 100,
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
    pub fn refresh<T>(&mut self, buffer: &Vec<T>) where T:Color {
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
}
