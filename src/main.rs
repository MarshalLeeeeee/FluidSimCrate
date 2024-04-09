use fluid_sim::parser;
use fluid_sim::canvas;

fn main() {
    let parser = parser::parse();

    let mut canvas = canvas::Canvas::new_by_parser(&parser);
    let width = canvas.get_width();
    let height = canvas.get_height();

    let t0 = std::time::SystemTime::now();
    let mut buffer: Vec<canvas::RGBAColor> = Vec::with_capacity(width*height);
    for _ in 0..width {
        for _ in 0..height {
            buffer.push(canvas::RGBAColor::new(0, 0, 0, 0));
        }
    }

    while canvas.is_valid() {
        let dt = std::time::SystemTime::now().duration_since(t0).expect("Expect delta time").subsec_nanos() as f64 / 1_000_000_000_f64;
        for i in 0..width {
            for j in 0..height {
                buffer[j * width + i] = canvas::RGBAColor::new(0, 0, (255.0 * dt) as u8, 0);
            }
        }
        canvas.refresh(&buffer);
    }
}
