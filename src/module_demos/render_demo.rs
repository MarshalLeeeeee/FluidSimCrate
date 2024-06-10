use fluid_sim::render;
use fluid_sim::canvas;

fn main() {
    let width = 640;
    let height = 480;
    let tick_dt = 10;
    let mut canvas = canvas::Canvas::new(width, height, tick_dt, 1);
    let mut shapes : Vec<Box<dyn render::Shape<canvas::RGBAColor>>> = Vec::new();
    for x in 0..width/5 {
        for y in 0..height/5 {
            let xx : f64 = 0.5_f64 + x as f64 * 5.0_f64;
            let yy : f64 = 0.5_f64 + y as f64 * 5.0_f64;
            shapes.push(Box::new(render::RectShape::new(
                canvas::RGBAColor::new_with_a(0, 255, 0, 255),
                xx, yy,
                1_f64, 1_f64,
            )));
        }
    }
    shapes.push(Box::new(render::CircleShape::new(
        canvas::RGBAColor::new_with_a(0, 0, 255, 255),
        240.0, 320.0,
        30.0,
    )));
    shapes.push(Box::new(render::CircleShape::new(
        canvas::RGBAColor::new_with_a(255, 0, 0, 120),
        120.0, 160.0,
        50.0,
    )));
    shapes.push(Box::new(render::CircleShape::new(
        canvas::RGBAColor::new_with_a(100, 120, 200, 80),
        360.0, 480.0,
        80.0,
    )));
    while canvas.is_valid() {
        canvas.draw(&shapes);
    }
}