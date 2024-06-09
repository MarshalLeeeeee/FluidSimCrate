use fluid_sim::render;
use fluid_sim::canvas;

fn main() {
    let height = 480;
    let width = 640;
    let tick_dt = 10;
    let mut canvas = canvas::Canvas::new(width, height, tick_dt, 1);
    // let mut renderer = render::Renderer::new(width, height, 1);
    let mut shapes : Vec<Box<dyn render::Shape<canvas::RGBAColor>>> = Vec::new();
    for x in 0..width/5 {
        for y in 0..height/5 {
            let xx : f64 = 0.5_f64 + x as f64 * 5.0_f64;
            let yy : f64 = 0.5_f64 + y as f64 * 5.0_f64;
            shapes.push(Box::new(render::RectShape::new(
                canvas::RGBAColor::new(0, 255, 0),
                xx, yy,
                1_f64, 1_f64,
            )));
        }
    }
    while canvas.is_valid() {
        // canvas.refresh(&renderer.draw(&shapes));
        canvas.draw(&shapes);
    }
}