use fluid_sim::parser;
use fluid_sim::canvas;
use fluid_sim::scene;
use fluid_sim::scene::Scene;

fn main() {
    let raw_parser_maps = vec!(
        canvas::Canvas::get_config_map(),
        scene::SingleSmokeGridScene::get_config_map(),
    );
    let parser_map = parser::register(raw_parser_maps);
    let parser = parser::parse(parser_map);

    let mut canvas = canvas::Canvas::new_by_parser(&parser);
    let mut scene = scene::SingleSmokeGridScene::new_by_parser(&parser);

    while canvas.is_valid() {
        scene.sim();
        canvas.refresh(&scene.visualize_density());
    }
}
