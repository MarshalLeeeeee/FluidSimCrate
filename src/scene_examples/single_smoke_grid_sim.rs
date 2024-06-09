use fluid_sim::parser;
use fluid_sim::canvas;
use fluid_sim::scene;

/// Run with 
/// cargo run --bin single_smoke_grid_sim -- json_file=./configs/single_smoke_grid.json [...]
///
/// Different configs can be put in [...], or change the corresponding json file
fn main() {
    let raw_parser_maps = vec!(
        canvas::Canvas::<canvas::RGBAColor>::get_config_map(),
        scene::single_smoke::SingleSmokeGridScene::<canvas::RGBAColor>::get_config_map(),
    );
    let parser_map = parser::register(raw_parser_maps);
    let parser = parser::parse(parser_map);

    let mut canvas = canvas::Canvas::<canvas::RGBAColor>::new_by_parser(&parser);
    let mut scene = scene::single_smoke::SingleSmokeGridScene::<canvas::RGBAColor>::new_by_parser(&parser);

    while canvas.is_valid() {
        scene.sim();
        canvas.draw(&scene.visualize_density());
    }
}
