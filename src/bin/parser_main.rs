use fluid_sim::parser;
use fluid_sim::utils::type_of;

// run with 
// cargo run --bin parser_main -- configs/debug.json --width 350 --height 300 --tick_dt 30
fn main() {
    let parser = parser::parse();
    println!("Parsed data: {}", type_of(&parser));
    println!("---------------");
    if let Some(w) = parser.get("width") {
        println!("width: {} {}", w, type_of(w));
    }
    println!("---------------");
    if let Some(h) = parser.get("height") {
        println!("height: {} {}", h, type_of(h));
    }
    println!("---------------");
    if let Some(dt) = parser.get("tick_dt") {
        println!("tick_dt: {} {}", dt, type_of(dt));
    }
}
