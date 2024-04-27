use fluid_sim::parser;
use fluid_sim::utils::type_of;
use serde_json::{Value, Map, Number};

// run with 
// cargo run --bin parser_main -- width=350 height=300 tick_dt=30
fn main() {
    let mut m = Map::new();
    m.insert(String::from("width"), Value::Number(Number::from(480_usize)));
    m.insert(String::from("height"), Value::Number(Number::from(640_usize)));
    m.insert(String::from("tick_dt"), Value::Number(Number::from(100_u64)));

    let parser = parser::parse(m);
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
