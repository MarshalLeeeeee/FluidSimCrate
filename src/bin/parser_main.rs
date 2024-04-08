use serde::{Deserialize, Serialize};
use serde_json::{self, Value};
use fluid_sim::parser;
use fluid_sim::utils::type_of;

fn main() {
    let parser = parser::parse();
    println!("Parsed data: {}", type_of(&parser));
    if let Some(w) = parser.get("width") {
        let mut ww: usize = 0;
        println!("width: {} {}", w, type_of(w));
        println!(
            "Conf: {} {}",
            type_of(&w.as_number().unwrap()),
            type_of(&w.as_u64().unwrap()),
        );
        ww = w.as_u64().unwrap() as usize;
        println!("w: {w}");
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
