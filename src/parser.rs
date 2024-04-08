/// Module for parse configuration

use std::env;
use std::fs::File;
use std::io::BufReader;
use serde::{Deserialize, Serialize};
use serde_json::Value;

/// RawParser defines all configs
///  - width(usize): The width of the canvas.
///  - height(usize): The height of the canvas.
///  - tick_dt(u64): The minimal refresh rate of the canvas.
#[derive(Deserialize, Serialize)]
struct RawParser {
    width: usize,
    height: usize,
    tick_dt: u64,
}

// TODO: encapsulate
struct Parser {
    json: Value,

}

/// The default values are provided
pub fn parse() -> Value {
    let raw_parser = RawParser {
        width: 480,
        height: 640,
        tick_dt: 100,
    };
    let parser_str = serde_json::to_string(&raw_parser).expect("Failed to serialize JSON");
    let parser = serde_json::from_str(&parser_str).expect("Failed to parse JSON");
    let (file_name, cmd_keys, cmd_values) = parse_from_cmd();

    // read from json
    if let Ok(file) = File::open(file_name) {
        let reader = BufReader::new(file);
        if let Ok(json_data) = serde_json::from_reader::<BufReader<File>, Value>(reader) {
            for (key, value) in json_data.as_object().unwrap() {
                println!("conf {:?} {:?}", key, value); // TODO
            }
        }
    }

    // read from cmd
    let k_v: Vec<_> = cmd_keys.iter().zip(cmd_values.iter()).collect();
    for (k, v) in k_v {
        println!("conf k {} v {}", k, v);
    }

    parser
}

// enum CmdTypeEnum {
//     CmdString,
//     CmdU64,
//     CmdUsize,
//     CmdD64,
// }

// const CMD_TYPE: HashMap<&str, CmdTypeEnum> =  {
//     let mut map = HashMap::new();
//     map.insert("config_json", CmdTypeEnum::CmdString);
//     map.insert("width", CmdTypeEnum::CmdUsize);
//     map
// };

fn parse_from_cmd() -> (String, Vec<String>, Vec<String>) {
    let cmd_args: Vec<String> = env::args().collect();
    let cmd_arg_cnt = cmd_args.len();
    let mut file_name = "";
    let mut offset = 1;
    if cmd_arg_cnt % 2 == 0 {
        if let Some(arg) = cmd_args.get(1) {
            file_name = arg;
        }
        offset = 2;
    }
    let mut vec_cmd_k: Vec<String> = Vec::new();
    let mut vec_cmd_v: Vec<String> = Vec::new();
    for i in 0..((cmd_arg_cnt-1) / 2) {
        let i1 = 2 * i + offset;
        let i2 = 2 * i + offset + 1;
        if let Some(k) = cmd_args.get(i1) {
            if let Some(v) = cmd_args.get(i2) {
                vec_cmd_k.push(String::from(&k[2..]));
                vec_cmd_v.push(String::from(&v[..]));
            }
        }
    }
    (String::from(file_name), vec_cmd_k, vec_cmd_v)
}
