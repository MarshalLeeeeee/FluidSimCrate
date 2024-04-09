/// Module for parse configuration

use std::env;
use std::fs::File;
use std::io::BufReader;
use serde_json::{Map, Value, json};

pub fn register(vec_map: Vec<Map<String, Value>>) -> Map<String, Value> {
    let mut json_map = Map::new();
    for m in vec_map {
        merge_map(&mut json_map, m)
    }
    json_map
}

fn merge_map(map_target: &mut Map<String, Value>, map_source: Map<String, Value>) {
    for (k, v) in map_source.iter() {
        map_target.insert(k.to_string(), v.clone());
    }
}

/// Return the configuration
///
/// The default configs are first provided
///
/// User can change the configs by command line in the following format
///
/// cargo run ... -- [json_file_path] --k1 v1 --k2 2 ...
///  - ```[json_file_path]``` is an optional cmd. When the length of user command is odd, the first user command will be considered as json_file_path
///  - ```--k1 v1``` should be paired, so that certain config with name k1 will be updated to v1
///
/// The override priority is: cmd_args > json > default
///
/// Args that are not defined in default configs will be not be added.
///
/// The return value is serde_json::Value
pub fn parse(map: Map<String, Value>) -> Value {
    // The default values are provided here
    // let mut json: Value = json!({
    //     "width": 480,
    //     "height": 640,
    //     "tick_dt": 100,
    // });
    let mut json = serde_json::to_value(map).expect("Should have collected default configs");
    let (file_name, cmd_keys, cmd_values) = parse_from_cmd();

    // read from json
    if let Ok(file) = File::open(file_name) {
        let reader = BufReader::new(file);
        if let Ok(json_data) = serde_json::from_reader::<BufReader<File>, Value>(reader) {
            json_value_update(&mut json, json_data);
        }
    }

    // read from cmd
    let k_v: Vec<_> = cmd_keys.iter().zip(cmd_values.iter()).collect();
    let mut json_map = Map::new();
    for (k, v) in k_v {
        if let Ok(num) = v.parse::<u64>() {
            json_map.insert(k.clone(), Value::Number(serde_json::Number::from(num)));
        }
        else {
            json_map.insert(k.clone(), Value::String(v.clone()));
        }
    }
    if let Ok(json_cmd) = serde_json::to_value(json_map) {
        json_value_update(&mut json, json_cmd);
    }

    json
}

/// Simple parse of the user commands
///
/// Return 
/// - json_file_path, empty when no json_file_path is provided according to the rules (in method parse)
/// - vector of key names
/// - vector of value (before parsed)
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

/// Update value of v2 to v1, if the key of v2 exists in v1
fn json_value_update(v1: &mut Value, v2: Value) {
    for (key, value) in v2.as_object().unwrap() {
        if let Some(v) = v1.get_mut(key) {
            *v = value.clone();
        }
    }
}
