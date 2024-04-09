/// Module for some command utils which are clearly belong to any functional module

pub fn type_of<T>(_: & T) -> &'static str {
    std::any::type_name::<T>()
}
