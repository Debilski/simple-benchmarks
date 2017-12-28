#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {}
}

extern crate time;

pub struct AutoTimer {
    start: time::SteadyTime,
}
impl AutoTimer {
    pub fn new() -> AutoTimer {
        AutoTimer {
            start: time::SteadyTime::now(),
        }
    }
}
impl Drop for AutoTimer {
    fn drop(&mut self) {
        let end = time::SteadyTime::now();
        println!("{} seconds.", end - self.start);
    }
}

use std::ops::Index;
use std::ops::IndexMut;

pub struct Array2D {
    size_x: usize,
    size_y: usize,
    data: Vec<f64>,
}

impl Array2D {
    pub fn new(size_x: usize, size_y: usize) -> Array2D {
        Array2D {
            size_x: size_x,
            size_y: size_y,
            data: vec![0.0; size_x * size_y],
        }
    }
    #[inline]
    pub fn size_x(&self) -> usize {
        self.size_x
    }
    #[inline]
    pub fn size_y(&self) -> usize {
        self.size_y
    }
}

type Tuple2d = (usize, usize);

impl Index<Tuple2d> for Array2D {
    type Output = f64;
    #[inline]
    fn index(&self, idx: Tuple2d) -> &f64 {
        let (i, j) = idx;
        &self.data[i + self.size_x * j]
    }
}

impl Index<usize> for Array2D {
    type Output = f64;
    #[inline]
    fn index(&self, idx: usize) -> &f64 {
        &self.data[idx]
    }
}

impl IndexMut<Tuple2d> for Array2D {
    #[inline]
    fn index_mut(&mut self, idx: Tuple2d) -> &mut f64 {
        let (i, j) = idx;
        &mut self.data[i + self.size_x * j]
    }
}

impl IndexMut<usize> for Array2D {
    #[inline]
    fn index_mut(&mut self, idx: usize) -> &mut f64 {
        &mut self.data[idx]
    }
}
