extern crate rd2d;
use rd2d::Array2D;

#[derive(Debug, Copy, Clone)]
struct Config {
    lambda: f64,
    kappa: f64,
    tau: f64,
    sigma: f64,
}

fn r_step_all(dt: f64, uu: &mut Array2D, vv: &mut Array2D, c: &Config) {
    for idx in 0..(uu.size_x() * uu.size_y()) {
        let (u_up, v_up) = r_step(dt, uu[idx], vv[idx], c);
        uu[idx] = u_up;
        vv[idx] = v_up;
    }
}

#[inline(always)]
fn r_step(dt: f64, u: f64, v: f64, c: &Config) -> (f64, f64) {
    let lambda = c.lambda;
    let kappa = c.kappa;
    let tau = c.tau;
    let sigma = c.sigma;

    let du = (lambda * u - u * u * u + kappa - sigma * v) * dt;
    let dv = (1. / tau * (u - v)) * dt;

    (u + du, v + dv)
}

#[inline]
fn positive_modulo(i: usize, diff: i8, n: usize) -> usize {
    if i == n - 1 && diff == 1 {
        0
    } else if i == 0 && diff == -1 {
        n - 1
    } else {
        (i as isize + diff as isize) as usize
    }
}

fn diff_step(
    dt: f64,
    uu: &mut Array2D,
    vv: &mut Array2D,
    c: &Config,
    tmp_uu: &mut Array2D,
    tmp_vv: &mut Array2D,
) {
    let size_x = uu.size_x();
    let size_y = uu.size_y();

    for i in 0..size_x {
        for j in 0..size_y {
            tmp_uu[(i, j)] = (-4. * uu[(i, j)] + uu[(positive_modulo(i, -1, size_x), j)]
                + uu[(positive_modulo(i, 1, size_x), j)]
                + uu[(i, positive_modulo(j, -1, size_y))]
                + uu[(i, positive_modulo(j, 1, size_y))]);
            tmp_vv[(i, j)] = (-4. * vv[(i, j)] + vv[(positive_modulo(i, -1, size_x), j)]
                + vv[(positive_modulo(i, 1, size_x), j)]
                + vv[(i, positive_modulo(j, -1, size_y))]
                + vv[(i, positive_modulo(j, 1, size_y))]);
        }
    }

    for idx in 0..(uu.size_x() * uu.size_y()) {
        uu[idx] = uu[idx] + 0.0028 * tmp_uu[idx] * dt;
        vv[idx] = vv[idx] + 0.5 * tmp_vv[idx] * dt;
    }
}

use std::error::Error;

fn write_to_file(ar: Array2D, filename: String) -> Result<(), Box<Error>> {
    use std::fs::File;
    use std::io::prelude::*;
    use std::path::Path;
    use std::error::Error;
    let path = Path::new(&filename);

    let display = path.display();

    // Open a file in write-only mode, returns `io::Result<File>`
    let mut file = File::create(&path)?;

    for i in 0..ar.size_x() {
        for j in 0..ar.size_y() {
            file.write(format!("{} ", ar[(i, j)]).as_bytes());
        }
        file.write(b"\n");
    }

    Ok(())
}

#[macro_use]
extern crate clap;
use clap::{App, Arg};

extern crate indicatif;
use indicatif::{ProgressBar, ProgressStyle};

fn main() {
    let config = Config {
        lambda: 1.,
        kappa: -0.05,
        tau: 10.,
        sigma: 1.,
    };

    let matches = App::new("rd2d")
        .version("1.0")
        .arg(
            Arg::with_name("size_x")
                .long("size_x")
                .value_name("NUM")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("size_y")
                .long("size_y")
                .value_name("NUM")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("dt")
                .long("dt")
                .value_name("NUM")
                .takes_value(true)
                .required(true),
        )
        .arg(
            Arg::with_name("steps")
                .long("steps")
                .value_name("NUM")
                .takes_value(true)
                .required(true),
        )
        .get_matches();

    let dt = value_t!(matches, "dt", f64).unwrap_or_else(|e| e.exit());
    let size_x = value_t!(matches, "size_x", usize).unwrap_or_else(|e| e.exit());
    let size_y = value_t!(matches, "size_y", usize).unwrap_or_else(|e| e.exit());
    let steps = value_t!(matches, "steps", usize).unwrap_or_else(|e| e.exit());

    let mut uu = Array2D::new(size_x, size_y);
    let mut vv = Array2D::new(size_x, size_y);

    let mut tmp_uu = Array2D::new(uu.size_x(), uu.size_y());
    let mut tmp_vv = Array2D::new(vv.size_x(), vv.size_y());
    uu[(1, 1)] = 1.0;

    {
        let _timeit = rd2d::AutoTimer::new();

        let bar = ProgressBar::new(steps as u64);
        bar.set_style(
            ProgressStyle::default_bar()
                .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos:>7}/{len:7} {msg}")
                .progress_chars("##-"),
        );
        for _i in 0..steps {
            r_step_all(dt, &mut uu, &mut vv, &config);
            diff_step(dt, &mut uu, &mut vv, &config, &mut tmp_uu, &mut tmp_vv);
            bar.inc(1);
        }
        bar.finish();
    }

    write_to_file(uu, "uu.data".to_string());
    write_to_file(vv, "vv.data".to_string());
}
