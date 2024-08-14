use ndarray::{Array1, ArrayView1, s};
use std::f64::consts::PI;
use num_complex::*;
use num_traits::{Zero, One};

fn main() {
    println!("Hello, world!");
    let x = Array1::from_vec(vec![1f64, 2.0, 3.0, 4.0]);
    println!("{:?}",fft(x.view()));
    let p1 = Array1::from_vec(vec![-1.,1.]).pow(4);
    let p2 = Array1::from_vec(vec![-1.,1.]);
    let fft = p1.fft_mul(&p2, 10);
    // assert_eq!(p1.multiply(&p2), fft);
    println!("{:?}", fft)
}

fn fft(p:ArrayView1<f64>) -> Array1<Complex64> {
    // returns the value representation of p,
    let n = p.len();
    if n == 1 {
        return Array1::from_elem(1, Complex::new(p[0], 0.))
    }
    let deg = 2.0 * PI / (n as f64);
    let w_n = Complex64::new(deg.cos(), deg.sin());
    let even = p.slice(s![..;2]);
    let odd = p.slice(s![1..p.len();2]);
    let y_e = fft(even);
    let y_o = fft(odd);
    let mut y = Array1::from_elem(n, Complex64::zero());
    let mut w = Complex64::one();
    let half = n/2;
    for j in 0..half {
        let odd_term = w * y_o[j];
        y[j] = y_e[j] + odd_term;
        y[j + half] = y_e[j] - odd_term;
        w *= w_n;
    }
    y
}

fn inverse_fft(p:ArrayView1<Complex64>) -> Array1<Complex64> {

        let n = p.len();
        if n == 1 {
            return Array1::from_elem(1, p[0])
        }

        let deg = -2.0 * PI / (n as f64);
        let w_n = Complex64::new(deg.cos(), deg.sin());
        let even = p.slice(s![..;2]);
        let odd = p.slice(s![1..p.len();2]);
        let y_e = inverse_fft(even);
        let y_o = inverse_fft(odd);
        let mut y = Array1::from_elem(n, Complex64::zero());
        let mut w = Complex64::one();
        let half = n/2;
        for j in 0..half {
            let odd_term = w * y_o[j];
            y[j] = y_e[j] + odd_term;
            y[j + half] = y_e[j] - odd_term;
            w *= w_n;
        }
        y
    }

fn no_leading_zeros(mut c:Vec<f64>) -> Array1<f64> {
    if c.len() == 0 {
        panic!("Cannot generate polynomial from empty vec.")
    }
    while c.len() > 1 {
        let last = c.last().unwrap();
        if *last == 0.0 {
            c.pop();
        } else {
            break
        }
    }
    Array1::from_vec(c)
}

fn smallest_pow_2(n:usize) -> usize {
    // smallest power of 2 that is >= n
    let mut new_val:usize = 1;
    while new_val < n {
        new_val <<= 1;
    }
    new_val
}

fn with_leading_zeros(q:ArrayView1<f64>, target_len:usize) -> Array1<f64>{
    let mut new_array = Array1::from_elem(target_len, 0.);
    for i in 0..q.len() {
        new_array[i] = q[i];
    }
    new_array
}

fn fft_mul(p:ArrayView1<f64>, q:ArrayView1<f64>, decimal_places:usize) -> Array1<f64> {
        let q_deg = q.len();
        let p_deg = p.len();
        /*
        need fix this
        if q_deg == 0 || p_deg == 0 {
            return self.multiply(q)
        }
        */
        // q.deg + p.deg + 1 = length of the output
        let target_len = smallest_pow_2(p_deg + q_deg + 1);
        //
        let new_p = with_leading_zeros(p, target_len);
        let new_q = with_leading_zeros(q, target_len);
        // pointwise multiplication, then apply inverse
        let product_ptwise = fft(new_p.view()) * fft(new_q.view());
        let new_coeff_array = inverse_fft(product_ptwise.view())/(target_len as f64);
        // extract real parts and return, round real parts to 5 decimals.
        let rounding_factor = (10.).powi(decimal_places as i32);
        no_leading_zeros(new_coeff_array.map(|z| ((z.re()*rounding_factor).trunc())/rounding_factor).to_vec())
    }
