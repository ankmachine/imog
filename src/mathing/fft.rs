use ndarray::{Array1, ArrayView1, s};
use std::f64::consts::PI;
use num_complex::*;
use num_traits::{Zero, One};


pub fn fft_mul(first_list:Vec<f64>, second_list:Vec<f64>, decimal_places:usize) -> Vec<f64> {
        let first_list_array = Array1::from_vec(first_list);
        let second_list_array = Array1::from_vec(second_list);
        let first_list_deg = first_list_array.len();
        let second_list_deg = second_list_array.len();
        /*
        need fix this
        if q_deg == 0 || p_deg == 0 {
            return self.multiply(q)
        }
        */
        // q.deg + p.deg + 1 = length of the output
        let target_len = smallest_pow_2(first_list_deg + second_list_deg + 1);
        //
        let new_first_list_array = with_leading_zeros(first_list_array.view(), target_len);
        let new_second_list_array = with_leading_zeros(second_list_array.view(), target_len);
        // pointwise multiplication, then apply inverse
        let product_ptwise = fft(new_first_list_array.view()) * fft(new_second_list_array.view());
        let new_coeff_array = inverse_fft(product_ptwise.view())/(target_len as f64);
        // extract real parts and return, round real parts to 5 decimals.
        let rounding_factor = (10.).powi(decimal_places as i32);
        no_leading_zeros(new_coeff_array.map(|z| ((z.re()*rounding_factor).trunc())/rounding_factor).to_vec())
    }


fn fft(input_list:ArrayView1<f64>) -> Array1<Complex64> {
    // returns the value representation of p,
    let list_len = input_list.len();
    if list_len == 1 {
        return Array1::from_elem(1, Complex::new(input_list[0], 0.))
    }
    let deg = 2.0 * PI / (list_len as f64);
    let w_n = Complex64::new(deg.cos(), deg.sin());
    let even = input_list.slice(s![..;2]);
    let odd = input_list.slice(s![1..input_list.len();2]);
    let y_e = fft(even);
    let y_o = fft(odd);
    let mut y = Array1::from_elem(list_len, Complex64::zero());
    let mut w = Complex64::one();
    let half = list_len/2;
    for j in 0..half {
        let odd_term = w * y_o[j];
        y[j] = y_e[j] + odd_term;
        y[j + half] = y_e[j] - odd_term;
        w *= w_n;
    }
    y
}

fn inverse_fft(input_list:ArrayView1<Complex64>) -> Array1<Complex64> {

        let list_len = input_list.len();
        if list_len == 1 {
            return Array1::from_elem(1, input_list[0])
        }

        let deg = -2.0 * PI / (list_len as f64);
        let w_n = Complex64::new(deg.cos(), deg.sin());
        let even = input_list.slice(s![..;2]);
        let odd = input_list.slice(s![1..input_list.len();2]);
        let y_e = inverse_fft(even);
        let y_o = inverse_fft(odd);
        let mut y = Array1::from_elem(list_len, Complex64::zero());
        let mut w = Complex64::one();
        let half = list_len/2;
        for j in 0..half {
            let odd_term = w * y_o[j];
            y[j] = y_e[j] + odd_term;
            y[j + half] = y_e[j] - odd_term;
            w *= w_n;
        }
        y
    }

fn no_leading_zeros(mut input_list:Vec<f64>) -> Vec<f64> {
    if input_list.len() == 0 {
        panic!("Cannot generate polynomial from empty vec.")
    }
    while input_list.len() > 1 {
        let last = input_list.last().unwrap();
        if *last == 0.0 {
            input_list.pop();
        } else {
            break
        }
    }
    input_list
}

fn smallest_pow_2(n:usize) -> usize {
    // smallest power of 2 that is >= n
    let mut new_val:usize = 1;
    while new_val < n {
        new_val <<= 1;
    }
    new_val
}

fn with_leading_zeros(input_list:ArrayView1<f64>, target_len:usize) -> Array1<f64>{
    let mut new_array = Array1::from_elem(target_len, 0.);
    for i in 0..input_list.len() {
        new_array[i] = input_list[i];
    }
    new_array
}
