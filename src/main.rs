mod mathing;
use crate::mathing::fft;

fn main() {
    println!("Hello, world!");
    let p1 = vec![-1.,1.];
    let p2 = vec![-1.,1.];
    let fft = fft::fft_mul(p1, p2, 10);
    assert_eq!(vec![1.0, -2.0, 1.0], fft);
    println!("{:?}", fft)
}
