use rug::Integer;

use MPQS::{check_is_divisor, time};
use MPQS::message_MPQS::mpqs;

fn main() {
    let n = "676292275716558246502605230897191366469551764092181362779759"
        .parse::<Integer>()
        .unwrap();
    let r = time(|| mpqs(&n));
    check_is_divisor(n, r);
}
