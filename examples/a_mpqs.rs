use rug::Integer;

use MPQS::message_MPQS::mpqs;
use MPQS::time;

fn main() {
    time(|| {
        mpqs(
            &"676292275716558246502605230897191366469551764092181362779759"
                .parse::<Integer>()
                .unwrap(),
        )
    });
}
