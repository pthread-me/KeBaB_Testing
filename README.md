# Benchmark Results of SIMD-ed AES vs MultiplyHash

## ASM 
- input the contents of comp_exp.txt to compiler explorer with the follwing flags:
    -   -march=native
	-   -maes
	-   -msse4
	-   -mvaes
	-   -O3
    -   -std=c++23

And the latest gcc compiler


## Prelim:
I think the main bottleneck for the aes hashing are the actual rounds. If true,
then this would be good, however i'm almost certain that less rounds => more concentrated
hash values => more collisions.

I've tested collisions with 2 and 10 rounds with the later being much better.
However the performance of the aes hash is worse than MulHash by ~10% when using 10 rounds
and is only better with 3 rounds. So i'll need to test the collisions with 3 aesenc rounds.

There is also the sha intrinsic, but i dont have access to hardware support for it yet.



