# Blocked Bloom Filter

## Reqs 
- Compiler flags for g++:
	-   -msse4
	-   -O3
    -   -std=c++23

All the work is in src/main BTW

## Prelim:
Probably not going to implement the SIMD version of blocked Bloom filter, why?
-   It requires a complete refactor of the BloomFilter class in KeBaB which i dont have the auth to do nor want to rn
-   IDK the size of the dataset its author tested it on
-   Different Licence which idk how to deal with when copying

## Good news:
I'm "designing" the register filter myself :P (Going great btw). 
figured out the structure I want to follow and all the intermidiary processes (aka bit reversal, indexing etc etc)

## Bad news:
Spent 7 hours figuring out intel's BS pshufb and how to get it to work then tested it.
Turns out GCC with -O3 optimized away the byte swapping using bswap and is more than twice as fast (2000ns vs 600ns).
There is a \_\_builtin\_bitreverse64 but im poor and my laptop doesnt have it ig :(

On other bad news:
I bombed an RBC Capital markets interview soooo bad i just wanna die :(


Okkk thats it, I'll hopefully finish this by christmas and will try to not bother nate with it till after. oh
also im getting cats next week (for 2 weeks YaY)
