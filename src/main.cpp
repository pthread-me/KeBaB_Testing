#include <climits>
#include <print>
#include <cstdint>
#include <vector>
#include <tmmintrin.h>
#include <bloom_filter.hpp>



using namespace std;
using word_t = uint64_t;

const word_t MAX_WORD = ULLONG_MAX;
const word_t WORD_LEN = ULLONG_WIDTH;
const __m64 INTEL_PSHUFB_FLAG = __m64(0x1020304050607); 

template<typename Hash = kebab::MultiplyShift>
void init(size_t elements, double error_rate, size_t num_hashes, kebab::MultiplyShift& hash, size_t& bits, size_t& set_bits, vector<word_t>& filter) {
    bits = 256;
    set_bits = 0;

    // The filter should have a protected boundary word which can not be indexed to
    // but can be part of a block
    filter = std::vector<word_t>(bits/WORD_LEN + WORD_LEN, 0ULL);
  
    hash = Hash(bits); 
}



// TESTED :P
word_t get_word(vector<word_t>& filter, size_t hash){
  if (hash % WORD_LEN == 0){
    return filter[hash/WORD_LEN];
  } 

  size_t b1 = filter[hash/WORD_LEN];
  size_t b2 = filter[hash/WORD_LEN + 1];
  size_t i = hash%WORD_LEN;
  size_t j = WORD_LEN - i;

  size_t mask1 = MAX_WORD << i; 
  size_t mask2 = MAX_WORD >> j; 

  return ((b2 & mask2)<<j) | ((b1 & mask1) >> i);
  
}




word_t reverse_bits(uint64_t word){

// Unused but keep for ref
# if defined __builin_bitreverse64
  return __builtin_bitreverse64(word);
#else
  word = ((word & 0xaaaaaaaaaaaaaaaa)>>1) | ((word & 0x5555555555555555) << 1); 
  word = ((word & 0xcccccccccccccccc)>>2) | ((word & 0x3333333333333333) << 2); 
  word = ((word & 0xf0f0f0f0f0f0f0f0)>>4) | ((word & 0x0f0f0f0f0f0f0f0f) << 4); 
  
  /*
   All instructions bellow SHOULD be optimized out with -O3, and are left
   to illustrate intent to the compiler. If not optimized add the -msse4 flag
   and replace with:

    const __m64 INTEL_PSHUFB_FLAG = __m64(0x1020304050607); 
    word = (uint64_t)( _mm_shuffle_pi8(__m64(word), INTEL_PSHUFB_FLAG));
  */
  word = ((word & 0xff00ff00ff00ff00)>>8) | ((word & 0x00ff00ff00ff00ff) << 8);  
  word = ((word & 0xffff0000ffff0000)>>16) | ((word & 0x0000ffff0000ffff) << 16);
  return (word >> 32) | (word << 32);
#endif
}




size_t set_word(vector<word_t>& filter, size_t hash, word_t word){
  if (!word) return 0; 

  if (hash % WORD_LEN == 0){
    size_t prev_bitcount = std::popcount(filter[hash/WORD_LEN + 1]);
    filter[hash/WORD_LEN + 1] |= word;
    return std::popcount(filter[hash/WORD_LEN + 1]) - prev_bitcount;
  } 

  size_t i = hash%WORD_LEN;
  return 0;
  /*
  size_t mask1 = word << i;
  size_t mask2 = __builtin_bitreverse64(word << (WORD_LEN - i));
  
  
  size_t prev_setbits = std::popcount(filter[hash/WORD_LEN + 1]) + std::popcount(filter[WORD_LEN]);

  filter[hash/WORD_LEN +1] |= (word & mask1); 
  filter[hash/WORD_LEN] |= (word & mask2); 

  return std::popcount(filter[hash/WORD_LEN + 1]) + std::popcount(filter[WORD_LEN]) - prev_setbits;
  */
}


int main(){

  size_t bits;
  size_t set_bits;
  vector<word_t> filter;

  kebab::MultiplyShift hash;
  init(300, 0.1, 3, hash, bits, set_bits, filter);

  filter[0] = 0xa000000000000000;
  filter[1] = 0xE000000000000003;
  

//  set_bits += set_word(filter, 66, MAX_WORD);
  size_t w = get_word(filter, 61);
  println("{}", w);

//  uint64_t a = 4;
//  reverse_bits(a);
}
