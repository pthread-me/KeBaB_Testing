#include <climits>
#include <print>
#include <cstdint>
#include <vector>
#include <tmmintrin.h>
#include <bloom_filter.hpp>
#include <random>
#include <fstream>

using namespace std;
using word_t = uint64_t;

const word_t MAX_WORD = ULLONG_MAX;
const word_t WORD_LEN = ULLONG_WIDTH;
const __m64 INTEL_PSHUFB_FLAG = __m64(0x1020304050607); 




namespace seeds{

// Supports at most 32 hash functions
constexpr uint64_t SEED[] = {          
    0x153C67147CEBD9C1, 0xE9E9221977E2486E,
    0xBD2A5DE364F86CEC, 0xF53E63242C7C96CA,
    0xEA71F713607B8025, 0xDA1DC2E81860AC93,
    0x700FC578B9B89EFC, 0x7ED09A9433D0F542,
    0xED43BDEDBCF69432, 0x1D322B028A861DAA,
    0x6E8CDB8F04EE5FFD, 0xEC53221EFD3A5C53,
    0x01EE14F09892D967, 0xD6382ACCCBCF0420,
    0xD448F78598D09FBE, 0x922AA2623D2BF77A,
    0x4AF98D70BD02F4D9, 0xBE9A532696D539D9,
    0x57CB1CF8FA6F105D, 0x4347990C105CF57C,
    0xD5E6B9B31C51D5D6, 0x2196C4CF3D467371,
    0x78BD99C62BA864CD, 0x0B747BD60B9F2FB4,
    0xE636A63B15DC2C60, 0xE3D4C1379D7C2FF0,
    0x2B5C7FAF45C1B370, 0xFE0247B305095328,
    0xE4F3205AADABEA31, 0xD631A450CF4BA7BA,
    0x7E0034EEC6C9E610, 0xCAF71C56BB5D4B4D
};
}




template<typename Hash = kebab::MultiplyShift>
void init(size_t elements, double error_rate, size_t num_hashes, kebab::MultiplyShift& hash, size_t& bits, size_t& set_bits, vector<word_t>& filter) {
    set_bits = 0;
    num_hashes = 3;

    // The filter should have a protected boundary word which can not be indexed to
    // but can be part of a block
    filter = std::vector<word_t>(elements/WORD_LEN + WORD_LEN, 0ULL);
  
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



// TESTED :P
size_t set_word(vector<word_t>& filter, size_t hash, word_t word){
  if (!word) return 0; 
  size_t prev_bitcount;

  if (hash % WORD_LEN == 0){
    prev_bitcount = std::popcount(filter[hash/WORD_LEN]);
    filter[hash/WORD_LEN] |= word;
    return std::popcount(filter[hash/WORD_LEN]) - prev_bitcount;
  }

  prev_bitcount = std::popcount(filter[hash/WORD_LEN]) + std::popcount(filter[hash/WORD_LEN + 1]);
  size_t i = hash%WORD_LEN;
  size_t j = WORD_LEN - i;

  size_t mask1 = word << i;
  size_t mask2 = word >> j;
  
  filter[hash/WORD_LEN] |= mask1;
  filter[hash/WORD_LEN + 1] |= mask2;

  return (std::popcount(filter[hash/WORD_LEN]) + std::popcount(filter[hash/WORD_LEN + 1])) - prev_bitcount;
}


bool contains(vector<word_t>& filter, size_t hash, size_t num_hashs){
  size_t word = get_word(filter, hash);
  for (auto i=0; i<num_hashs; ++i){
    size_t pos = (hash * seeds::SEED[i]) % WORD_LEN;
    if (!((word >> pos) & static_cast<word_t>(1))) return false; 
  }
  return true;
}

void add(vector<word_t>& filter, size_t hash, size_t num_hashs){
  size_t b1 = hash/WORD_LEN; 
  if(hash % WORD_LEN == 0){
    size_t mask = 0;
    for (size_t i=0; i<num_hashs; ++i){
      size_t pos = (hash * seeds::SEED[i]) % WORD_LEN;
      size_t cur = 1ULL<<pos;
      mask |= cur;
    }
    println("Block: {}, Block before: {}", b1, filter[b1]);
    filter[b1] |= mask;
    println("mask: {}, block after: {}", mask, filter[b1]);


  }else{
    size_t start = hash%WORD_LEN;
    size_t b2 = b1 + 1; 
    size_t mask1 = 0, mask2 = 0;

    for (size_t i=0; i<num_hashs; ++i){
      size_t pos = (hash * seeds::SEED[i]) % WORD_LEN;
      println("start: {}, pos: {}", start, pos);
      if(start+pos < WORD_LEN){
        println("Added to mask1 at {}", start+pos);
        mask1 |= (1ULL<<(start+pos));
        println("Mask1: {}", mask1);
      }else{
        println("Added to mask2 at {}", pos-(WORD_LEN-start));
        mask2 |= (1ULL<<(pos-(WORD_LEN-start)));
        println("Mask2: {}", mask2);
      }
    }
    println("start: {}", start);
    println("Block1: {}, Block1 before: {}", b1, filter[b1]);
    println("Block2: {}, Block2 before: {}", b2, filter[b2]);
    filter[b1] |= mask1;
    filter[b2] |= mask2;

    println("mask1: {}, block1 after: {}", mask1, filter[b1]);
    println("mask2: {}, block2 after: {}", mask2, filter[b2]);
    println("Final word: {}", get_word(filter, hash));
  }
}


void check(){
  std::random_device rd;
  std::default_random_engine gen(rd());
  std::uniform_int_distribution<uint64_t> dist(0, 99999);
  std::uniform_int_distribution<uint64_t> dist2(0, 1);
  
  size_t bits = 0;
  size_t set_bits = 0;
  size_t num_hashs = 0;
  vector<word_t> filter;

  kebab::MultiplyShift hash;
  init(99999, 0.1, num_hashs, hash, bits, set_bits, filter);
  
  vector<size_t> hashes;
  vector<size_t> checks;
  


  for(size_t i=0; i<1000; ++i){
    size_t val = dist(gen);
    hashes.push_back(val);
    add(filter, val, 3);
  }

  for(size_t i=0; i<1000; ++i){
    if(dist2(gen)){
      checks.push_back(hashes[i]);
    }else{
      checks.push_back(dist(gen));
    }
  }
  
 
  for(auto e: checks){
    bool in_hashes = std::find(hashes.begin(), hashes.end(), e) != hashes.end();
    bool in_filter = contains(filter, e, 3);

    if(in_filter and in_hashes){
      println("True positive: {}", e);
    }else if( !in_filter and in_hashes){
      println("False negative: {}", e);
    }else if(in_filter and !in_hashes){
      println("False positive: {}", e);
    }else{
      println("True negative: {}", e);
    } 
  }
}

int main(){

  check(); 

}
