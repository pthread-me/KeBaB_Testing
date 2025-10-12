#include <algorithm>
#include <pmmintrin.h>
#include <smmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>
#include <cstdint>
#include <cmath>
#include <array>
#include <wmmintrin.h>

class MurmurHash2 {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const {
        uint64_t h = seed ^ (len * m);
        uint64_t k = x;
        
        k *= m;
        k ^= k >> r;
        k *= m;
        
        h ^= k;
        h *= m;
        
        h ^= h >> r;
        h *= m;
        h ^= h >> r;
        
        return h;
    }   
private:
    static constexpr uint64_t m = 0xc6a4a7935bd1e995ULL;
    static constexpr uint8_t r = 47;
    static constexpr uint8_t len = 8;
};



class NtManyHash {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const {
        x *= seed;
        x ^= x >> shift;
        return x;
    }

    uint64_t operator()(uint64_t x) const {
        return operator()(x, seed);
    }
private:
    static constexpr uint8_t shift = 27;
    static constexpr uint64_t seed = 0x90b45d39fb6da1fa;
};

class MultiplyHash {
public:
    uint64_t operator()(uint64_t x, uint64_t seed) const {
        return x * seed;
    }
		std::array<uint64_t, 2> acc_hash(uint64_t x, uint64_t seed){
			return std::array<uint64_t, 2>();
		}
};



class AesHash{
	public:		
		static constexpr uint8_t rounds = 3;
		static constexpr uint8_t sets = 1;

		AesHash(){
			__m128i seed;
			
			// key calculation for rounds requires constant int, thus is hardcoded 
			// for each aes hashing operation (4 op for 512 bits)
			seed  = _mm_set_epi64x((long long)seeds[0], (long long)seeds[1]); 
			m128i_keyset[0] = seed;
			m128i_keyset[1] = _mm_aeskeygenassist_si128(seed, 1);
			m128i_keyset[2] = _mm_aeskeygenassist_si128(seed, 2);
			m128i_keyset[3] = _mm_aeskeygenassist_si128(seed, 3);
			m128i_keyset[4] = _mm_aeskeygenassist_si128(seed, 4);
			m128i_keyset[5] = _mm_aeskeygenassist_si128(seed, 5);
			m128i_keyset[6] = _mm_aeskeygenassist_si128(seed, 6);
			m128i_keyset[7] = _mm_aeskeygenassist_si128(seed, 7);
			m128i_keyset[8] = _mm_aeskeygenassist_si128(seed, 8);
			m128i_keyset[9] = _mm_aeskeygenassist_si128(seed, 9);
		}



		inline std::array<uint64_t, 2> acc_hash(uint64_t mers, const std::array<__m128i, rounds> &round_keys) const{
			uint64_t padded_mer[2] = {mers, mers};
			__m128i hash = _mm_lddqu_si128((__m128i*)padded_mer);

			for(size_t i=0; i<rounds; ++i) hash = _mm_aesenc_si128(hash, round_keys[i]);

			return std::array<uint64_t, 2>{_mm_extract_epi64(hash, 0), _mm_extract_epi64(hash, 1)};
		}


		inline uint64_t operator()(uint64_t mer, const std::array<__m128i, rounds> &round_keys) const{
			uint64_t padded_mer[2] = {mer, mer};
			__m128i hash = _mm_lddqu_si128((__m128i*)padded_mer);
			for(size_t i=0; i<rounds; ++i) hash = _mm_aesenc_si128(hash, round_keys[i]);
			return _mm_extract_epi64(hash, 0);
		}


		inline uint64_t operator()(uint64_t mer) const{
			return this->operator()(mer, round_keys);
		}
		inline uint64_t operator()(uint64_t mer, uint64_t seed) const{
			return this->operator()(mer, round_keys);
		}
		inline std::array<uint64_t, sets*2> acc_hash(uint64_t mer, uint64_t seed) const {
			return acc_hash(mer, round_keys);
		}

	
	private:
		// useful for the 512 bit implementation
		std::array<__m128i, rounds*sets> m128i_keyset; 
		std::array<__m128i, rounds> round_keys; 

		static constexpr std::array<uint64_t, 8> seeds{
			0xc6a4a7935bd1e995LL, 0x90b45d39fb6da1faLL,
			0xc914aab35bd1e995LL, 0x80c45a37fb6dadfaLL,
			0xce9a914ab35195bdLL, 0x37fadf80c5ab6da4LL,
			0xc7f6aadca67facadLL, 0x67facadeababab00LL,
		};
};


class ShiftReducer {
public:
    explicit ShiftReducer(size_t domain_size) 
        : shift(64 - std::floor(std::log2(domain_size))) {}
    
    uint64_t operator()(uint64_t hash) const {
        return hash >> shift;
    }
private:
    uint8_t shift;
};




template<typename Hash, typename Reducer>
class DomainHashFunction {
public:
    DomainHashFunction()
        : hash_(Hash())
        , reducer_(Reducer(0)) {}

    DomainHashFunction(size_t domain_size)
        : hash_(Hash())
        , reducer_(Reducer(domain_size)) {}

    uint64_t operator()(uint64_t x, uint64_t seed) const {
        return reducer_(hash_(x, seed));
    }
		
    uint64_t hash(uint64_t x, uint64_t seed) const {
        return hash_(x, seed);
    }
		
		std::array<uint64_t, 2> acc_hash(uint64_t x, uint64_t seed){
			std::array<uint64_t, 2> h = hash_.acc_hash(x, seed);
			std::array<uint64_t, 2> res;
			std::transform(h.begin(), h.end(), res.begin(), reducer_);
			return res;
		}
    uint64_t reduce(uint64_t hash) const {
        return reducer_(hash);
    }

private:
    Hash hash_;
    Reducer reducer_;
};

using AesShift = DomainHashFunction<AesHash, ShiftReducer>;
using MulShift = DomainHashFunction<MultiplyHash, ShiftReducer>;
using ManShift = DomainHashFunction<NtManyHash, ShiftReducer>;
using MurShift = DomainHashFunction<MurmurHash2, ShiftReducer>;
