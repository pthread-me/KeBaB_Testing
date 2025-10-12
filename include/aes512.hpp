#include <smmintrin.h>
#include <tmmintrin.h>
#include <immintrin.h>
#include <cstdint>
#include <cmath>
#include <array>

class AesHash{
	public:		
		static constexpr uint8_t rounds = 10;
		static constexpr uint8_t sets = 4;

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

			seed  = _mm_set_epi64x((long long)seeds[2], (long long)seeds[3]); 
			m128i_keyset[10] = seed;
			m128i_keyset[11] = _mm_aeskeygenassist_si128(seed, 1);
			m128i_keyset[12] = _mm_aeskeygenassist_si128(seed, 2);
			m128i_keyset[13] = _mm_aeskeygenassist_si128(seed, 3);
			m128i_keyset[14] = _mm_aeskeygenassist_si128(seed, 4);
			m128i_keyset[15] = _mm_aeskeygenassist_si128(seed, 5);
			m128i_keyset[16] = _mm_aeskeygenassist_si128(seed, 6);
			m128i_keyset[17] = _mm_aeskeygenassist_si128(seed, 7);
			m128i_keyset[18] = _mm_aeskeygenassist_si128(seed, 8);
			m128i_keyset[19] = _mm_aeskeygenassist_si128(seed, 9);

			seed  = _mm_set_epi64x((long long)seeds[2], (long long)seeds[3]); 
			m128i_keyset[20] = seed;
			m128i_keyset[21] = _mm_aeskeygenassist_si128(seed, 1);
			m128i_keyset[22] = _mm_aeskeygenassist_si128(seed, 2);
			m128i_keyset[23] = _mm_aeskeygenassist_si128(seed, 3);
			m128i_keyset[24] = _mm_aeskeygenassist_si128(seed, 4);
			m128i_keyset[25] = _mm_aeskeygenassist_si128(seed, 5);
			m128i_keyset[26] = _mm_aeskeygenassist_si128(seed, 6);
			m128i_keyset[27] = _mm_aeskeygenassist_si128(seed, 7);
			m128i_keyset[28] = _mm_aeskeygenassist_si128(seed, 8);
			m128i_keyset[29] = _mm_aeskeygenassist_si128(seed, 9);

			seed  = _mm_set_epi64x((long long)seeds[2], (long long)seeds[3]); 
			m128i_keyset[30] = seed;
			m128i_keyset[31] = _mm_aeskeygenassist_si128(seed, 1);
			m128i_keyset[32] = _mm_aeskeygenassist_si128(seed, 2);
			m128i_keyset[33] = _mm_aeskeygenassist_si128(seed, 3);
			m128i_keyset[34] = _mm_aeskeygenassist_si128(seed, 4);
			m128i_keyset[35] = _mm_aeskeygenassist_si128(seed, 5);
			m128i_keyset[36] = _mm_aeskeygenassist_si128(seed, 6);
			m128i_keyset[37] = _mm_aeskeygenassist_si128(seed, 7);
			m128i_keyset[38] = _mm_aeskeygenassist_si128(seed, 8);
			m128i_keyset[39] = _mm_aeskeygenassist_si128(seed, 9);


			// aligning keys then loading them to the 512 bits expected by aes
			for(size_t i=0; i<rounds; ++i){
				__m128i hash_keys[4];
				hash_keys[0] = m128i_keyset[i];
				hash_keys[1] = m128i_keyset[i+10];
				hash_keys[1] = m128i_keyset[i+20];
				hash_keys[1] = m128i_keyset[i+30];
				
				round_keys[i] = _mm512_loadu_si512((void*) hash_keys);
			}
		}

		inline uint64_t operator()(uint64_t mer, const std::array<__m512i, rounds> &round_keys) const{
			uint64_t padded_mer[8] = {mer, mer, mer, mer, mer, mer, mer, mer};
			__m512i hash = _mm512_loadu_si512((void*) padded_mer);

			for(size_t i=0; i<rounds; ++i){
				hash = _mm512_aesenc_epi128(hash, round_keys[i]);
			}

			return _mm_extract_epi64(
					_mm512_extracti64x2_epi64(hash, 0), 0);
		}

		inline uint64_t operator()(uint64_t mer) const{
			return this->operator()(mer, round_keys);
		}
		inline uint64_t operator()(uint64_t mer, uint64_t _seed) const {
			return this->operator()(mer, round_keys);
		}
	
	private:
		std::array<__m128i, rounds*sets> m128i_keyset; 
		std::array<__m512i, rounds> round_keys; 

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




template<typename Hash=AesHash, typename Reducer=ShiftReducer>
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

    uint64_t reduce(uint64_t hash) const {
        return reducer_(hash);
    }

private:
    Hash hash_;
    Reducer reducer_;
};


using AesShift = DomainHashFunction<AesHash, ShiftReducer>;

