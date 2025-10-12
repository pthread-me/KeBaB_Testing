#include <aes.hpp>
#include <array>
#include <nt_hash.hpp>
#include <fstream>
#include <string>
#include <print>
#include <chrono>

using namespace std;
using namespace kebab;

int main(){
	fstream input("/home/Cutie/hub/datasets/kebab/Ciona_savignyi.CSAV2.0.dna_sm.toplevel.fa");	

	{
		string discard;
		getline(input, discard);
		print("line1: {}\n", discard);
	}

	string line;
	string dataset;
	while(getline(input, line)) dataset.append(line);
	
	
	size_t kmer_size = 5;	
	NtHash<uint64_t> mul_nt(kmer_size);
	NtHash<uint64_t> acc_nt(kmer_size);
	NtHash<uint64_t> man_nt(kmer_size);
	NtHash<uint64_t> mur_nt(kmer_size);
	mul_nt.set_sequence(dataset.c_str(), dataset.length());
	acc_nt.set_sequence(dataset.c_str(), dataset.length());
	man_nt.set_sequence(dataset.c_str(), dataset.length());
	mur_nt.set_sequence(dataset.c_str(), dataset.length());

	AesShift aes;	
	MulShift mul;
	NtManyHash nt;
	MurmurHash2 mur;

	uint64_t seed = 0xc6a4a7935bd1e995LL;
	uint64_t seed2 = 0xc6a4a7944bd1e995LL;
	uint64_t seed3 = 0x36a4a7944bd1e995LL;
	uint64_t seed4 = 0x44a4a7944bd1e995LL;
	
	uint64_t mul_sum = 0;
	uint64_t aes_sum = 0;
	uint64_t man_sum = 0;
	uint64_t mur_sum = 0;
		
	auto mul_start = chrono::high_resolution_clock::now();
	while(true){
		uint64_t h = mul(mul_nt.hash(), seed);
		uint64_t h2 = mul(mul_nt.hash(), seed2);
		uint64_t h3 = mul(mul_nt.hash(), seed3);
		uint64_t h4 = mul(mul_nt.hash(), seed4);
		mul_sum += (h+h2+h3+h4);
		if(!mul_nt.roll()) break;
	}
	auto mul_end = chrono::high_resolution_clock::now();
	
	
	auto aes_start = chrono::high_resolution_clock::now();
	while(true){
		std::array<uint64_t, 2> h =  aes.acc_hash(acc_nt.hash(), seed);
		aes_sum += (h[0] + h[1]);
		if(!acc_nt.roll()) break;
	}
	auto aes_end = chrono::high_resolution_clock::now();

	auto man_start = chrono::high_resolution_clock::now();
	while(true){
		uint64_t h = nt(man_nt.hash(), seed);
		uint64_t h2 = nt(man_nt.hash(), seed2);
		uint64_t h3 = nt(man_nt.hash(), seed3);
		uint64_t h4 = nt(man_nt.hash(), seed4);
		mul_sum += (h+h2+h3+h4);
		if(!man_nt.roll()) break;
	}
	auto man_end = chrono::high_resolution_clock::now();
	

	auto mur_start = chrono::high_resolution_clock::now();
	while(true){
		uint64_t h = mur(mur_nt.hash(), seed);
		uint64_t h2 = mur(mur_nt.hash(), seed2);
		uint64_t h3 = mur(mur_nt.hash(), seed3);
		uint64_t h4 = mur(mur_nt.hash(), seed4);
		mul_sum += (h+h2+h3+h4);
		if(!mur_nt.roll()) break;
	}
	auto mur_end = chrono::high_resolution_clock::now();
	

	println("mul time {}", chrono::duration_cast<chrono::milliseconds>(mul_end - mul_start));
	println("aes time {}", chrono::duration_cast<chrono::milliseconds>(aes_end - aes_start));
	println("many nt time {}", chrono::duration_cast<chrono::milliseconds>(man_end - man_start));
	println("mur time {}", chrono::duration_cast<chrono::milliseconds>(mur_end - mur_start));

	println("sums: {} {} {} {}", mul_sum, aes_sum, man_sum, mur_sum);
}
