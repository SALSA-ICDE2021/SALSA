#pragma once

#ifndef COUNT_MIN_CC
#define COUNT_MIN_CC

#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>

#include "RngFast.hpp"
#include "BobHash.hpp"
#include "AEE_Defs.hpp"

using namespace std;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class CountMinCC {

	int seed;

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	cc_counter_t **cc_cms;

	int* indices;

	mt19937 generator_arr;
	uniform_real_distribution<> dis_arr;

	rng::tsc_seed seed_arr;
	rng::rng128 gen_arr;

	uint32_t w, w1, w2;

	uint32_t width_over_four;

	uint8_t minus_log_p = 0;
	uint8_t minus_log_p_bound = 6;

	double LogOneMinusP;

	uint32_t interval;

	uint64_t curr_sum;

	uint8_t rand_pos = 0;
	__m256i rand_bits;
	char * RandByteArray;

	uint8_t weighted_rand_pos = 0;
	__m256i weighted_rand_bits;
	uint16_t * WeightedRandByteArray;

	/* mask for parallel devision */
	uint64_t division_mask = 0xFFFEFFFEFFFEFFFE;

	inline unsigned accelerated_geometric(double LogOneMinusP);

	void divide_array_counters();

	void increment_h(const char * str);

	void add_h(const char * str, int c);

public:

	CountMinCC(int width, int height, int seed);
	~CountMinCC();

	// Unweighted
	void increment(const char * str);

	// Weighted
	void add(const char * str, int c);

	uint64_t query(const char * str);

};

class CountMinCCSanity {

	uint64_t numPackets;

	int seed;

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	cc_counter_t **cc_cms;

	int* indices;

	mt19937 generator_arr;
	uniform_real_distribution<> dis_arr;

	rng::tsc_seed seed_arr;
	rng::rng128 gen_arr;

	uint32_t w, w1, w2;

	uint32_t width_over_four;

	uint8_t minus_log_p = 0;
	uint8_t minus_log_p_bound = 6;

	double LogOneMinusP;

	uint32_t interval;

	uint64_t curr_sum;

	uint8_t rand_pos = 0;
	__m256i rand_bits;
	char * RandByteArray;

	uint8_t weighted_rand_pos = 0;
	__m256i weighted_rand_bits;
	uint16_t * WeightedRandByteArray;

	/* mask for parallel devision */
	uint64_t division_mask = 0xFFFEFFFEFFFEFFFE;

	inline unsigned accelerated_geometric(double LogOneMinusP);

	void divide_array_counters();

	void increment_h(const char * str);

	void add_h(const char * str, int c);

public:

	CountMinCCSanity(int width, int height, int seed);
	~CountMinCCSanity();

	// Unweighted
	void increment(const char * str);

	// Weighted
	void add(const char * str, int c);

	uint64_t query(const char * str);

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class ConservativeUpdateCC {

	int seed;

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	cc_counter_t **cc_cus;

	int* indices;

	mt19937 generator_arr;
	uniform_real_distribution<> dis_arr;

	rng::tsc_seed seed_arr;
	rng::rng128 gen_arr;

	uint32_t w, w1, w2;

	uint32_t width_over_four;

	uint8_t minus_log_p = 0;
	uint8_t minus_log_p_bound = 6;

	double LogOneMinusP;

	uint32_t interval;

	uint64_t curr_sum;

	uint8_t rand_pos = 0;
	__m256i rand_bits;
	char * RandByteArray;

	/* mask for parallel devision */
	uint64_t division_mask = 0xFFFEFFFEFFFEFFFE;

	inline unsigned accelerated_geometric(double LogOneMinusP);

	void divide_array_counters();

	void increment_h(const char * str);

public:

	ConservativeUpdateCC(int width, int height, int seed);
	~ConservativeUpdateCC();

	void increment(const char * str);
	uint64_t query(const char * str);

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class CountMinCC_MaxSpeed {

	uint64_t N_current, N_prime, Shifted_N_prime;

	int seed;

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	cc_counter_t **cc_cms;

	int* indices;

	mt19937 generator_arr;
	uniform_real_distribution<> dis_arr;

	rng::tsc_seed seed_arr;
	rng::rng128 gen_arr;

	uint32_t w, w1, w2;

	uint32_t width_over_four;

	uint8_t minus_log_p = 0;
	uint8_t minus_log_p_bound = 6;

	double LogOneMinusP;

	uint32_t interval;

	uint64_t curr_sum;

	uint8_t rand_pos = 0;
	__m256i rand_bits;
	char * RandByteArray;

	uint8_t weighted_rand_pos = 0;
	__m256i weighted_rand_bits;
	uint16_t * WeightedRandByteArray;

	/* mask for parallel devision */
	uint64_t division_mask = 0xFFFEFFFEFFFEFFFE;

	inline unsigned accelerated_geometric(double LogOneMinusP);

	void divide_array_counters();

	void increment_h(const char * str);

public:

	CountMinCC_MaxSpeed(int width, int height, int seed, double epsilon_cc, double delta_cc);
	~CountMinCC_MaxSpeed();

	void increment(const char * str);

	uint64_t query(const char * str);

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#endif