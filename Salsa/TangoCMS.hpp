#pragma once

#ifndef COUNT_MIN_TANGO_AEE
#define COUNT_MIN_TANGO_AEE

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
#include "Defs.hpp"

using namespace std;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class smarTango {

	int counterSize;
	int logCounterSize;
	int counters_per_word;
	int log_counters_per_word;
	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	uint64_t **baseline_cms_counters;
	uint64_t **baseline_cms_merges;

	inline pair<uint8_t, uint8_t> counter_size_and_offset(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	bool isBitSet(uint64_t merges_word, uint8_t bitIndex)
	{
		return !!(merges_word & (((uint64_t)1) << bitIndex));
	}

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	uint8_t largest_counter;

	uint8_t downsamplings_per_merge;
	uint8_t remaining_downsamplings_per_merge;

	int* counter_indices;
	uint8_t* counter_sizes;
	uint8_t* counter_offsets;
	uint64_t* counter_values;

	mt19937 generator_arr;
	uniform_real_distribution<> dis_arr;

	rng::tsc_seed seed_arr;
	rng::rng128 gen_arr;

	uint8_t minus_log_p = 0;
	uint8_t minus_log_p_bound = 6;

	double LogOneMinusP;

	uint32_t interval;

	uint8_t rand_pos = 0;
	__m256i rand_bits;
	char * RandByteArray;

	inline unsigned accelerated_geometric(double LogOneMinusP) {
		// Return a randomly drawn geometric variable with parameter p, such that LogOneMinusP = log(1-p)
		return log(dis_arr(generator_arr)) / LogOneMinusP + 1;
	}

	void increment_h(const char * str);
	void divide_counters();
	void handle_overflow_and_increment();
	void increment_and_potentially_merge();

	uint8_t what_would_be_the_merged_counter_size(uint row);

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////

public:

	smarTango(int counterSize, int width, int height, int seed, int downsamplings_per_merge = 3);
	~smarTango();

	// for debug
	inline uint64_t read_row_value(uint index, uint row)
	{
		return read_counter(index, row);
	}
	inline uint64_t read_row_value_ft(const char * str, uint row)
	{
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		return read_counter(index, row);
	}

	void increment(const char * str);
	uint64_t query(const char * str);

};

class Tango {

	int counterSize;
	int logCounterSize;
	int counters_per_word;
	int log_counters_per_word;
	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	uint64_t **baseline_cms_counters;
	uint64_t **baseline_cms_merges;

	inline pair<uint8_t, uint8_t> counter_size_and_offset(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	bool isBitSet(uint64_t merges_word, uint8_t bitIndex)
	{
		return !!(merges_word & (((uint64_t)1) << bitIndex));
	}

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	uint8_t largest_counter;

	uint8_t downsamplings_per_merge;
	uint8_t remaining_downsamplings_per_merge;

	int* counter_indices;
	uint8_t* counter_sizes;
	uint8_t* counter_offsets;
	uint64_t* counter_values;

	mt19937 generator_arr;
	uniform_real_distribution<> dis_arr;

	rng::tsc_seed seed_arr;
	rng::rng128 gen_arr;

	uint8_t minus_log_p = 0;
	uint8_t minus_log_p_bound = 6;

	double LogOneMinusP;

	uint32_t interval;

	uint8_t rand_pos = 0;
	__m256i rand_bits;
	char * RandByteArray;

	inline unsigned accelerated_geometric(double LogOneMinusP) {
		// Return a randomly drawn geometric variable with parameter p, such that LogOneMinusP = log(1-p)
		return log(dis_arr(generator_arr)) / LogOneMinusP + 1;
	}

	void increment_h(const char * str);
	void divide_counters();
	void handle_overflow_and_increment();
	void increment_and_potentially_merge();

	uint8_t what_would_be_the_merged_counter_size(uint row);

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////

public:

	Tango(int counterSize, int width, int height, int seed, int downsamplings_per_merge = 3);
	~Tango();

	// for debug
	inline uint64_t read_row_value(uint index, uint row)
	{
		return read_counter(index, row);
	}
	inline uint64_t read_row_value_ft(const char * str, uint row)
	{
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		return read_counter(index, row);
	}

	void increment(const char * str);
	uint64_t query(const char * str);

};

#endif // COUNT_MIN_TANGO_AEE
