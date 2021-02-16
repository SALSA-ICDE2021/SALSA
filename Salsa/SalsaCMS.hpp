#pragma once

#ifndef COUNT_MIN_AEE_SALSA
#define COUNT_MIN_AEE_SALSA

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
#include "topK.hpp"

using namespace std;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class SalsaCMSSanity {

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	uint8_t **baseline_cms_counters;

	// aliases
	uint16_t **baseline_cms_counters_16;
	uint32_t **baseline_cms_counters_32;
	uint64_t **baseline_cms_counters_64;

	uint16_t **baseline_cms_merges;

	uint16_t merges_lookup_table_8[16];
	uint16_t merges_lookup_table_16[16];
	uint16_t merges_lookup_table_32[16];

	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	void merge_8_bit_counters(uint index, uint row);
	void merge_16_bit_counters(uint index, uint row);
	void merge_32_bit_counters(uint index, uint row);

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	uint8_t largest_log_counter;

	uint8_t downsamplings_per_merge;
	uint8_t remaining_downsamplings_per_merge;

	int* counter_indices;
	uint8_t* log_counter_sizes;
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

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////

public:

	SalsaCMSSanity(int width, int height, int seed, int downsamplings_per_merge=3);
	~SalsaCMSSanity();

	void increment(const char * str);
	uint64_t query(const char * str);

};

class FineGrainedSalsaCMSSanity {

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

	//uint16_t merges_lookup_table_8[16];
	//uint16_t merges_lookup_table_16[16];
	//uint16_t merges_lookup_table_32[16];

	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	void set_non_overflow_counter(uint index, uint row, uint8_t log_counter_size_value, uint64_t new_value, int bitSize);

	void merge_size_1_counters(uint index, uint row);
	void merge_size_2_counters(uint index, uint row);
	void merge_size_4_counters(uint index, uint row);
	void merge_size_8_counters(uint index, uint row);
	void merge_size_16_counters(uint index, uint row);
	void merge_size_32_counters(uint index, uint row);

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	uint8_t largest_log_counter;

	uint8_t downsamplings_per_merge;
	uint8_t remaining_downsamplings_per_merge;

	int* counter_indices;
	uint8_t* log_counter_sizes;
	uint64_t* counter_values;
	uint8_t* counter_bitSizes;

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

	void split_counters();
	void split_counter(int row, int index, uint8_t log_counter_size_value, uint64_t counter_value);

	bool split_counters_flag;
	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////

public:

	FineGrainedSalsaCMSSanity(int counterSize, int width, int height, int seed, int downsamplings_per_merge, bool split_counters_flag);
	~FineGrainedSalsaCMSSanity();

	// for debug
	inline uint64_t read_row_value(uint index, uint row)
	{
		return read_counter(index, row);
	}
	inline uint64_t read_row_value_ft(const char * str, uint row)
	{
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
		return read_counter(index, row);
	}

	void increment(const char * str);
	uint64_t query(const char * str);

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class SalsaAnalyticalErrorCMS {

	double epsilonSketch;
	double delta;

	int downsamplings_before_merge;

	int numPackets;

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	uint8_t **baseline_cms_counters;

	// aliases
	uint16_t **baseline_cms_counters_16;
	uint32_t **baseline_cms_counters_32;
	uint64_t **baseline_cms_counters_64;

	uint16_t **baseline_cms_merges;

	uint16_t merges_lookup_table_8[16];
	uint16_t merges_lookup_table_16[16];
	uint16_t merges_lookup_table_32[16];

	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	void merge_8_bit_counters(uint index, uint row);
	void merge_16_bit_counters(uint index, uint row);
	void merge_32_bit_counters(uint index, uint row);

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	uint8_t largest_log_counter;

	int* counter_indices;
	uint8_t* log_counter_sizes;
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

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////

public:

	SalsaAnalyticalErrorCMS(int width, int height, int seed, double delta, int downsamplings_before_merge);
	~SalsaAnalyticalErrorCMS();

	void increment(const char * str);
	uint64_t query(const char * str);

};

class FineGrainedSalsaAnalyticalErrorCMSSanity {

	double epsilonSketch;
	double delta;

	int numPackets;

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

	//uint16_t merges_lookup_table_8[16];
	//uint16_t merges_lookup_table_16[16];
	//uint16_t merges_lookup_table_32[16];

	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	void set_non_overflow_counter(uint index, uint row, uint8_t log_counter_size_value, uint64_t new_value, int bitSize);

	void merge_size_1_counters(uint index, uint row);
	void merge_size_2_counters(uint index, uint row);
	void merge_size_4_counters(uint index, uint row);
	void merge_size_8_counters(uint index, uint row);
	void merge_size_16_counters(uint index, uint row);
	void merge_size_32_counters(uint index, uint row);

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////
	uint8_t largest_log_counter;

	int downsamplings_before_merge;

	int* counter_indices;
	uint8_t* log_counter_sizes;
	uint64_t* counter_values;
	uint8_t* counter_bitSizes;

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

	void split_counters();
	void split_counter(int row, int index, uint8_t log_counter_size_value, uint64_t counter_value);

	bool split_counters_flag;
	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////

public:

	FineGrainedSalsaAnalyticalErrorCMSSanity(int counterSize, int width, int height, int seed, int downsamplings_before_merge, bool split_counters_flag, double delta);
	~FineGrainedSalsaAnalyticalErrorCMSSanity();

	// for debug
	inline uint64_t read_row_value(uint index, uint row)
	{
		return read_counter(index, row);
	}
	inline uint64_t read_row_value_ft(const char * str, uint row)
	{
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
		return read_counter(index, row);
	}

	void increment(const char * str);
	uint64_t query(const char * str);

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#endif // COUNT_MIN_AEE_SALSA
