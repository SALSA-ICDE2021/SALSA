#pragma once

#ifndef COUNT_MIN_SALSA
#define COUNT_MIN_SALSA

#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>

#include "BobHash.hpp"
#include "Defs.hpp"

using namespace std;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class SalsaCMSBaseline {

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
	void merge_64_bit_counters(uint index, uint row);

public:

	SalsaCMSBaseline();
	~SalsaCMSBaseline();

	void initialize(int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

class MaximumSalsaCMSBaseline {

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

public:

	MaximumSalsaCMSBaseline();
	~MaximumSalsaCMSBaseline();

	void initialize(int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

class MaximumSalsaCMSBaselineSanity {

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

public:

	MaximumSalsaCMSBaselineSanity();
	~MaximumSalsaCMSBaselineSanity();

	// for debug
	inline uint64_t read_row_value(uint index, uint row)
	{
		return read_counter(index, row);
	}
	inline uint64_t read_row_value_ft(const char * str, uint row)
	{
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> 3;
		return read_counter(index, row);
	}

	void initialize(int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

class FineGrainedSalsaCMSBaseline {

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

public:

	FineGrainedSalsaCMSBaseline(int counterSize, int width, int height, int seed);
	~FineGrainedSalsaCMSBaseline();

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

class FineGrainedSalsaCMSBaselineSanity {

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

	// for count distinct 
	int mergedZeroValuedCounters[6];

	int **mergedZeroValuedCountersMatrix;

	// row count distinct estimate
	double countDistinctRowEstimate(uint row);
	double countDistinctRowEstimateAverage(uint row);

public:

	FineGrainedSalsaCMSBaselineSanity(int counterSize, int width, int height, int seed);
	~FineGrainedSalsaCMSBaselineSanity();

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

	double provable_query_count_distinct();
	double bound_average_query_count_distinct();
	double global_query_count_distinct();
	double best_guess_query_count_distinct();
	double global_best_guess_query_count_distinct();

	void increment(const char * str);
	uint64_t query(const char * str);

};


////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class MaximumSalsaCUSBaseline {

	int width;
	int height;

	int width_mask;

	BOBHash *bobhash;

	uint8_t **baseline_cus_counters;

	// aliases
	uint16_t **baseline_cus_counters_16;
	uint32_t **baseline_cus_counters_32;
	uint64_t **baseline_cus_counters_64;

	uint16_t **baseline_cus_merges;

	uint16_t merges_lookup_table_8[16];
	uint16_t merges_lookup_table_16[16];
	uint16_t merges_lookup_table_32[16];

	uint* incremet_indices;
	uint64_t* incremet_values;

	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	void merge_8_bit_counters(uint index, uint row);
	void merge_16_bit_counters(uint index, uint row);
	void merge_32_bit_counters(uint index, uint row);

public:

	MaximumSalsaCUSBaseline();
	~MaximumSalsaCUSBaseline();

	void initialize(int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

#endif // COUNT_MIN_SALSA
