#pragma once

#ifndef COUNT_SKETCH
#define COUNT_SKETCH

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
#include "RngFast.hpp"
#include "topK.hpp"

using namespace std;

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class CountSketchBaseline {

protected:

	int width;
	int height;

	int width_and_sign_mask;

	BOBHash *bobhash;

	int32_t **counters;

	uint64_t *counter_values;

	virtual uint64_t MedianCounterValue()
	{
		sort(counter_values, counter_values + height);
		return counter_values[height >> 1];
	}

public:

	CountSketchBaseline(int width, int height, int seed);
	~CountSketchBaseline();

	virtual void increment(const char * str);
	virtual uint64_t query(const char * str);

	// for debug
	inline int64_t read_row_value(uint index, uint row)
	{
		return counters[row][index];
	}
	inline int64_t read_row_value_ft(const char * str, uint row)
	{
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_and_sign_mask) >> 1;
		return counters[row][index];
	}

};

class CountSketchBaselineFiveRows : public CountSketchBaseline {

	uint64_t MedianCounterValue();

public:

	CountSketchBaselineFiveRows(int width, int seed);
	~CountSketchBaselineFiveRows();

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

class FineGrainedSalsaCountSketch {

	int inputSeed;

	int counterSize;
	int logCounterSize;
	int counters_per_word;
	int log_counters_per_word;

	bool maximum_merge;

	int width;
	int height;

	int width_and_sign_mask;

	BOBHash *bobhash;

	uint64_t **counters;
	uint64_t **merges_and_signs;
	uint32_t **merges_and_signs_32;

	
	inline int64_t read_counter(uint index, uint row);
	inline int64_t compute_counter_pair_sum(uint index, uint row, uint8_t log_counter_size_value);
	inline uint8_t log_counter_size(uint index, uint row);

	void set_non_overflow_counter(uint index, uint row, uint8_t log_counter_size_value, int64_t new_value);
	void merge_counters(uint index, uint row, uint8_t log_counter_size_value);

	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////

	uint8_t largest_log_counter;

	uint8_t downsamplings_per_merge;
	uint8_t remaining_downsamplings_per_merge;

	int* counter_indices;
	int* flow_signs;
	uint8_t* log_counter_sizes;
	int64_t* counter_values;
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

	bool split_counters_flag;
	///////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////

	uint64_t MedianCounterValue()
	{
		sort(counter_values, counter_values + height);
		return (uint64_t)counter_values[height >> 1];
	}

public:

	FineGrainedSalsaCountSketch(int counterSize, int width, int height, int seed, int downsamplings_per_merge, bool maximum_merge);
	~FineGrainedSalsaCountSketch();

	// for debug
	inline int64_t read_row_value(uint index, uint row)
	{
		return read_counter(index, row);
	}
	inline int64_t read_log_counter_size_value(uint index, uint row)
	{
		return log_counter_size(index, row);
	}
	inline int64_t read_row_value_ft(const char * str, uint row)
	{
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_and_sign_mask) >> 1;
		return read_counter(index, row);
	}

	virtual void increment(const char * str);
	uint64_t query(const char * str);

	long double l2_estimation();

	////

	int64_t signedQuery(const char * str);

	FineGrainedSalsaCountSketch* operator - (FineGrainedSalsaCountSketch& subtrahend);

};

class FineGrainedSalsaCountSketchHH : public FineGrainedSalsaCountSketch {

	//mapTopK< string, uint64_t> topK;
	orderedMapTopK< string, uint64_t> topK;
public:

	FineGrainedSalsaCountSketchHH(int counterSize, int width, int height, int seed, int downsamplings_per_merge, int heapSize, bool maximum_merge);
	~FineGrainedSalsaCountSketchHH();

	void increment(const char * str);

	vector<pair<string, uint64_t>> HH();

};

////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////

#endif // 
