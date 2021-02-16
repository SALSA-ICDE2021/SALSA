#pragma once

#ifndef COUNT_MIN_TANGO
#define COUNT_MIN_TANGO

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

class TangoBaseline {

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

public:

	TangoBaseline();
	~TangoBaseline();

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

	void initialize(int counterSize, int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

class TangoBaselineSanity {

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

public:

	TangoBaselineSanity();
	~TangoBaselineSanity();

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

	void initialize(int counterSize, int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

class smarTangoBaselineSanity {

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

public:

	smarTangoBaselineSanity();
	~smarTangoBaselineSanity();

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

	void initialize(int counterSize, int width, int height, int seed);
	void increment(const char * str);
	uint64_t query(const char * str);

};

#endif // COUNT_MIN_SALSA
