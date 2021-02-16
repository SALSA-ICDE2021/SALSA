#pragma once
#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string.h>
#include <iostream>
#include <iomanip>

#include "BobHash.hpp"
#include "Defs.hpp"


#ifndef SMART_SALSA_H
#define SMART_SALSA_H


/// <summary>
/// ///
/// </summary>

class SmartEncodingSALSA {

	int width;
	int height;

	int width_mask;
	int log_width;


#ifdef USE_BOBHASH
	BOBHash* bobhash;
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
	int seeds[12];
#endif // if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES	

#ifdef USE_TWO_XXHASHES
	int height_plus_1_over_2;
#endif // USE_TWO_XXHASHES	

	uint8_t** basic_counters;

	// aliases
	uint16_t** size_16_counter_alias;
	uint32_t** size_32_counter_alias;

	uint8_t** twelve_counter_encodings;


	uint8_t mod12_lookup[512];

	//uint8_t smart_encoding_lookup_table[125 * 12];

	//uint8_t smart_merging_lookup_table[125 * 12];

	uint8_t smart_encoding_lookup_table[256 * 12];

	uint8_t smarter_merging_lookup_table_0[256 * 12];
	uint8_t smarter_merging_lookup_table_1[256 * 12];

	//uint8_t smarter_merging_lookup_table_0[125 * 12];
	//uint8_t smarter_merging_lookup_table_1[125 * 12];


	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	//void merge_8_bit_counters(uint index, uint row);
	//void merge_16_bit_counters(uint index, uint row);
	//void merge_32_bit_counters(uint index, uint row);

public:

	SmartEncodingSALSA();
	~SmartEncodingSALSA();

	void initialize(int width, int height, int seed);
	void increment(const char* str);
	uint64_t query(const char* str);

};

/////////

class FastEncodingSALSA {

	int width;
	int height;

	int width_mask;
	int log_width;


#ifdef USE_BOBHASH
	BOBHash* bobhash;
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
	int seeds[12];
#endif // if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES	

#ifdef USE_TWO_XXHASHES
	int height_plus_1_over_2;
#endif // USE_TWO_XXHASHES	

	uint8_t** basic_counters;

	// aliases
	uint16_t** size_16_counter_alias;
	uint32_t** size_32_counter_alias;

	uint8_t** eight_counter_encodings;


	//uint8_t smart_encoding_lookup_table[125 * 12];

	//uint8_t smart_merging_lookup_table[125 * 12];

	//uint8_t smarter_merging_lookup_table_0[125 * 12];
	//uint8_t smarter_merging_lookup_table_1[125 * 12];

	uint8_t smart_encoding_lookup_table[64 * 5];

	uint8_t smarter_merging_lookup_table_0[64 * 5];
	uint8_t smarter_merging_lookup_table_1[64 * 5];


	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	//void merge_8_bit_counters(uint index, uint row);
	//void merge_16_bit_counters(uint index, uint row);
	//void merge_32_bit_counters(uint index, uint row);

public:

	inline uint64_t read_counter_temp(uint index, uint row) { return read_counter(index, row); }
#ifdef USE_BOBHASH
	inline uint64_t indexof(const char* str, uint row) {return (bobhash[row].run(str, FT_SIZE))& width_mask; }
#else //defined USE_ONE_XXHASH
	inline uint64_t indexof(const char* str, uint row) { return xxh::xxhash3<64>(str, FT_SIZE, seeds[0]) & width_mask; }
#endif

	FastEncodingSALSA();
	~FastEncodingSALSA();

	void initialize(int width, int height, int seed);
	void increment(const char* str);
	uint64_t query(const char* str);

};

/////////

class StupidEncodingSALSA {

	int width;
	int height;

	int width_mask;
	int log_width;


#ifdef USE_BOBHASH
	BOBHash* bobhash;
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
	int seeds[12];
#endif // if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES	

#ifdef USE_TWO_XXHASHES
	int height_plus_1_over_2;
#endif // USE_TWO_XXHASHES	

	uint8_t** basic_counters;

	// aliases
	uint16_t** size_16_counter_alias;
	uint32_t** size_32_counter_alias;

	uint8_t** eight_counter_encodings;

	//uint8_t smart_encoding_lookup_table[125 * 12];

	//uint8_t smart_merging_lookup_table[125 * 12];

	//uint8_t smarter_merging_lookup_table_0[125 * 12];
	//uint8_t smarter_merging_lookup_table_1[125 * 12];

	uint8_t smart_encoding_lookup_table[64 * 5];

	uint8_t smarter_merging_lookup_table_0[64 * 5];
	uint8_t smarter_merging_lookup_table_1[64 * 5];


	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	//void merge_8_bit_counters(uint index, uint row);
	//void merge_16_bit_counters(uint index, uint row);
	//void merge_32_bit_counters(uint index, uint row);

public:

	inline uint64_t read_counter_temp(uint index, uint row) { return read_counter(index, row); }
#ifdef USE_BOBHASH
	inline uint64_t indexof(const char* str, uint row) { return (bobhash[row].run(str, FT_SIZE)) & width_mask; }
#elif defined USE_ONE_XXHASH
	inline uint64_t indexof(const char* str, uint row) { return xxh::xxhash3<64>(str, FT_SIZE, seeds[0]) & width_mask; }
#endif
	StupidEncodingSALSA();
	~StupidEncodingSALSA();

	void initialize(int width, int height, int seed);
	void increment(const char* str);
	uint64_t query(const char* str);

};

/////////

class StupiderEncodingSALSA {

	int width;
	int height;

	int width_mask;
	int log_width;


#ifdef USE_BOBHASH
	BOBHash* bobhash;
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
	int seeds[12];
#endif // if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES	

#ifdef USE_TWO_XXHASHES
	int height_plus_1_over_2;
#endif // USE_TWO_XXHASHES	

	uint8_t** basic_counters;

	// aliases
	uint16_t** size_16_counter_alias;
	uint32_t** size_32_counter_alias;

	uint64_t** sixty_four_counter_encodings;

	//uint8_t smart_encoding_lookup_table[125 * 12];

	//uint8_t smart_merging_lookup_table[125 * 12];

	//uint8_t smarter_merging_lookup_table_0[125 * 12];
	//uint8_t smarter_merging_lookup_table_1[125 * 12];

	uint8_t smart_encoding_lookup_table[64 * 5];

	uint8_t smarter_merging_lookup_table_0[64 * 5];
	uint8_t smarter_merging_lookup_table_1[64 * 5];


	inline uint8_t log_counter_size(uint index, uint row);
	inline uint64_t read_counter(uint index, uint row);

	//void merge_8_bit_counters(uint index, uint row);
	//void merge_16_bit_counters(uint index, uint row);
	//void merge_32_bit_counters(uint index, uint row);

public:

	inline uint64_t read_counter_temp(uint index, uint row) { return read_counter(index, row); }
#ifdef USE_BOBHASH
	inline uint64_t indexof(const char* str, uint row) { return (bobhash[row].run(str, FT_SIZE)) & width_mask; }
#elif defined USE_ONE_XXHASH
	inline uint64_t indexof(const char* str, uint row) { return xxh::xxhash3<64>(str, FT_SIZE, seeds[0]) & width_mask; }
#endif
	StupiderEncodingSALSA();
	~StupiderEncodingSALSA();

	void initialize(int width, int height, int seed);
	void increment(const char* str);
	uint64_t query(const char* str);

};

#endif // SMART_SALSA_H