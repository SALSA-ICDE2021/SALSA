#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "RngFast.hpp"
#include "CountSketch.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

CountSketchBaseline::CountSketchBaseline(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_and_sign_mask = (width << 1) - 1;

	assert((height % 2) == 1 && "We assume that height is odd");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2");
	assert(width >= 4 && "We assume that (width % 4 == 0)");
	assert(width < (1 << 30) && "We assume that you are not a database lunatic");
	
	counters = new int32_t*[height];
	bobhash = new BOBHash[height];

	counter_values = new uint64_t[height];

	for (int i = 0; i < height; ++i)
	{
		counters[i] = new int32_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}
}

CountSketchBaseline::~CountSketchBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] counters[i];
	}
	delete[] bobhash;
	delete[] counters;

	delete[] counter_values;
}

void CountSketchBaseline::increment(const char * str)
{
	for (int row = 0; row < height; ++row) {
		uint index_and_sign = (bobhash[row].run(str, FT_SIZE)) & width_and_sign_mask;
		uint index = index_and_sign >> 1;
		int sign = 1 - 2 * (index_and_sign & 0b1);
		counters[row][index] += sign;
	}
}

uint64_t CountSketchBaseline::query(const char * str)
{
	for (int row = 0; row < height; ++row) {
		uint index_and_sign = (bobhash[row].run(str, FT_SIZE)) & width_and_sign_mask;
		uint index = index_and_sign >> 1;
		int sign = 1 - 2 * (index_and_sign & 0b1);

		int64_t singned_counter_value = sign * counters[row][index];
		counter_values[row] = singned_counter_value > 0 ? singned_counter_value : 0;
	}
	return MedianCounterValue();
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

uint64_t CountSketchBaselineFiveRows::MedianCounterValue()
{
	return counter_values[1] < counter_values[0] ? counter_values[3] < counter_values[2] ? counter_values[1] < counter_values[3] ? counter_values[0] < counter_values[4] ? counter_values[0] < counter_values[3] ? counter_values[4] < counter_values[3] ? counter_values[4] : counter_values[3]
		 : counter_values[2] < counter_values[0] ? counter_values[2] : counter_values[0]
		 : counter_values[4] < counter_values[3] ? counter_values[0] < counter_values[3] ? counter_values[0] : counter_values[3]
		 : counter_values[2] < counter_values[4] ? counter_values[2] : counter_values[4]
		 : counter_values[2] < counter_values[4] ? counter_values[1] < counter_values[2] ? counter_values[0] < counter_values[2] ? counter_values[0] : counter_values[2]
		 : counter_values[4] < counter_values[1] ? counter_values[4] : counter_values[1]
		 : counter_values[1] < counter_values[4] ? counter_values[0] < counter_values[4] ? counter_values[0] : counter_values[4]
		 : counter_values[2] < counter_values[1] ? counter_values[2] : counter_values[1]
		 : counter_values[1] < counter_values[2] ? counter_values[0] < counter_values[4] ? counter_values[0] < counter_values[2] ? counter_values[4] < counter_values[2] ? counter_values[4] : counter_values[2]
		 : counter_values[3] < counter_values[0] ? counter_values[3] : counter_values[0]
		 : counter_values[4] < counter_values[2] ? counter_values[0] < counter_values[2] ? counter_values[0] : counter_values[2]
		 : counter_values[3] < counter_values[4] ? counter_values[3] : counter_values[4]
		 : counter_values[3] < counter_values[4] ? counter_values[1] < counter_values[3] ? counter_values[0] < counter_values[3] ? counter_values[0] : counter_values[3]
		 : counter_values[4] < counter_values[1] ? counter_values[4] : counter_values[1]
		 : counter_values[1] < counter_values[4] ? counter_values[0] < counter_values[4] ? counter_values[0] : counter_values[4]
		 : counter_values[3] < counter_values[1] ? counter_values[3] : counter_values[1]
		 : counter_values[3] < counter_values[2] ? counter_values[0] < counter_values[3] ? counter_values[1] < counter_values[4] ? counter_values[1] < counter_values[3] ? counter_values[4] < counter_values[3] ? counter_values[4] : counter_values[3]
		 : counter_values[2] < counter_values[1] ? counter_values[2] : counter_values[1]
		 : counter_values[4] < counter_values[3] ? counter_values[1] < counter_values[3] ? counter_values[1] : counter_values[3]
		 : counter_values[2] < counter_values[4] ? counter_values[2] : counter_values[4]
		 : counter_values[2] < counter_values[4] ? counter_values[0] < counter_values[2] ? counter_values[1] < counter_values[2] ? counter_values[1] : counter_values[2]
		 : counter_values[4] < counter_values[0] ? counter_values[4] : counter_values[0]
		 : counter_values[0] < counter_values[4] ? counter_values[1] < counter_values[4] ? counter_values[1] : counter_values[4]
		 : counter_values[2] < counter_values[0] ? counter_values[2] : counter_values[0]
		 : counter_values[0] < counter_values[2] ? counter_values[1] < counter_values[4] ? counter_values[1] < counter_values[2] ? counter_values[4] < counter_values[2] ? counter_values[4] : counter_values[2]
		 : counter_values[3] < counter_values[1] ? counter_values[3] : counter_values[1]
		 : counter_values[4] < counter_values[2] ? counter_values[1] < counter_values[2] ? counter_values[1] : counter_values[2]
		 : counter_values[3] < counter_values[4] ? counter_values[3] : counter_values[4]
		 : counter_values[3] < counter_values[4] ? counter_values[0] < counter_values[3] ? counter_values[1] < counter_values[3] ? counter_values[1] : counter_values[3]
		 : counter_values[4] < counter_values[0] ? counter_values[4] : counter_values[0]
		 : counter_values[0] < counter_values[4] ? counter_values[1] < counter_values[4] ? counter_values[1] : counter_values[4]
		 : counter_values[3] < counter_values[0] ? counter_values[3] : counter_values[0];
}

CountSketchBaselineFiveRows::CountSketchBaselineFiveRows(int width, int seed) : CountSketchBaseline(width, 5, seed)
{
}

CountSketchBaselineFiveRows::~CountSketchBaselineFiveRows()
{
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

FineGrainedSalsaCountSketch::FineGrainedSalsaCountSketch(int counterSize, int width, int height, int seed, int downsamplings_per_merge, bool maximum_merge) :
	////
	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX)),
	RandByteArray((char *)&rand_bits)
	////
{

	inputSeed = seed;

	this->counterSize = counterSize;
	this->width = width;
	this->height = height;
	this->maximum_merge = maximum_merge;

	counters_per_word = 64 / counterSize;

	logCounterSize = 0;
	int index = counterSize;
	while (index >>= 1) ++logCounterSize;

	log_counters_per_word = 0;
	index = counters_per_word;
	while (index >>= 1) ++log_counters_per_word;

	width_and_sign_mask = (width << 1) - 1;

	assert((counterSize & (counterSize - 1)) == 0 && "We assume that counterSize is a power of 2!");
	//assert(counterSize > 1 && "We assume!");
	assert(height == 5 && "We assume!");
	assert(width != 0 && "We assume too much!");
	assert(width % 64 == 0 && "We assume that (w % 64 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");
	assert(width < (1 << 30) && "We assume that you are not a database lunatic");

	counters = new uint64_t*[height];
	merges_and_signs = new uint64_t*[height];
	merges_and_signs_32 = new uint32_t*[height];

	bobhash = new BOBHash[height];

	for (int row = 0; row < height; ++row)
	{
		counters[row] = new uint64_t[width / counters_per_word]();
		merges_and_signs[row] = new uint64_t[width >> 5]();
		merges_and_signs_32[row] = (uint32_t*)merges_and_signs[row];

		bobhash[row].initialize(seed*(7 + row) + row + 100);
	}

	////
	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	counter_indices = new int[height];
	flow_signs = new int[height];
	log_counter_sizes = new uint8_t[height];
	counter_values = new int64_t[height];
	counter_bitSizes = new uint8_t[height];

	largest_log_counter = 0;
	this->downsamplings_per_merge = downsamplings_per_merge;
	remaining_downsamplings_per_merge = downsamplings_per_merge;

	this->split_counters_flag = split_counters_flag;
	////
}

FineGrainedSalsaCountSketch::~FineGrainedSalsaCountSketch()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] counters[i];
		delete[] merges_and_signs[i];
	}

	delete[] bobhash;
	delete[] counters;
	delete[] merges_and_signs;

	////
	delete[] counter_indices;
	delete[] counter_values;
	delete[] log_counter_sizes;
	delete[] counter_bitSizes;
	delete[] flow_signs;
	////
}

inline uint8_t FineGrainedSalsaCountSketch::log_counter_size(uint index, uint row)
{

	// 32 counters in each merge word  
	uint32_t merge_word_index = (index >> 5) << 1;

	// index % 32
	uint16_t index_mod_32 = index & 0b11111;

	uint32_t mod_2_bit_mask = 1 << (index_mod_32 & 0b11110);
	uint32_t mod_4_bit_mask = 1 << ((index_mod_32 & 0b11100) | 0b1);
	uint32_t mod_8_bit_mask = 1 << ((index_mod_32 & 0b11000) | 0b11);
	uint32_t mod_16_bit_mask = 1 << ((index_mod_32 & 0b10000) | 0b111);
	uint32_t mod_32_bit_mask = 1 << ((index_mod_32 & 0b00000) | 0b1111);

	uint32_t &merges_word = merges_and_signs_32[row][merge_word_index];

	bool mm2 = ((merges_word & mod_2_bit_mask) != 0);
	bool mm4 = ((merges_word & mod_4_bit_mask) != 0);
	bool mm8 = ((merges_word & mod_8_bit_mask) != 0);
	bool mm16 = ((merges_word & mod_16_bit_mask) != 0);
	bool mm32 = ((merges_word & mod_32_bit_mask) != 0);

	return mm2 + mm4 + mm8 + mm16 + mm32;

}

void FineGrainedSalsaCountSketch::merge_counters(uint index, uint row, uint8_t log_counter_size_value)
{
	uint32_t merge_word_index = (index >> 5) << 1;

	uint32_t merge_bit_mask = (index & 0b11111) & (~((1 << (log_counter_size_value + 1)) - 1));
	merge_bit_mask |= ((1 << (log_counter_size_value)) - 1);
	merge_bit_mask = 1 << merge_bit_mask;
	
	/*
	uint32_t merge_bit_mask;
	switch (log_counter_size_value) {
	case 0:
		merge_bit_mask = 1 << (index & 0b11110);
		break;
	case 1:
		merge_bit_mask = 1 << ((index & 0b11100) | 0b1);
		break;
	case 2:
		merge_bit_mask = 1 << ((index & 0b11000) | 0b11);
		break;
	case 3:
		merge_bit_mask = 1 << ((index & 0b10000) | 0b111);
		break;
	case 4:
		merge_bit_mask = 1 << ((index & 0b00000) | 0b1111);
		break;
	default:
		cout << "Bob!" << endl;
		exit(1);
	}
	*/

	uint32_t &merges_word = merges_and_signs_32[row][merge_word_index];

	// log_counter_size_value == 0: 

	if ((merges_word & merge_bit_mask))
	{
		// counters already merged
	}
	else
	{
		if (log_counter_size_value > 0)
		{
			uint32_t shifted_index = index >> log_counter_size_value;
			uint32_t paired_shifted_index = shifted_index ^ 0b1;
			uint32_t paired_index = paired_shifted_index << log_counter_size_value;
			merge_counters(paired_index, row, log_counter_size_value - 1);
			merge_counters(index, row, log_counter_size_value - 1);
		}
		int64_t counter_pair_sum = compute_counter_pair_sum(index, row, log_counter_size_value);
		merges_word |= merge_bit_mask;
		set_non_overflow_counter(index, row, log_counter_size_value + 1, counter_pair_sum);
	}
}

void FineGrainedSalsaCountSketch::set_non_overflow_counter(uint index, uint row, uint8_t log_counter_size_value, int64_t new_value)
{
	int counter_word_index = index >> log_counters_per_word;
	int offset = ((index >> log_counter_size_value) << log_counter_size_value) & (counters_per_word - 1);
	int bitOffset = offset * counterSize;
	int bitSize = counterSize << log_counter_size_value;
	int negativeBitOffset = 64 - bitOffset - bitSize;
	uint64_t ones = ((uint64_t)1 << (counterSize*(1 << log_counter_size_value))) - 1;
	uint64_t counter_mask = ones << negativeBitOffset;

	counters[row][counter_word_index] &= ~counter_mask;
	counters[row][counter_word_index] |= ((new_value >= 0 ? new_value : -new_value) << negativeBitOffset);

	int zero_one_sign = (new_value >= 0 ? 1 : 0);

	uint32_t sign_word_index = ((index >> 5) << 1) + 1;
	uint32_t &sign_word = merges_and_signs_32[row][sign_word_index];

	//uint32_t first_counter_index = index & (~((1 << log_counter_size_value) - 1));
	uint32_t first_counter_index = (index >> log_counter_size_value) << log_counter_size_value;
	uint32_t first_counter_index_bit_mask = (1 << (0b11111 & first_counter_index));

	if (zero_one_sign == 1) // positive
	{
		sign_word |= first_counter_index_bit_mask;
	}
	else // negative
	{
		sign_word &= ~first_counter_index_bit_mask;
	}
}

void FineGrainedSalsaCountSketch::increment(const char * str)
{
	if (minus_log_p == 0)
	{
		increment_h(str);
	}
	else if (minus_log_p < minus_log_p_bound)
	{
		if (rand_pos < 255)
		{
			++rand_pos;
		}
		else
		{
			rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
			for (int j = 0; j < minus_log_p - 1; ++j)
			{
				rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
			}
			RandByteArray = (char *)&rand_bits;
			rand_pos = 0;
		}

		if (RandByteArray[rand_pos >> 3] & (1 << ((rand_pos & 0x7) - 1))) // TODO: lookup table 
		{
			increment_h(str);
		}
		else
		{
		}
	}
	else
	{
		if (--interval)
		{
		}
		else
		{
			increment_h(str);
			interval = accelerated_geometric(LogOneMinusP);
		}
	}
}

void FineGrainedSalsaCountSketch::increment_h(const char * str)
{

	// remember all counter indices and the sizes
	bool is_overflow = false;
	for (int row = 0; row < height; ++row)
	{

		uint index_and_sign = (bobhash[row].run(str, FT_SIZE)) & width_and_sign_mask;

		flow_signs[row] = 1 - 2 * (index_and_sign & 0b1);
		counter_indices[row] = index_and_sign >> 1;
		counter_values[row] = read_counter(counter_indices[row], row);
		log_counter_sizes[row] = log_counter_size(counter_indices[row], row);
		counter_bitSizes[row] = (1 << ((6 - log_counters_per_word) + log_counter_sizes[row]));

		if ((counter_values[row] > 0) && (flow_signs[row] > 0) && (counter_values[row] == (((uint64_t)1 << counter_bitSizes[row]) - 1)))
		{
			is_overflow = true;
		}
		else if ((counter_values[row] < 0) && (flow_signs[row] < 0) && (-counter_values[row] == (((uint64_t)1 << counter_bitSizes[row]) - 1)))
		{
			is_overflow = true;
		}
	}
	// we are here!
	// safe to increment counters
	if (!is_overflow)
	{
		for (int row = 0; row < height; ++row)
		{
			
			set_non_overflow_counter(counter_indices[row], row, log_counter_sizes[row], counter_values[row] + flow_signs[row]);
		}
	}
	// NOT safe to increment counters
	else
	{
		handle_overflow_and_increment();
	}

}

void FineGrainedSalsaCountSketch::increment_and_potentially_merge()
{
	for (int row = 0; row < height; ++row)
	{

		int index = counter_indices[row];
		int64_t counter_value = counter_values[row];
		int flow_sign = flow_signs[row];
		uint8_t &log_counter_size_value = log_counter_sizes[row];
		int bitSize = counter_bitSizes[row];

		bool is_overflow = (flow_sign*counter_value + 1) == ((uint64_t)1 << bitSize);

		if (!is_overflow)
		{
		}
		else
		{
			merge_counters(index, row, log_counter_size_value);
			++log_counter_size_value;
		}
		set_non_overflow_counter(index, row, log_counter_size_value, read_counter(index, row) + flow_sign);

	}
}

void FineGrainedSalsaCountSketch::handle_overflow_and_increment()
{

	bool is_max_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		if ((flow_signs[row]*counter_values[row] == (((uint64_t)1 << counter_bitSizes[row]) - 1)) && (log_counter_sizes[row] == largest_log_counter))
		{
			is_max_overflow = true;
		}
	}

	if (!is_max_overflow) // merge
	{
	}
	else if (remaining_downsamplings_per_merge--) // do downsampling
	{
		++minus_log_p;
		LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
		divide_counters();
		cout << "divide_counters()" << endl;
		rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
	}
	else
	{
		remaining_downsamplings_per_merge = downsamplings_per_merge;
		++largest_log_counter;
	}
	increment_and_potentially_merge();

}

inline int64_t FineGrainedSalsaCountSketch::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);
	int counter_word_index = index >> log_counters_per_word;

	uint32_t first_counter_index = (index >> log_counter_size_value) << log_counter_size_value;
	uint32_t first_counter_index_bit_mask = (1 << (0b11111 & first_counter_index));

	int offset = first_counter_index & (counters_per_word - 1);
	int bitSize = counterSize * ((uint64_t)1 << log_counter_size_value);
	int bitOffset = offset * counterSize;
	int negativeBitOffset = 64 - bitOffset - bitSize;
	uint64_t ones = ((uint64_t)1 << (counterSize*(1 << log_counter_size_value))) - 1;
	uint64_t unsigned_counter_value = (counters[row][counter_word_index] >> negativeBitOffset) & ones;

	uint32_t sign_word_index = ((index >> 5) << 1) + 1;
	uint32_t &sign_word = merges_and_signs_32[row][sign_word_index];

	// set bit means positive // we are here!
	int counter_sign = 1 - 2 * ((sign_word & first_counter_index_bit_mask) == 0);

	return (int64_t)counter_sign*(int64_t)unsigned_counter_value;
}

inline int64_t FineGrainedSalsaCountSketch::compute_counter_pair_sum(uint index, uint row, uint8_t log_counter_size_value)
{
	int counter_word_index = index >> log_counters_per_word;
	
	uint32_t first_counter_index = (index >> log_counter_size_value) << log_counter_size_value;
	uint32_t first_counter_index_bit_mask = (1 << (0b11111 & first_counter_index));

	int offset = first_counter_index & (counters_per_word - 1);
	int bitSize = counterSize * ((uint64_t)1 << log_counter_size_value);
	int bitOffset = offset * counterSize;
	int negativeBitOffset = 64 - bitOffset - bitSize;

	uint64_t counters_word = counters[row][counter_word_index];

	uint64_t ones = ((uint64_t)1 << (counterSize*(1 << log_counter_size_value))) - 1;
	uint64_t unsigned_counter_value = (counters_word >> negativeBitOffset) & ones;

	uint32_t merge_word_index = (index >> 5) << 1;
	uint32_t sign_word_index = merge_word_index + 1;

	uint32_t shifted_index = index >> log_counter_size_value;
	uint32_t paired_shifted_index = shifted_index ^ 0b1;
	uint32_t paired_index = paired_shifted_index << log_counter_size_value;

	uint32_t paired_first_counter_index_bit_mask = (1 << (paired_index & 0b11111));

	int paired_offset = ((paired_index >> log_counter_size_value) << log_counter_size_value) & (counters_per_word - 1);
	int paired_bitOffset = paired_offset * counterSize;
	int paired_negativeBitOffset = 64 - paired_bitOffset - bitSize;

	uint64_t unsigned_paired_counter_value = (counters_word >> paired_negativeBitOffset) & ones;
	
	uint32_t &sign_word = merges_and_signs_32[row][sign_word_index];

	// set bit means positive 
	int counter_sign = 1 - 2 * ((sign_word & first_counter_index_bit_mask) == 0);
	int paired_counter_sign = 1 - 2 * ((sign_word & paired_first_counter_index_bit_mask) == 0);

	if ((maximum_merge) && (counter_sign == paired_counter_sign))
	{
		return unsigned_counter_value > unsigned_paired_counter_value ? ((int64_t)counter_sign*(int64_t)unsigned_counter_value) : ((int64_t)paired_counter_sign*(int64_t)unsigned_paired_counter_value);
	}

	return ((int64_t)counter_sign*(int64_t)unsigned_counter_value) + ((int64_t)paired_counter_sign*(int64_t)unsigned_paired_counter_value);
}

uint64_t FineGrainedSalsaCountSketch::query(const char * str)
{
	for (int row = 0; row < height; ++row) {
		uint index_and_sign = (bobhash[row].run(str, FT_SIZE)) & width_and_sign_mask;
		uint index = index_and_sign >> 1;
		int sign = 1 - 2 * (index_and_sign & 0b1);

		int64_t singned_counter_value = sign*read_counter(index, row);
		counter_values[row] = singned_counter_value > 0 ? singned_counter_value : 0;
	}
	return MedianCounterValue() << minus_log_p;
}

void FineGrainedSalsaCountSketch::divide_counters()
{
	//cout << "divide_counters" << endl;
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{

			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			int64_t counter_value = read_counter(current_counter_index, row);

			// not symmetric with respect to the sign (shift works differently for positive and negative values)
			set_non_overflow_counter(current_counter_index, row, log_counter_size_value, counter_value >> 1);

			current_counter_index += (1 << log_counter_size_value);
		}
	}
}

long double FineGrainedSalsaCountSketch::l2_estimation()
{
	uint64_t* l2_sos_array = new uint64_t[height]();

	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			int64_t counter_value = read_counter(current_counter_index, row);
			l2_sos_array[row] += counter_value * counter_value;

			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			current_counter_index += (1 << log_counter_size_value);
		}
	}

	sort(l2_sos_array, l2_sos_array + height);
	long double result = sqrtl((long double)l2_sos_array[height >> 1]);
	delete[] l2_sos_array;
	return result;
}

int64_t FineGrainedSalsaCountSketch::signedQuery(const char * str)
{
	int64_t* signedCounterValues = new int64_t[height]();
	for (int row = 0; row < height; ++row) {
		uint index_and_sign = (bobhash[row].run(str, FT_SIZE)) & width_and_sign_mask;
		uint index = index_and_sign >> 1;
		int sign = 1 - 2 * (index_and_sign & 0b1);

		int64_t singned_counter_value = sign * read_counter(index, row);
		signedCounterValues[row] = singned_counter_value;
	}
	sort(signedCounterValues, signedCounterValues + height);
	int64_t estimate = signedCounterValues[height >> 1];
	delete[] signedCounterValues;
	return estimate << minus_log_p;
}

FineGrainedSalsaCountSketch* FineGrainedSalsaCountSketch::operator-(FineGrainedSalsaCountSketch & subtrahend)
{

	FineGrainedSalsaCountSketch* result = new FineGrainedSalsaCountSketch(counterSize, width, height, inputSeed, downsamplings_per_merge, maximum_merge);

	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			
			uint8_t self_log_counter_size_value = log_counter_size(current_counter_index, row);
			int64_t self_counter_value = read_counter(current_counter_index, row);

			uint8_t subtrahend_log_counter_size_value = subtrahend.log_counter_size(current_counter_index, row);
			int64_t subtrahend_counter_value = subtrahend.read_counter(current_counter_index, row);

			uint8_t log_counter_size_value = self_log_counter_size_value > subtrahend_log_counter_size_value ? self_log_counter_size_value : subtrahend_log_counter_size_value;

			int64_t diff_value;

			if (self_log_counter_size_value > subtrahend_log_counter_size_value)
			{
				diff_value = self_counter_value;

				int current_subtrahend_counter_index = current_counter_index;
				int stopping_counter_index = current_counter_index + (1 << self_log_counter_size_value);

				while (current_subtrahend_counter_index < stopping_counter_index)
				{
					diff_value -= subtrahend.read_counter(current_subtrahend_counter_index, row);
					current_subtrahend_counter_index += (1 << subtrahend.log_counter_size(current_subtrahend_counter_index, row));
				}
			}
			else if (self_log_counter_size_value < subtrahend_log_counter_size_value)
			{
				diff_value = -subtrahend_counter_value;

				int current_self_counter_index = current_counter_index;
				int stopping_counter_index = current_counter_index + (1 << subtrahend_log_counter_size_value);

				while (current_self_counter_index < stopping_counter_index)
				{
					diff_value += read_counter(current_self_counter_index, row);
					current_self_counter_index += (1 << log_counter_size(current_self_counter_index, row));
				}
			}
			else
			{
				diff_value = self_counter_value - subtrahend_counter_value;
			}

			///////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////

			bool log_counter_size_value_is_large_enough = abs(diff_value) < ((uint64_t)1 << (counterSize * (1 << log_counter_size_value)));

			if (log_counter_size_value_is_large_enough)
			{
				if (log_counter_size_value > 0)
				{
					result->merge_counters(current_counter_index, row, log_counter_size_value - 1);
				}				
				result->set_non_overflow_counter(current_counter_index, row, log_counter_size_value, diff_value);
			}
			else
			{
				++log_counter_size_value;

				int index_from = (current_counter_index >> log_counter_size_value) << log_counter_size_value;
				int index_to = index_from + (1 << log_counter_size_value);

				int64_t new_diff_value = 0;

				int new_diff_current_self_counter_index = index_from;
				while (new_diff_current_self_counter_index < index_to)
				{

					uint8_t new_diff_self_log_counter_size_value = log_counter_size(new_diff_current_self_counter_index, row);
					int64_t new_diff_self_counter_value = read_counter(new_diff_current_self_counter_index, row);

					new_diff_value += new_diff_self_counter_value;

					new_diff_current_self_counter_index += (1 << new_diff_self_log_counter_size_value);
				}

				int new_diff_current_subtrahend_counter_index = index_from;
				while (new_diff_current_subtrahend_counter_index < index_to)
				{

					uint8_t new_diff_subtrahend_log_counter_size_value = log_counter_size(new_diff_current_subtrahend_counter_index, row);
					int64_t new_diff_subtrahend_counter_value = read_counter(new_diff_current_subtrahend_counter_index, row);

					new_diff_value -= new_diff_subtrahend_counter_value;

					new_diff_current_subtrahend_counter_index += (1 << new_diff_subtrahend_log_counter_size_value);
				}

				result->merge_counters(current_counter_index, row, log_counter_size_value);
				result->set_non_overflow_counter(current_counter_index, row, log_counter_size_value, new_diff_value);
			}

			///////////////////////////////////////////////////////////////////////
			///////////////////////////////////////////////////////////////////////

			/*
			uint8_t res_log_counter_size_value = result->log_counter_size(current_counter_index, row);
			int64_t res_current_value = result->read_counter(current_counter_index, row);

			// overflow during diff
			if (abs(diff_value + res_current_value) >= (1 << (counterSize * (1 << res_log_counter_size_value))))
			{
				//int shifted_current_counter_index = current_counter_index >> log_counter_size_value;
				//int paired_shifted_current_counter_index = shifted_current_counter_index ^ 0b1;
				//int paired_current_counter_index = paired_shifted_current_counter_index << log_counter_size_value;
				uint8_t max_log_counter_size_value = res_log_counter_size_value > log_counter_size_value ? res_log_counter_size_value : log_counter_size_value;

				result->merge_counters(current_counter_index, row, max_log_counter_size_value + 1);
				int64_t current_res_value = result->read_counter(current_counter_index, row);
				result->set_non_overflow_counter(current_counter_index, row, log_counter_size_value + 1, diff_value + current_res_value);
				
			}
			else
			{
				result->set_non_overflow_counter(current_counter_index, row, res_log_counter_size_value, diff_value + res_current_value);
			}			
			*/
			current_counter_index += (1 << log_counter_size_value);
		}
	}

	return result;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

FineGrainedSalsaCountSketchHH::FineGrainedSalsaCountSketchHH(int counterSize, int width, int height, int seed, int downsamplings_per_merge, int HH_num, bool maximum_merge) :
	FineGrainedSalsaCountSketch(counterSize, width, height, seed, downsamplings_per_merge, maximum_merge),
	topK(HH_num)
{
}

FineGrainedSalsaCountSketchHH::~FineGrainedSalsaCountSketchHH()
{
}

void FineGrainedSalsaCountSketchHH::increment(const char * str)
{
	FineGrainedSalsaCountSketch::increment(str);
	string s;
	s.assign(str, FT_SIZE);
	topK.update(s, query(str));
}

vector<pair<string, uint64_t>> FineGrainedSalsaCountSketchHH::HH()
{
	return topK.items();
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
