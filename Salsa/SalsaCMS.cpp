#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "SalsaCMS.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

SalsaCMSSanity::SalsaCMSSanity(int width, int height, int seed, int downsamplings_per_merge):
	////
	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX)),
	RandByteArray((char *)&rand_bits)
	////
{
	this->width = width;
	this->height = height;

	//width_mask = width - 1;
	width_mask = (width - 1) << 3;

	assert(width > 0 && "We assume too much!");
	assert(width % 16 == 0 && "We assume that (w % 16 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cms_counters = new uint8_t*[height];

	baseline_cms_counters_16 = new uint16_t*[height];
	baseline_cms_counters_32 = new uint32_t*[height];
	baseline_cms_counters_64 = new uint64_t*[height];

	baseline_cms_merges = new uint16_t*[height];
	bobhash = new BOBHash[height];

	for (int row = 0; row < height; ++row)
	{
		baseline_cms_counters[row] = new uint8_t[width]();

		baseline_cms_counters_16[row] = (uint16_t*)baseline_cms_counters[row];
		baseline_cms_counters_32[row] = (uint32_t*)baseline_cms_counters[row];
		baseline_cms_counters_64[row] = (uint64_t*)baseline_cms_counters[row];

		baseline_cms_merges[row] = new uint16_t[width >> 4]();
		bobhash[row].initialize(seed*(7 + row) + row + 100);
	}

	for (int index_mod_16 = 0; index_mod_16 < 16; ++index_mod_16)
	{
		uint16_t mod_2_bit_mask = 1 << (index_mod_16 & 0b1110);
		uint16_t mod_4_bit_mask = 1 << ((index_mod_16 & 0b1100) | 0b1);
		uint16_t mod_8_bit_mask = 1 << ((index_mod_16 & 0b1000) | 0b11);

		merges_lookup_table_8[index_mod_16] = mod_2_bit_mask;
		merges_lookup_table_16[index_mod_16] = mod_4_bit_mask;
		merges_lookup_table_32[index_mod_16] = mod_8_bit_mask;
	}

	////
	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	counter_indices = new int[height];
	log_counter_sizes = new uint8_t[height];
	counter_values = new uint64_t[height];

	largest_log_counter = 0;
	this->downsamplings_per_merge = downsamplings_per_merge;
	remaining_downsamplings_per_merge = downsamplings_per_merge;
	////
}

SalsaCMSSanity::~SalsaCMSSanity()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms_counters[i];
		delete[] baseline_cms_merges[i];
	}

	delete[] bobhash;
	delete[] baseline_cms_counters;

	delete[] baseline_cms_counters_16;
	delete[] baseline_cms_counters_32;
	delete[] baseline_cms_counters_64;

	delete[] baseline_cms_merges;

	////
	delete[] counter_indices;
	delete[] counter_values;
	delete[] log_counter_sizes;
	////
}

inline uint8_t SalsaCMSSanity::log_counter_size(uint index, uint row)
{

	// 16 counters in each index  
	uint32_t byte_index = index >> 4;

	// index % 16
	uint16_t index_mod_16 = index & 0b1111;

	uint16_t mod_2_bit_mask = merges_lookup_table_8[index_mod_16]; // 1 << (index_mod_16 & 0b1110);
	uint16_t mod_4_bit_mask = merges_lookup_table_16[index_mod_16]; // 1 << ((index_mod_16 & 0b1100) | 0b1);
	uint16_t mod_8_bit_mask = merges_lookup_table_32[index_mod_16]; // 1 << ((index_mod_16 & 0b1000) | 0b11);

	uint16_t merges_word = baseline_cms_merges[row][byte_index];

	//uint16_t mm2 = merges_word & mod_2_bit_mask;
	//uint16_t mm4 = merges_word & mod_4_bit_mask;
	//uint16_t mm8 = merges_word & mod_8_bit_mask;

	return (merges_word & mod_2_bit_mask) != 0 ? ((merges_word & mod_4_bit_mask) == 0 ? 1 : (merges_word & mod_8_bit_mask) == 0 ? 2 : 3) : 0;

	// we assume that counter sizes cannot exceed 8X the original
	// if the original is 8 bits then we have a 64-bit limit

}

void SalsaCMSSanity::merge_8_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_2_bit_mask = 1 << (index & 0b1110);

	if ((baseline_cms_merges[row][byte_index] & mod_2_bit_mask))
	{
		// counters already merged
	}
	else
	{
		baseline_cms_merges[row][byte_index] |= mod_2_bit_mask;

		// happens in increment
		//uint16_t* p = (uint16_t*)&baseline_cms_counters[row][index & 0xFFFFFFFE];
		//*p = max(baseline_cms_counters[row][index], baseline_cms_counters[row][index ^ 0b1]);
	}
}

void SalsaCMSSanity::merge_16_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_4_bit_mask = 1 << ((index & 0b1100) | 0b1);

	if ((baseline_cms_merges[row][byte_index] & mod_4_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint index_16 = index >> 1;
		uint paired_index_16 = index_16 ^ 0b1;

		// merge other counters if needed
		merge_8_bit_counters(index, row);
		merge_8_bit_counters(paired_index_16 << 1, row);

		baseline_cms_merges[row][byte_index] |= mod_4_bit_mask;

		// happens in increment
		//uint32_t* p = (uint32_t*)&baseline_cms_counters_16[row][index_16 & 0xFFFFFFFE];		
		//*p = max(baseline_cms_counters_16[row][index_16], baseline_cms_counters_16[row][index_16 ^ 0b1]);
	}
}

void SalsaCMSSanity::merge_32_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_8_bit_mask = 1 << ((index & 0b1000) | 0b11);

	if ((baseline_cms_merges[row][byte_index] & mod_8_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint index_32 = index >> 2;
		uint paired_index_32 = index_32 ^ 0b1;

		// merge other counters if needed
		merge_16_bit_counters(index, row);
		merge_16_bit_counters(paired_index_32 << 2, row);

		baseline_cms_merges[row][byte_index] |= mod_8_bit_mask;

		// happens in increment
		//uint64_t* p = (uint64_t*)&baseline_cms_counters_16[row][index_32 & 0xFFFFFFFE];		
		//*p = max(baseline_cms_counters_32[row][index_32], baseline_cms_counters_32[row][index_32 ^ 0b1]);
	}
}

inline uint64_t SalsaCMSSanity::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);

	switch (log_counter_size_value) {
	case 0:
		return baseline_cms_counters[row][index];
	case 1:
		return baseline_cms_counters_16[row][index >> 1];
	case 2:
		return baseline_cms_counters_32[row][index >> 2];
	case 3:
		return baseline_cms_counters_64[row][index >> 3];
	case 4:
		// not implemented
		return -1;
	default:
		// BOB!
		return -1;
	}
}

uint64_t SalsaCMSSanity::query(const char * str)
{
	int row = 0;

	//uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
	uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> 3;
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		//uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> 3;
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return ((uint64_t)min) << minus_log_p;
}

void SalsaCMSSanity::increment_h(const char * str)
{

	// remember all counter indices and the sizes
	bool is_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		//counter_indices[row] = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		counter_indices[row] = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> 3;
		counter_values[row] = read_counter(counter_indices[row], row);
		log_counter_sizes[row] = log_counter_size(counter_indices[row], row);

		if (counter_values[row] == (((uint64_t)1 << (1 << (3 + log_counter_sizes[row]))) - 1))
		{
			is_overflow = true;
		}
	}

	// safe to increment counters
	if (!is_overflow)
	{
		for (int row = 0; row < height; ++row)
		{
			switch (log_counter_sizes[row])
			{
			case 0:
				++baseline_cms_counters[row][counter_indices[row]];
				break;
			case 1:
				++baseline_cms_counters_16[row][counter_indices[row] >> 1];
				break;
			case 2:
				++baseline_cms_counters_32[row][counter_indices[row] >> 2];
				break;
			default:
				++baseline_cms_counters_64[row][counter_indices[row] >> 3];
			}
		}
	}
	// NOT safe to increment counters
	else
	{
		handle_overflow_and_increment();
	}

}

void SalsaCMSSanity::increment(const char * str)
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

void SalsaCMSSanity::divide_counters()
{
	//cout << "divide_counters" << endl;
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			switch (log_counter_size_value) 
			{
			case 0:
				baseline_cms_counters[row][current_counter_index] >>= 1;
				break;
			case 1:
				baseline_cms_counters_16[row][current_counter_index >> 1] >>= 1;
				break;
			case 2:
				baseline_cms_counters_32[row][current_counter_index >> 2] >>= 1;
				break;
			default:
				baseline_cms_counters_64[row][current_counter_index >> 3] >>= 1;
			}
			current_counter_index += (1 << log_counter_size_value);
		}
	}

}

void SalsaCMSSanity::increment_and_potentially_merge()
{
	for (int row = 0; row < height; ++row) {

		int index = counter_indices[row];
		uint64_t value = counter_values[row];
		uint8_t log_counter_size_value = log_counter_sizes[row];

		// we are here... 
		if (log_counter_size_value == 0)
		{
			if (++baseline_cms_counters[row][index])
			{
				// no overflow
			}
			else
			{
				merge_8_bit_counters(index, row);

				// the current counter became 0 while its real value is 256...
				uint index_16 = index >> 1;
				baseline_cms_counters_16[row][index_16] = 1 << 8;
			}
		}
		else if (log_counter_size_value == 1)
		{
			uint index_16 = index >> 1;
			if (++baseline_cms_counters_16[row][index_16])
			{
				// no overflow
			}
			else
			{
				merge_16_bit_counters(index, row);

				// the current counter became 0 while its real value is 2^16...
				uint index_32 = index >> 2;
				baseline_cms_counters_32[row][index_32] = 1 << 16;
			}
		}
		else if (log_counter_size_value == 2)
		{
			uint index_32 = index >> 2;
			if (++baseline_cms_counters_32[row][index_32])
			{
				// no overflow
			}
			else
			{
				merge_32_bit_counters(index, row);

				// the current counter became 0 while its real value is 2^32...
				uint index_64 = index >> 3;
				baseline_cms_counters_64[row][index_64] = ((uint64_t)1) << 32;
			}
		}
		else // if (log_counter_size_value == 3)
		{
			uint index_64 = index >> 3;
			++baseline_cms_counters_64[row][index_64];
		}
	}
}

void SalsaCMSSanity::handle_overflow_and_increment()
{

	bool is_max_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		if ((counter_values[row] == (((uint64_t)1 << (1 << (3 + log_counter_sizes[row]))) - 1)) && (log_counter_sizes[row] == largest_log_counter))
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

		for (int row = 0; row < height; ++row)
		{
			counter_values[row] >>= 1;
		}

		rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
	}
	else
	{
		remaining_downsamplings_per_merge = downsamplings_per_merge;
		++largest_log_counter;
	}
	increment_and_potentially_merge();

}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

FineGrainedSalsaCMSSanity::FineGrainedSalsaCMSSanity(int counterSize, int width, int height, int seed, int downsamplings_per_merge, bool split_counters_flag) :
	////
	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX)),
	RandByteArray((char *)&rand_bits)
	////
{

	this->counterSize = counterSize;
	this->width = width;
	this->height = height;

	counters_per_word = 64 / counterSize;

	logCounterSize = 0;
	int index = counterSize;
	while (index >>= 1) ++logCounterSize;

	log_counters_per_word = 0;
	index = counters_per_word;
	while (index >>= 1) ++log_counters_per_word;

	width_mask = (width - 1) << (6 - log_counters_per_word);

	assert((counterSize & (counterSize - 1)) == 0 && "We assume that counterSize is a power of 2!");

	assert(width > 0 && "We assume too much!");
	assert(width % 64 == 0 && "We assume that (w % 64 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cms_counters = new uint64_t*[height];
	baseline_cms_merges = new uint64_t*[height];

	bobhash = new BOBHash[height];

	for (int row = 0; row < height; ++row)
	{
		baseline_cms_counters[row] = new uint64_t[width / counters_per_word]();
		baseline_cms_merges[row] = new uint64_t[width >> 6]();

		bobhash[row].initialize(seed*(7 + row) + row + 100);
	}

	/*
	for (int index_mod_16 = 0; index_mod_16 < 16; ++index_mod_16)
	{
		uint16_t mod_2_bit_mask = 1 << (index_mod_16 & 0b1110);
		uint16_t mod_4_bit_mask = 1 << ((index_mod_16 & 0b1100) | 0b1);
		uint16_t mod_8_bit_mask = 1 << ((index_mod_16 & 0b1000) | 0b11);

		merges_lookup_table_8[index_mod_16] = mod_2_bit_mask;
		merges_lookup_table_16[index_mod_16] = mod_4_bit_mask;
		merges_lookup_table_32[index_mod_16] = mod_8_bit_mask;
	}
	*/

	////
	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	counter_indices = new int[height];
	log_counter_sizes = new uint8_t[height];
	counter_values = new uint64_t[height];
	counter_bitSizes = new uint8_t[height];

	largest_log_counter = 0;
	this->downsamplings_per_merge = downsamplings_per_merge;
	remaining_downsamplings_per_merge = downsamplings_per_merge;

	this->split_counters_flag = split_counters_flag;
	////
}

FineGrainedSalsaCMSSanity::~FineGrainedSalsaCMSSanity()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms_counters[i];
		delete[] baseline_cms_merges[i];
	}

	delete[] bobhash;
	delete[] baseline_cms_counters;
	delete[] baseline_cms_merges;

	////
	delete[] counter_indices;
	delete[] counter_values;
	delete[] log_counter_sizes;
	delete[] counter_bitSizes;
	////
}

inline uint8_t FineGrainedSalsaCMSSanity::log_counter_size(uint index, uint row)
{

	// 64 counters in each merge word  
	uint32_t word_index = index >> 6;

	// index % 64
	uint16_t index_mod_64 = index & 0b111111;

	//uint16_t mod_2_bit_mask = merges_lookup_table_8[index_mod_16];
	//uint16_t mod_4_bit_mask = merges_lookup_table_16[index_mod_16];
	//uint16_t mod_8_bit_mask = merges_lookup_table_32[index_mod_16];

	uint64_t mod_2_bit_mask = (uint64_t)1 << (index_mod_64 & 0b111110);
	uint64_t mod_4_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b111100) | 0b1);
	uint64_t mod_8_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b111000) | 0b11);
	uint64_t mod_16_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b110000) | 0b111);
	uint64_t mod_32_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b100000) | 0b1111);
	uint64_t mod_64_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b000000) | 0b11111);

	uint64_t merges_word = baseline_cms_merges[row][word_index];

	bool mm2 = ((merges_word & mod_2_bit_mask) != 0);
	bool mm4 = ((merges_word & mod_4_bit_mask) != 0);
	bool mm8 = ((merges_word & mod_8_bit_mask) != 0);
	bool mm16 = ((merges_word & mod_16_bit_mask) != 0);
	bool mm32 = ((merges_word & mod_32_bit_mask) != 0);
	bool mm64 = ((merges_word & mod_64_bit_mask) != 0);

	return mm2 + mm4 + mm8 + mm16 + mm32 + mm64;

}

void FineGrainedSalsaCMSSanity::merge_size_1_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_2_bit_mask = ((uint64_t)1 << (index & 0b111110));

	if ((baseline_cms_merges[row][word_index] & mod_2_bit_mask))
	{
		// counters already merged
	}
	else
	{
		baseline_cms_merges[row][word_index] |= mod_2_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSSanity::merge_size_2_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_4_bit_mask = ((uint64_t)1 << ((index & 0b111100) | 0b1));

	if ((baseline_cms_merges[row][word_index] & mod_4_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_1_index = index >> 1;
		uint paired_shifted_1_index = shifted_1_index ^ 0b1;

		// merge other counters if needed
		merge_size_1_counters(index, row);
		merge_size_1_counters(paired_shifted_1_index << 1, row);

		baseline_cms_merges[row][word_index] |= mod_4_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSSanity::merge_size_4_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_8_bit_mask = ((uint64_t)1 << ((index & 0b111000) | 0b11));

	if ((baseline_cms_merges[row][word_index] & mod_8_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_2_index = index >> 2;
		uint paired_shifted_2_index = shifted_2_index ^ 0b1;

		// merge other counters if needed
		merge_size_2_counters(index, row);
		merge_size_2_counters(paired_shifted_2_index << 2, row);

		baseline_cms_merges[row][word_index] |= mod_8_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSSanity::merge_size_8_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_16_bit_mask = ((uint64_t)1 << ((index & 0b110000) | 0b111));

	if ((baseline_cms_merges[row][word_index] & mod_16_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_3_index = index >> 3;
		uint paired_shifted_3_index = shifted_3_index ^ 0b1;

		// merge other counters if needed
		merge_size_4_counters(index, row);
		merge_size_4_counters(paired_shifted_3_index << 3, row);

		baseline_cms_merges[row][word_index] |= mod_16_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSSanity::merge_size_16_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_32_bit_mask = ((uint64_t)1 << ((index & 0b100000) | 0b1111));

	if ((baseline_cms_merges[row][word_index] & mod_32_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_4_index = index >> 4;
		uint paired_shifted_4_index = shifted_4_index ^ 0b1;

		// merge other counters if needed
		merge_size_8_counters(index, row);
		merge_size_8_counters(paired_shifted_4_index << 4, row);

		baseline_cms_merges[row][word_index] |= mod_32_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSSanity::merge_size_32_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_64_bit_mask = ((uint64_t)1 << ((index & 0b000000) | 0b11111));

	if ((baseline_cms_merges[row][word_index] & mod_64_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_5_index = index >> 5;
		uint paired_shifted_5_index = shifted_5_index ^ 0b1;

		// merge other counters if needed
		merge_size_16_counters(index, row);
		merge_size_16_counters(paired_shifted_5_index << 5, row);

		baseline_cms_merges[row][word_index] |= mod_64_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSSanity::set_non_overflow_counter(uint index, uint row, uint8_t log_counter_size_value, uint64_t new_value, int bitSize)
{
	int counter_word_index = index >> log_counters_per_word;
	int offset = ((index >> log_counter_size_value) << log_counter_size_value) & (counters_per_word - 1);
	int bitOffset = offset * counterSize;
	int negativeBitOffset = 64 - bitOffset - bitSize;
	uint64_t ones = ((uint64_t)1 << (counterSize*(1 << log_counter_size_value))) - 1;
	uint64_t counter_mask = ones << negativeBitOffset;

	baseline_cms_counters[row][counter_word_index] &= ~counter_mask;
	baseline_cms_counters[row][counter_word_index] |= (new_value << negativeBitOffset);
}

void FineGrainedSalsaCMSSanity::increment(const char * str)
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

void FineGrainedSalsaCMSSanity::increment_h(const char * str)
{

	// remember all counter indices and the sizes
	bool is_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		counter_indices[row] = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
		counter_values[row] = read_counter(counter_indices[row], row);
		log_counter_sizes[row] = log_counter_size(counter_indices[row], row);
		counter_bitSizes[row] = (1 << ((6 - log_counters_per_word) + log_counter_sizes[row]));

		if (counter_values[row] == (((uint64_t)1 << counter_bitSizes[row]) - 1))
		{
			is_overflow = true;
		}
	}

	// safe to increment counters
	if (!is_overflow)
	{
		for (int row = 0; row < height; ++row)
		{
			set_non_overflow_counter(counter_indices[row], row, log_counter_sizes[row], counter_values[row] + 1, counter_bitSizes[row]);
		}		
	}
	// NOT safe to increment counters
	else
	{
		handle_overflow_and_increment();
	}

}

void FineGrainedSalsaCMSSanity::increment_and_potentially_merge()
{
	for (int row = 0; row < height; ++row)
	{

		int index = counter_indices[row];
		uint64_t counter_value = counter_values[row];
		uint8_t log_counter_size_value = log_counter_sizes[row];
		int bitSize = counter_bitSizes[row];

		bool is_overflow = (counter_value + 1) == ((uint64_t)1 << bitSize);

		if (!is_overflow)
		{
		}
		else
		{
			if (log_counter_size_value == 0)
			{
				merge_size_1_counters(index, row);
			}
			else if (log_counter_size_value == 1)
			{
				merge_size_2_counters(index, row);
			}
			else if (log_counter_size_value == 2)
			{
				merge_size_4_counters(index, row);
			}
			else if (log_counter_size_value == 3)
			{
				merge_size_8_counters(index, row);
			}
			else if (log_counter_size_value == 4)
			{
				merge_size_16_counters(index, row);
			}
			else if (log_counter_size_value == 5)
			{
				merge_size_32_counters(index, row);
			}
			++log_counter_size_value;
			bitSize <<= 1;
		}
		set_non_overflow_counter(index, row, log_counter_size_value, counter_value + 1, bitSize);

	}
}

void FineGrainedSalsaCMSSanity::handle_overflow_and_increment()
{

	bool is_max_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		if ((counter_values[row] == (((uint64_t)1 << counter_bitSizes[row]) - 1)) && (log_counter_sizes[row] == largest_log_counter))
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

		for (int row = 0; row < height; ++row)
		{
			counter_values[row] >>= 1;
		}
		
		if (split_counters_flag)
		{
			split_counters();
		}

		rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
	}
	else
	{
		remaining_downsamplings_per_merge = downsamplings_per_merge;
		++largest_log_counter;
	}
	increment_and_potentially_merge();

}

inline uint64_t FineGrainedSalsaCMSSanity::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);
	int counter_word_index = index >> log_counters_per_word;
	int offset = ((index >> log_counter_size_value) << log_counter_size_value) & (counters_per_word - 1);
	int bitSize = counterSize * ((uint64_t)1 << log_counter_size_value);
	int bitOffset = offset * counterSize;
	int negativeBitOffset = 64 - bitOffset - bitSize;
	uint64_t ones = ((uint64_t)1 << (counterSize*(1 << log_counter_size_value))) - 1;

	return (baseline_cms_counters[row][counter_word_index] >> negativeBitOffset) & ones;
}

uint64_t FineGrainedSalsaCMSSanity::query(const char * str)
{
	int row = 0;

	uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return ((uint64_t)min) << minus_log_p;
}

void FineGrainedSalsaCMSSanity::divide_counters()
{
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			
			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			uint64_t counter_value = read_counter(current_counter_index, row);
			int bitSize = counterSize * (1 << log_counter_size_value);
			
			set_non_overflow_counter(current_counter_index, row, log_counter_size_value, counter_value >> 1, bitSize);

			current_counter_index += (1 << log_counter_size_value);
		}
	}
}

void FineGrainedSalsaCMSSanity::split_counters()
{
	//cout << "split_counters" << endl;
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{

			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			uint64_t counter_value = read_counter(current_counter_index, row);
			int bitSize = counterSize * (1 << log_counter_size_value);

			uint64_t threshold_value = (1 << (bitSize >> 1));
			if ((counter_value < threshold_value) && (log_counter_size_value > 0))
			{
				split_counter(row, current_counter_index, log_counter_size_value, counter_value);
			}
			current_counter_index += (1 << log_counter_size_value);
		}
	}

}

void FineGrainedSalsaCMSSanity::split_counter(int row, int index, uint8_t log_counter_size_value, uint64_t counter_value)
{
	uint64_t bit_mask;

	switch (log_counter_size_value)
	{
	case 0:
		// minimal size counter - cannot split
		return;
	case 1:
		bit_mask = ((uint64_t)1 << (index & 0b111110));
		break;
	case 2:
		bit_mask = ((uint64_t)1 << ((index & 0b111100) | 0b1));
		break;
	case 3:
		bit_mask = ((uint64_t)1 << ((index & 0b111000) | 0b11));
		break;
	case 4:
		bit_mask = ((uint64_t)1 << ((index & 0b110000) | 0b111));
		break;
	case 5:
		bit_mask = ((uint64_t)1 << ((index & 0b100000) | 0b1111));
		break;
	case 6:
		bit_mask = ((uint64_t)1 << ((index & 0b000000) | 0b11111));
		break;
	default:
		cout << "Error in FineGrainedSalsaCMSSanity::split_counter" << endl;
		exit(1);
	}

	int word_index = (index >> 6);
	baseline_cms_merges[row][word_index] &= ~bit_mask;

	uint shifted_index = index >> (log_counter_size_value - 1);
	uint paired_shifted_index = shifted_index ^ 0b1;

	int bitSize = (1 << ((6 - log_counters_per_word) + (log_counter_size_value - 1)));

	set_non_overflow_counter(index, row, log_counter_size_value - 1, counter_value, bitSize);
	set_non_overflow_counter(paired_shifted_index << (log_counter_size_value - 1), row, log_counter_size_value - 1, counter_value, bitSize);
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

SalsaAnalyticalErrorCMS::SalsaAnalyticalErrorCMS(int width, int height, int seed, double delta, int downsamplings_before_merge) :
	////
	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX)),
	RandByteArray((char *)&rand_bits)
	////
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 16 == 0 && "We assume that (w % 16 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cms_counters = new uint8_t*[height];

	baseline_cms_counters_16 = new uint16_t*[height];
	baseline_cms_counters_32 = new uint32_t*[height];
	baseline_cms_counters_64 = new uint64_t*[height];

	baseline_cms_merges = new uint16_t*[height];
	bobhash = new BOBHash[height];

	for (int row = 0; row < height; ++row)
	{
		baseline_cms_counters[row] = new uint8_t[width]();

		baseline_cms_counters_16[row] = (uint16_t*)baseline_cms_counters[row];
		baseline_cms_counters_32[row] = (uint32_t*)baseline_cms_counters[row];
		baseline_cms_counters_64[row] = (uint64_t*)baseline_cms_counters[row];

		baseline_cms_merges[row] = new uint16_t[width >> 4]();
		bobhash[row].initialize(seed*(7 + row) + row + 100);
	}

	for (int index_mod_16 = 0; index_mod_16 < 16; ++index_mod_16)
	{
		uint16_t mod_2_bit_mask = 1 << (index_mod_16 & 0b1110);
		uint16_t mod_4_bit_mask = 1 << ((index_mod_16 & 0b1100) | 0b1);
		uint16_t mod_8_bit_mask = 1 << ((index_mod_16 & 0b1000) | 0b11);

		merges_lookup_table_8[index_mod_16] = mod_2_bit_mask;
		merges_lookup_table_16[index_mod_16] = mod_4_bit_mask;
		merges_lookup_table_32[index_mod_16] = mod_8_bit_mask;
	}

	////
	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	counter_indices = new int[height];
	log_counter_sizes = new uint8_t[height];
	counter_values = new uint64_t[height];

	largest_log_counter = 0;
	////
	this->delta = delta;
	epsilonSketch = (pow(1.0 / delta, 1.0 / height)) / width;
	numPackets = 0;
	this->downsamplings_before_merge = downsamplings_before_merge;

}

SalsaAnalyticalErrorCMS::~SalsaAnalyticalErrorCMS()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms_counters[i];
		delete[] baseline_cms_merges[i];
	}

	delete[] bobhash;
	delete[] baseline_cms_counters;

	delete[] baseline_cms_counters_16;
	delete[] baseline_cms_counters_32;
	delete[] baseline_cms_counters_64;

	delete[] baseline_cms_merges;

	////
	delete[] counter_indices;
	delete[] counter_values;
	delete[] log_counter_sizes;
	////
}

inline uint8_t SalsaAnalyticalErrorCMS::log_counter_size(uint index, uint row)
{

	// 16 counters in each index  
	uint32_t byte_index = index >> 4;

	// index % 16
	uint16_t index_mod_16 = index & 0b1111;

	uint16_t mod_2_bit_mask = merges_lookup_table_8[index_mod_16]; // 1 << (index_mod_16 & 0b1110);
	uint16_t mod_4_bit_mask = merges_lookup_table_16[index_mod_16]; // 1 << ((index_mod_16 & 0b1100) | 0b1);
	uint16_t mod_8_bit_mask = merges_lookup_table_32[index_mod_16]; // 1 << ((index_mod_16 & 0b1000) | 0b11);

	uint16_t merges_word = baseline_cms_merges[row][byte_index];

	//uint16_t mm2 = merges_word & mod_2_bit_mask;
	//uint16_t mm4 = merges_word & mod_4_bit_mask;
	//uint16_t mm8 = merges_word & mod_8_bit_mask;

	return (merges_word & mod_2_bit_mask) != 0 ? ((merges_word & mod_4_bit_mask) == 0 ? 1 : (merges_word & mod_8_bit_mask) == 0 ? 2 : 3) : 0;

	// we assume that counter sizes cannot exceed 8X the original
	// if the original is 8 bits then we have a 64-bit limit

}

void SalsaAnalyticalErrorCMS::merge_8_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_2_bit_mask = 1 << (index & 0b1110);

	if ((baseline_cms_merges[row][byte_index] & mod_2_bit_mask))
	{
		// counters already merged
	}
	else
	{
		baseline_cms_merges[row][byte_index] |= mod_2_bit_mask;

		// happens in increment
		//uint16_t* p = (uint16_t*)&baseline_cms_counters[row][index & 0xFFFFFFFE];
		//*p = max(baseline_cms_counters[row][index], baseline_cms_counters[row][index ^ 0b1]);
	}
}

void SalsaAnalyticalErrorCMS::merge_16_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_4_bit_mask = 1 << ((index & 0b1100) | 0b1);

	if ((baseline_cms_merges[row][byte_index] & mod_4_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint index_16 = index >> 1;
		uint paired_index_16 = index_16 ^ 0b1;

		// merge other counters if needed
		merge_8_bit_counters(index, row);
		merge_8_bit_counters(paired_index_16 << 1, row);

		baseline_cms_merges[row][byte_index] |= mod_4_bit_mask;

		// happens in increment
		//uint32_t* p = (uint32_t*)&baseline_cms_counters_16[row][index_16 & 0xFFFFFFFE];		
		//*p = max(baseline_cms_counters_16[row][index_16], baseline_cms_counters_16[row][index_16 ^ 0b1]);
	}
}

void SalsaAnalyticalErrorCMS::merge_32_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_8_bit_mask = 1 << ((index & 0b1000) | 0b11);

	if ((baseline_cms_merges[row][byte_index] & mod_8_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint index_32 = index >> 2;
		uint paired_index_32 = index_32 ^ 0b1;

		// merge other counters if needed
		merge_16_bit_counters(index, row);
		merge_16_bit_counters(paired_index_32 << 2, row);

		baseline_cms_merges[row][byte_index] |= mod_8_bit_mask;

		// happens in increment
		//uint64_t* p = (uint64_t*)&baseline_cms_counters_16[row][index_32 & 0xFFFFFFFE];		
		//*p = max(baseline_cms_counters_32[row][index_32], baseline_cms_counters_32[row][index_32 ^ 0b1]);
	}
}

inline uint64_t SalsaAnalyticalErrorCMS::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);

	switch (log_counter_size_value) {
	case 0:
		return baseline_cms_counters[row][index];
	case 1:
		return baseline_cms_counters_16[row][index >> 1];
	case 2:
		return baseline_cms_counters_32[row][index >> 2];
	case 3:
		return baseline_cms_counters_64[row][index >> 3];
	case 4:
		// not implemented
		return -1;
	default:
		// BOB!
		return -1;
	}
}

uint64_t SalsaAnalyticalErrorCMS::query(const char * str)
{
	int row = 0;

	//uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
	uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		//uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return ((uint64_t)min) << minus_log_p;
}

void SalsaAnalyticalErrorCMS::increment_h(const char * str)
{

	// remember all counter indices and the sizes
	bool is_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		counter_indices[row] = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		counter_values[row] = read_counter(counter_indices[row], row);
		log_counter_sizes[row] = log_counter_size(counter_indices[row], row);

		if (counter_values[row] == (((uint64_t)1 << (1 << (3 + log_counter_sizes[row]))) - 1))
		{
			is_overflow = true;
		}
	}

	// safe to increment counters
	if (!is_overflow)
	{
		for (int row = 0; row < height; ++row)
		{
			switch (log_counter_sizes[row])
			{
			case 0:
				++baseline_cms_counters[row][counter_indices[row]];
				break;
			case 1:
				++baseline_cms_counters_16[row][counter_indices[row] >> 1];
				break;
			case 2:
				++baseline_cms_counters_32[row][counter_indices[row] >> 2];
				break;
			default:
				++baseline_cms_counters_64[row][counter_indices[row] >> 3];
			}
		}
	}
	// NOT safe to increment counters
	else
	{
		handle_overflow_and_increment();
	}

}

void SalsaAnalyticalErrorCMS::increment(const char * str)
{
	++numPackets;
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

void SalsaAnalyticalErrorCMS::divide_counters()
{
	//cout << "divide_counters" << endl;
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			switch (log_counter_size_value)
			{
			case 0:
				baseline_cms_counters[row][current_counter_index] >>= 1;
				break;
			case 1:
				baseline_cms_counters_16[row][current_counter_index >> 1] >>= 1;
				break;
			case 2:
				baseline_cms_counters_32[row][current_counter_index >> 2] >>= 1;
				break;
			default:
				baseline_cms_counters_64[row][current_counter_index >> 3] >>= 1;
			}
			current_counter_index += (1 << log_counter_size_value);
		}
	}

}

void SalsaAnalyticalErrorCMS::increment_and_potentially_merge()
{
	for (int row = 0; row < height; ++row) {

		int index = counter_indices[row];
		uint64_t value = counter_values[row];
		uint8_t log_counter_size_value = log_counter_sizes[row];

		// we are here... 
		if (log_counter_size_value == 0)
		{
			if (++baseline_cms_counters[row][index])
			{
				// no overflow
			}
			else
			{
				merge_8_bit_counters(index, row);

				// the current counter became 0 while its real value is 256...
				uint index_16 = index >> 1;
				baseline_cms_counters_16[row][index_16] = 1 << 8;
			}
		}
		else if (log_counter_size_value == 1)
		{
			uint index_16 = index >> 1;
			if (++baseline_cms_counters_16[row][index_16])
			{
				// no overflow
			}
			else
			{
				merge_16_bit_counters(index, row);

				// the current counter became 0 while its real value is 2^16...
				uint index_32 = index >> 2;
				baseline_cms_counters_32[row][index_32] = 1 << 16;
			}
		}
		else if (log_counter_size_value == 2)
		{
			uint index_32 = index >> 2;
			if (++baseline_cms_counters_32[row][index_32])
			{
				// no overflow
			}
			else
			{
				merge_32_bit_counters(index, row);

				// the current counter became 0 while its real value is 2^32...
				uint index_64 = index >> 3;
				baseline_cms_counters_64[row][index_64] = ((uint64_t)1) << 32;
			}
		}
		else // if (log_counter_size_value == 3)
		{
			uint index_64 = index >> 3;
			++baseline_cms_counters_64[row][index_64];
		}
	}
}

void SalsaAnalyticalErrorCMS::handle_overflow_and_increment()
{

	bool is_max_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		if ((counter_values[row] == (((uint64_t)1 << (1 << (3 + log_counter_sizes[row]))) - 1)) && (log_counter_sizes[row] == largest_log_counter))
		{
			is_max_overflow = true;
		}
	}

	if (!is_max_overflow) // merge
	{
	}
	else
	{
		double diffEpsilonSketch = epsilonSketch;

		double newNumPacketsPrime = (double)numPackets / (double)(1 << (minus_log_p + 1));
		double newEpsilonEstimator = sqrt((2.0*log(2.0 / delta)) / newNumPacketsPrime);

		double currentNumPacketsPrime = (double)numPackets / (double)(1 << minus_log_p);
		double currentEpsilonEstimator = sqrt((2.0*log(2.0 / delta)) / currentNumPacketsPrime);

		double diffEpsilonEstimator = newEpsilonEstimator - currentEpsilonEstimator;

		//cout << "diffEpsilonSketch " << diffEpsilonSketch << endl;
		//cout << "diffEpsilonEstimator " << diffEpsilonEstimator << endl;

		if ((diffEpsilonSketch > diffEpsilonEstimator) || (downsamplings_before_merge > 0))
		{
			--downsamplings_before_merge;
			
			++minus_log_p;
			LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
			divide_counters();

			for (int row = 0; row < height; ++row)
			{
				counter_values[row] >>= 1;
			}

			rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));

			//cout << "downsampling with probability: " << 1.0 / (1 << (minus_log_p)) << endl;
		}
		else
		{
			epsilonSketch *= 2;

			++largest_log_counter;

			//cout << "merging, largest_log_counter: " << (int)largest_log_counter << endl;
		}
	}
	increment_and_potentially_merge();

}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

FineGrainedSalsaAnalyticalErrorCMSSanity::FineGrainedSalsaAnalyticalErrorCMSSanity(int counterSize, int width, int height, int seed, int downsamplings_before_merge, bool split_counters_flag, double delta) :
	////
	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX)),
	RandByteArray((char *)&rand_bits)
	////
{

	// for sanity with AEE
	gen_arr();
	gen_arr();
	gen_arr();
	gen_arr();

	this->counterSize = counterSize;
	this->width = width;
	this->height = height;

	counters_per_word = 64 / counterSize;

	logCounterSize = 0;
	int index = counterSize;
	while (index >>= 1) ++logCounterSize;

	log_counters_per_word = 0;
	index = counters_per_word;
	while (index >>= 1) ++log_counters_per_word;

	width_mask = (width - 1) << (6 - log_counters_per_word);

	assert((counterSize & (counterSize - 1)) == 0 && "We assume that counterSize is a power of 2!");

	assert(width > 0 && "We assume too much!");
	assert(width % 64 == 0 && "We assume that (w % 64 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cms_counters = new uint64_t*[height];
	baseline_cms_merges = new uint64_t*[height];

	bobhash = new BOBHash[height];

	for (int row = 0; row < height; ++row)
	{
		baseline_cms_counters[row] = new uint64_t[width / counters_per_word]();
		baseline_cms_merges[row] = new uint64_t[width >> 6]();

		bobhash[row].initialize(seed*(7 + row) + row + 100);
	}

	/*
	for (int index_mod_16 = 0; index_mod_16 < 16; ++index_mod_16)
	{
		uint16_t mod_2_bit_mask = 1 << (index_mod_16 & 0b1110);
		uint16_t mod_4_bit_mask = 1 << ((index_mod_16 & 0b1100) | 0b1);
		uint16_t mod_8_bit_mask = 1 << ((index_mod_16 & 0b1000) | 0b11);

		merges_lookup_table_8[index_mod_16] = mod_2_bit_mask;
		merges_lookup_table_16[index_mod_16] = mod_4_bit_mask;
		merges_lookup_table_32[index_mod_16] = mod_8_bit_mask;
	}
	*/

	////
	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	counter_indices = new int[height];
	log_counter_sizes = new uint8_t[height];
	counter_values = new uint64_t[height];
	counter_bitSizes = new uint8_t[height];

	largest_log_counter = 0;
	this->downsamplings_before_merge = downsamplings_before_merge;

	this->split_counters_flag = split_counters_flag;
	////
	this->delta = delta;
	epsilonSketch = (pow(1.0 / delta, 1.0 / height)) / width;
	numPackets = 0;
}

FineGrainedSalsaAnalyticalErrorCMSSanity::~FineGrainedSalsaAnalyticalErrorCMSSanity()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms_counters[i];
		delete[] baseline_cms_merges[i];
	}

	delete[] bobhash;
	delete[] baseline_cms_counters;
	delete[] baseline_cms_merges;

	////
	delete[] counter_indices;
	delete[] counter_values;
	delete[] log_counter_sizes;
	delete[] counter_bitSizes;
	////
}

inline uint8_t FineGrainedSalsaAnalyticalErrorCMSSanity::log_counter_size(uint index, uint row)
{

	// 64 counters in each merge word  
	uint32_t word_index = index >> 6;

	// index % 64
	uint16_t index_mod_64 = index & 0b111111;

	//uint16_t mod_2_bit_mask = merges_lookup_table_8[index_mod_16];
	//uint16_t mod_4_bit_mask = merges_lookup_table_16[index_mod_16];
	//uint16_t mod_8_bit_mask = merges_lookup_table_32[index_mod_16];

	uint64_t mod_2_bit_mask = (uint64_t)1 << (index_mod_64 & 0b111110);
	uint64_t mod_4_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b111100) | 0b1);
	uint64_t mod_8_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b111000) | 0b11);
	uint64_t mod_16_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b110000) | 0b111);
	uint64_t mod_32_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b100000) | 0b1111);
	uint64_t mod_64_bit_mask = (uint64_t)1 << ((index_mod_64 & 0b000000) | 0b11111);

	uint64_t merges_word = baseline_cms_merges[row][word_index];

	bool mm2 = ((merges_word & mod_2_bit_mask) != 0);
	bool mm4 = ((merges_word & mod_4_bit_mask) != 0);
	bool mm8 = ((merges_word & mod_8_bit_mask) != 0);
	bool mm16 = ((merges_word & mod_16_bit_mask) != 0);
	bool mm32 = ((merges_word & mod_32_bit_mask) != 0);
	bool mm64 = ((merges_word & mod_64_bit_mask) != 0);

	return mm2 + mm4 + mm8 + mm16 + mm32 + mm64;

}

void FineGrainedSalsaAnalyticalErrorCMSSanity::merge_size_1_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_2_bit_mask = ((uint64_t)1 << (index & 0b111110));

	if ((baseline_cms_merges[row][word_index] & mod_2_bit_mask))
	{
		// counters already merged
	}
	else
	{
		baseline_cms_merges[row][word_index] |= mod_2_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::merge_size_2_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_4_bit_mask = ((uint64_t)1 << ((index & 0b111100) | 0b1));

	if ((baseline_cms_merges[row][word_index] & mod_4_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_1_index = index >> 1;
		uint paired_shifted_1_index = shifted_1_index ^ 0b1;

		// merge other counters if needed
		merge_size_1_counters(index, row);
		merge_size_1_counters(paired_shifted_1_index << 1, row);

		baseline_cms_merges[row][word_index] |= mod_4_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::merge_size_4_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_8_bit_mask = ((uint64_t)1 << ((index & 0b111000) | 0b11));

	if ((baseline_cms_merges[row][word_index] & mod_8_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_2_index = index >> 2;
		uint paired_shifted_2_index = shifted_2_index ^ 0b1;

		// merge other counters if needed
		merge_size_2_counters(index, row);
		merge_size_2_counters(paired_shifted_2_index << 2, row);

		baseline_cms_merges[row][word_index] |= mod_8_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::merge_size_8_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_16_bit_mask = ((uint64_t)1 << ((index & 0b110000) | 0b111));

	if ((baseline_cms_merges[row][word_index] & mod_16_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_3_index = index >> 3;
		uint paired_shifted_3_index = shifted_3_index ^ 0b1;

		// merge other counters if needed
		merge_size_4_counters(index, row);
		merge_size_4_counters(paired_shifted_3_index << 3, row);

		baseline_cms_merges[row][word_index] |= mod_16_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::merge_size_16_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_32_bit_mask = ((uint64_t)1 << ((index & 0b100000) | 0b1111));

	if ((baseline_cms_merges[row][word_index] & mod_32_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_4_index = index >> 4;
		uint paired_shifted_4_index = shifted_4_index ^ 0b1;

		// merge other counters if needed
		merge_size_8_counters(index, row);
		merge_size_8_counters(paired_shifted_4_index << 4, row);

		baseline_cms_merges[row][word_index] |= mod_32_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::merge_size_32_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_64_bit_mask = ((uint64_t)1 << ((index & 0b000000) | 0b11111));

	if ((baseline_cms_merges[row][word_index] & mod_64_bit_mask))
	{
		// counters already merged
	}
	else
	{
		uint shifted_5_index = index >> 5;
		uint paired_shifted_5_index = shifted_5_index ^ 0b1;

		// merge other counters if needed
		merge_size_16_counters(index, row);
		merge_size_16_counters(paired_shifted_5_index << 5, row);

		baseline_cms_merges[row][word_index] |= mod_64_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::set_non_overflow_counter(uint index, uint row, uint8_t log_counter_size_value, uint64_t new_value, int bitSize)
{
	int counter_word_index = index >> log_counters_per_word;
	int offset = ((index >> log_counter_size_value) << log_counter_size_value) & (counters_per_word - 1);
	int bitOffset = offset * counterSize;
	int negativeBitOffset = 64 - bitOffset - bitSize;
	uint64_t ones = ((uint64_t)1 << (counterSize*(1 << log_counter_size_value))) - 1;
	uint64_t counter_mask = ones << negativeBitOffset;

	baseline_cms_counters[row][counter_word_index] &= ~counter_mask;
	baseline_cms_counters[row][counter_word_index] |= (new_value << negativeBitOffset);
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::increment(const char * str)
{
	++numPackets;
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

void FineGrainedSalsaAnalyticalErrorCMSSanity::increment_h(const char * str)
{

	// remember all counter indices and the sizes
	bool is_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		counter_indices[row] = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
		counter_values[row] = read_counter(counter_indices[row], row);
		log_counter_sizes[row] = log_counter_size(counter_indices[row], row);
		counter_bitSizes[row] = (1 << ((6 - log_counters_per_word) + log_counter_sizes[row]));

		if (counter_values[row] == (((uint64_t)1 << counter_bitSizes[row]) - 1))
		{
			is_overflow = true;
		}
	}

	// safe to increment counters
	if (!is_overflow)
	{
		for (int row = 0; row < height; ++row)
		{
			set_non_overflow_counter(counter_indices[row], row, log_counter_sizes[row], counter_values[row] + 1, counter_bitSizes[row]);
		}
	}
	// NOT safe to increment counters
	else
	{
		handle_overflow_and_increment();
	}

}

void FineGrainedSalsaAnalyticalErrorCMSSanity::increment_and_potentially_merge()
{
	for (int row = 0; row < height; ++row)
	{

		int index = counter_indices[row];
		uint64_t counter_value = counter_values[row];
		uint8_t log_counter_size_value = log_counter_sizes[row];
		int bitSize = counter_bitSizes[row];

		bool is_overflow = (counter_value + 1) == ((uint64_t)1 << bitSize);

		if (!is_overflow)
		{
		}
		else
		{
			if (log_counter_size_value == 0)
			{
				merge_size_1_counters(index, row);
			}
			else if (log_counter_size_value == 1)
			{
				merge_size_2_counters(index, row);
			}
			else if (log_counter_size_value == 2)
			{
				merge_size_4_counters(index, row);
			}
			else if (log_counter_size_value == 3)
			{
				merge_size_8_counters(index, row);
			}
			else if (log_counter_size_value == 4)
			{
				merge_size_16_counters(index, row);
			}
			else if (log_counter_size_value == 5)
			{
				merge_size_32_counters(index, row);
			}
			++log_counter_size_value;
			bitSize <<= 1;
		}
		set_non_overflow_counter(index, row, log_counter_size_value, counter_value + 1, bitSize);

	}
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::handle_overflow_and_increment()
{

	bool is_max_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		if ((counter_values[row] == (((uint64_t)1 << counter_bitSizes[row]) - 1)) && (log_counter_sizes[row] == largest_log_counter))
		{
			is_max_overflow = true;
		}
	}

	if (!is_max_overflow) // merge
	{
	}
	else
	{
		double diffEpsilonSketch = epsilonSketch;

		double newNumPacketsPrime = (double)numPackets / (double)(1 << (minus_log_p + 1));
		double newEpsilonEstimator = sqrt((2.0*log(2.0 / delta)) / newNumPacketsPrime);

		double currentNumPacketsPrime = (double)numPackets / (double)(1 << minus_log_p);
		double currentEpsilonEstimator = sqrt((2.0*log(2.0 / delta)) / currentNumPacketsPrime);

		double diffEpsilonEstimator = newEpsilonEstimator - currentEpsilonEstimator;

		//cout << "newEpsilonEstimator " << newEpsilonEstimator << endl;
		//cout << "currentEpsilonEstimator " << currentEpsilonEstimator << endl;

		//cout << "diffEpsilonSketch " << diffEpsilonSketch << endl;
		//cout << "diffEpsilonEstimator " << diffEpsilonEstimator << endl;

		if ((diffEpsilonSketch > diffEpsilonEstimator) || (downsamplings_before_merge > 0))
		{
			--downsamplings_before_merge;

			++minus_log_p;
			LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
			divide_counters();

			for (int row = 0; row < height; ++row)
			{
				counter_values[row] >>= 1;
			}

			if (split_counters_flag)
			{
				split_counters();
			}

			rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));

			//cout << "downsampling with probability: " << 1.0 / (1 << (minus_log_p)) << endl;
		}
		else
		{
			epsilonSketch *= 2;

			++largest_log_counter;

			//cout << "merging, largest_log_counter: " << (int)largest_log_counter << endl;
		}
	}
	increment_and_potentially_merge();

}

inline uint64_t FineGrainedSalsaAnalyticalErrorCMSSanity::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);
	int counter_word_index = index >> log_counters_per_word;
	int offset = ((index >> log_counter_size_value) << log_counter_size_value) & (counters_per_word - 1);
	int bitSize = counterSize * ((uint64_t)1 << log_counter_size_value);
	int bitOffset = offset * counterSize;
	int negativeBitOffset = 64 - bitOffset - bitSize;
	uint64_t ones = ((uint64_t)1 << (counterSize*(1 << log_counter_size_value))) - 1;

	return (baseline_cms_counters[row][counter_word_index] >> negativeBitOffset) & ones;
}

uint64_t FineGrainedSalsaAnalyticalErrorCMSSanity::query(const char * str)
{
	int row = 0;

	uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return ((uint64_t)min) << minus_log_p;
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::divide_counters()
{
	//cout << "divide_counters" << endl;
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{

			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			uint64_t counter_value = read_counter(current_counter_index, row);
			int bitSize = counterSize * (1 << log_counter_size_value);

			set_non_overflow_counter(current_counter_index, row, log_counter_size_value, counter_value >> 1, bitSize);

			current_counter_index += (1 << log_counter_size_value);
		}
	}
}

void FineGrainedSalsaAnalyticalErrorCMSSanity::split_counters()
{
	//cout << "split_counters" << endl;
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{

			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			uint64_t counter_value = read_counter(current_counter_index, row);
			int bitSize = counterSize * (1 << log_counter_size_value);

			uint64_t threshold_value = (1 << (bitSize >> 1));
			if ((counter_value < threshold_value) && (log_counter_size_value > 0))
			{
				split_counter(row, current_counter_index, log_counter_size_value, counter_value);
			}
			current_counter_index += (1 << log_counter_size_value);
		}
	}

}

void FineGrainedSalsaAnalyticalErrorCMSSanity::split_counter(int row, int index, uint8_t log_counter_size_value, uint64_t counter_value)
{
	uint64_t bit_mask;

	switch (log_counter_size_value)
	{
	case 0:
		// minimal size counter - cannot split
		return;
	case 1:
		bit_mask = ((uint64_t)1 << (index & 0b111110));
		break;
	case 2:
		bit_mask = ((uint64_t)1 << ((index & 0b111100) | 0b1));
		break;
	case 3:
		bit_mask = ((uint64_t)1 << ((index & 0b111000) | 0b11));
		break;
	case 4:
		bit_mask = ((uint64_t)1 << ((index & 0b110000) | 0b111));
		break;
	case 5:
		bit_mask = ((uint64_t)1 << ((index & 0b100000) | 0b1111));
		break;
	case 6:
		bit_mask = ((uint64_t)1 << ((index & 0b000000) | 0b11111));
		break;
	default:
		cout << "Error in FineGrainedSalsaCMSSanity::split_counter" << endl;
		exit(1);
	}

	int word_index = (index >> 6);
	baseline_cms_merges[row][word_index] &= ~bit_mask;

	uint shifted_index = index >> (log_counter_size_value - 1);
	uint paired_shifted_index = shifted_index ^ 0b1;

	int bitSize = (1 << ((6 - log_counters_per_word) + (log_counter_size_value - 1)));

	set_non_overflow_counter(index, row, log_counter_size_value - 1, counter_value, bitSize);
	set_non_overflow_counter(paired_shifted_index << (log_counter_size_value - 1), row, log_counter_size_value - 1, counter_value, bitSize);
}







