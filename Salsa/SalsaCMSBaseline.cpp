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

#include "SalsaCMSBaseline.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

SalsaCMSBaseline::SalsaCMSBaseline()
{
}

SalsaCMSBaseline::~SalsaCMSBaseline()
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
}

void SalsaCMSBaseline::initialize(int width, int height, int seed)
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
		bobhash[row].initialize(seed*(7+row) + row + 100);
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
}

inline uint8_t SalsaCMSBaseline::log_counter_size(uint index, uint row)
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

void SalsaCMSBaseline::merge_8_bit_counters(uint index, uint row)
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
		uint16_t* p = (uint16_t*)&baseline_cms_counters[row][index & 0xFFFFFFFE];
		
		*p = baseline_cms_counters[row][index] + baseline_cms_counters[row][index ^ 0b1];
	}	
}

void SalsaCMSBaseline::merge_16_bit_counters(uint index, uint row)
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
		uint32_t* p = (uint32_t*)&baseline_cms_counters_16[row][index_16 & 0xFFFFFFFE];
		
		*p = baseline_cms_counters_16[row][index_16] + baseline_cms_counters_16[row][index_16 ^ 0b1];
	}
}

void SalsaCMSBaseline::merge_32_bit_counters(uint index, uint row)
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
		uint64_t* p = (uint64_t*)&baseline_cms_counters_16[row][index_32 & 0xFFFFFFFE];
		
		*p = baseline_cms_counters_32[row][index_32] + baseline_cms_counters_32[row][index_32 ^ 0b1];
	}
}

void SalsaCMSBaseline::merge_64_bit_counters(uint index, uint row)
{
	// maybe in the future...
}

void SalsaCMSBaseline::increment(const char * str)
{
	for (int row = 0; row < height; ++row) {
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint8_t log_counter_size_value = log_counter_size(index, row);

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
				baseline_cms_counters_16[row][index_16] += 1 << 8;
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
				baseline_cms_counters_32[row][index_32] += 1 << 16;
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
				baseline_cms_counters_64[row][index_64] += ((uint64_t) 1) << 32;
			}
		}
		else if (log_counter_size_value == 3)
		{
			uint index_64 = index >> 3;
			if (++baseline_cms_counters_64[row][index_64])
			{
				// no overflow
			}
			else
			{
				// merge_64_bit_counters - not implemented
				//merge_64_bit_counters(index, row);

				// the current counter became 0 while its real value is 2^64...
				//uint index_128 = index >> 4;
				//baseline_cms_counters_128[row][index_128] += 1 << 64;
			}
		}
	}
}

inline uint64_t SalsaCMSBaseline::read_counter(uint index, uint row)
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

uint64_t SalsaCMSBaseline::query(const char * str)
{
	int row = 0;

	uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

MaximumSalsaCMSBaseline::MaximumSalsaCMSBaseline()
{
}

MaximumSalsaCMSBaseline::~MaximumSalsaCMSBaseline()
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
}

void MaximumSalsaCMSBaseline::initialize(int width, int height, int seed)
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

}

inline uint8_t MaximumSalsaCMSBaseline::log_counter_size(uint index, uint row)
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

void MaximumSalsaCMSBaseline::merge_8_bit_counters(uint index, uint row)
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

void MaximumSalsaCMSBaseline::merge_16_bit_counters(uint index, uint row)
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

void MaximumSalsaCMSBaseline::merge_32_bit_counters(uint index, uint row)
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

void MaximumSalsaCMSBaseline::increment(const char * str)
{
	for (int row = 0; row < height; ++row) {
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint8_t log_counter_size_value = log_counter_size(index, row);

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

inline uint64_t MaximumSalsaCMSBaseline::read_counter(uint index, uint row)
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

uint64_t MaximumSalsaCMSBaseline::query(const char * str)
{
	int row = 0;

	uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

MaximumSalsaCMSBaselineSanity::MaximumSalsaCMSBaselineSanity()
{
}

MaximumSalsaCMSBaselineSanity::~MaximumSalsaCMSBaselineSanity()
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
}

void MaximumSalsaCMSBaselineSanity::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

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

}

inline uint8_t MaximumSalsaCMSBaselineSanity::log_counter_size(uint index, uint row)
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

void MaximumSalsaCMSBaselineSanity::merge_8_bit_counters(uint index, uint row)
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

void MaximumSalsaCMSBaselineSanity::merge_16_bit_counters(uint index, uint row)
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

void MaximumSalsaCMSBaselineSanity::merge_32_bit_counters(uint index, uint row)
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

void MaximumSalsaCMSBaselineSanity::increment(const char * str)
{
	for (int row = 0; row < height; ++row) {

		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> 3;
		uint8_t log_counter_size_value = log_counter_size(index, row);

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

inline uint64_t MaximumSalsaCMSBaselineSanity::read_counter(uint index, uint row)
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

uint64_t MaximumSalsaCMSBaselineSanity::query(const char * str)
{
	int row = 0;

	uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> 3;
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> 3;
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

FineGrainedSalsaCMSBaseline::FineGrainedSalsaCMSBaseline(int counterSize, int width, int height, int seed)
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

	width_mask = width - 1;

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

}

FineGrainedSalsaCMSBaseline::~FineGrainedSalsaCMSBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms_counters[i];
		delete[] baseline_cms_merges[i];
	}

	delete[] bobhash;
	delete[] baseline_cms_counters;
	delete[] baseline_cms_merges;
}

inline uint8_t FineGrainedSalsaCMSBaseline::log_counter_size(uint index, uint row)
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

void FineGrainedSalsaCMSBaseline::merge_size_1_counters(uint index, uint row)
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

void FineGrainedSalsaCMSBaseline::merge_size_2_counters(uint index, uint row)
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

void FineGrainedSalsaCMSBaseline::merge_size_4_counters(uint index, uint row)
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

void FineGrainedSalsaCMSBaseline::merge_size_8_counters(uint index, uint row)
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

void FineGrainedSalsaCMSBaseline::merge_size_16_counters(uint index, uint row)
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

void FineGrainedSalsaCMSBaseline::merge_size_32_counters(uint index, uint row)
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

void FineGrainedSalsaCMSBaseline::set_non_overflow_counter(uint index, uint row, uint8_t log_counter_size_value, uint64_t new_value, int bitSize)
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

void FineGrainedSalsaCMSBaseline::increment(const char * str)
{
	for (int row = 0; row < height; ++row)
	{

		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint8_t log_counter_size_value = log_counter_size(index, row);
		uint64_t counter_value = read_counter(index, row);
		int bitSize = counterSize * (1 << log_counter_size_value);
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

inline uint64_t FineGrainedSalsaCMSBaseline::read_counter(uint index, uint row)
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

uint64_t FineGrainedSalsaCMSBaseline::query(const char * str)
{
	int row = 0;

	uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

FineGrainedSalsaCMSBaselineSanity::FineGrainedSalsaCMSBaselineSanity(int counterSize, int width, int height, int seed)
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

	mergedZeroValuedCountersMatrix = new int*[height];
	for (int row = 0; row < height; ++row)
	{
		mergedZeroValuedCountersMatrix[row] = new int[6]();
	}

	for (int i = 0; i < 6; ++i)
	{
		mergedZeroValuedCounters[i] = 0;
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

}

FineGrainedSalsaCMSBaselineSanity::~FineGrainedSalsaCMSBaselineSanity()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms_counters[i];
		delete[] baseline_cms_merges[i];

		delete[] mergedZeroValuedCountersMatrix[i];
	}

	delete[] mergedZeroValuedCountersMatrix;

	delete[] bobhash;
	delete[] baseline_cms_counters;
	delete[] baseline_cms_merges;
}

inline uint8_t FineGrainedSalsaCMSBaselineSanity::log_counter_size(uint index, uint row)
{

	// 64 counters in each merge word  
	uint32_t word_index = index >> 6;

	// index % 64
	uint16_t index_mod_64 = index & 0b111111;

	//uint16_t mod_2_bit_mask = merges_lookup_table_8[index_mod_16];
	//uint16_t mod_4_bit_mask = merges_lookup_table_16[index_mod_16];
	//uint16_t mod_8_bit_mask = merges_lookup_table_32[index_mod_16];

	uint64_t mod_2_bit_mask  = (uint64_t)1 << ( index_mod_64 & 0b111110);
	uint64_t mod_4_bit_mask  = (uint64_t)1 << ((index_mod_64 & 0b111100) | 0b1);
	uint64_t mod_8_bit_mask  = (uint64_t)1 << ((index_mod_64 & 0b111000) | 0b11);
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

void FineGrainedSalsaCMSBaselineSanity::merge_size_1_counters(uint index, uint row)
{
	int word_index = (index >> 6);
	uint64_t mod_2_bit_mask = ((uint64_t)1 << (index & 0b111110));

	if ((baseline_cms_merges[row][word_index] & mod_2_bit_mask))
	{
		// counters already merged
	}
	else
	{
		// for count distinct
		if (read_counter(index, row) == 0)
		{
			++mergedZeroValuedCounters[0];
			++mergedZeroValuedCountersMatrix[row][0];
		}
		if (read_counter(index ^ 0b1, row) == 0)
		{
			++mergedZeroValuedCounters[0];
			++mergedZeroValuedCountersMatrix[row][0];
		}
		
		baseline_cms_merges[row][word_index] |= mod_2_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSBaselineSanity::merge_size_2_counters(uint index, uint row)
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

		// for count distinct
		if (read_counter(index, row) == 0)
		{
			++mergedZeroValuedCounters[1];
			++mergedZeroValuedCountersMatrix[row][1];
		}
		if (read_counter(paired_shifted_1_index << 1, row) == 0)
		{
			++mergedZeroValuedCounters[1];
			++mergedZeroValuedCountersMatrix[row][1];
		}

		baseline_cms_merges[row][word_index] |= mod_4_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSBaselineSanity::merge_size_4_counters(uint index, uint row)
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

		// for count distinct
		if (read_counter(index, row) == 0)
		{
			++mergedZeroValuedCounters[2];
			++mergedZeroValuedCountersMatrix[row][2];
		}
		if (read_counter(paired_shifted_2_index << 2, row) == 0)
		{
			++mergedZeroValuedCounters[2];
			++mergedZeroValuedCountersMatrix[row][2];
		}

		baseline_cms_merges[row][word_index] |= mod_8_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSBaselineSanity::merge_size_8_counters(uint index, uint row)
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

		// for count distinct
		if (read_counter(index, row) == 0)
		{
			++mergedZeroValuedCounters[3];
			++mergedZeroValuedCountersMatrix[row][3];
		}
		if (read_counter(paired_shifted_3_index << 3, row) == 0)
		{
			++mergedZeroValuedCounters[3];
			++mergedZeroValuedCountersMatrix[row][3];
		}

		baseline_cms_merges[row][word_index] |= mod_16_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSBaselineSanity::merge_size_16_counters(uint index, uint row)
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

		// for count distinct
		if (read_counter(index, row) == 0)
		{
			++mergedZeroValuedCounters[4];
			++mergedZeroValuedCountersMatrix[row][4];
		}
		if (read_counter(paired_shifted_4_index << 4, row) == 0)
		{
			++mergedZeroValuedCounters[4];
			++mergedZeroValuedCountersMatrix[row][4];
		}

		baseline_cms_merges[row][word_index] |= mod_32_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSBaselineSanity::merge_size_32_counters(uint index, uint row)
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

		// for count distinct
		if (read_counter(index, row) == 0)
		{
			++mergedZeroValuedCounters[5];
			++mergedZeroValuedCountersMatrix[row][5];
		}
		if (read_counter(paired_shifted_5_index << 5, row) == 0)
		{
			++mergedZeroValuedCounters[5];
			++mergedZeroValuedCountersMatrix[row][5];
		}

		baseline_cms_merges[row][word_index] |= mod_64_bit_mask;
		// value update (taking the maximum) happens in increment.
	}
}

void FineGrainedSalsaCMSBaselineSanity::set_non_overflow_counter(uint index, uint row, uint8_t log_counter_size_value, uint64_t new_value, int bitSize)
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

double FineGrainedSalsaCMSBaselineSanity::countDistinctRowEstimate(uint row)
{
	double UpperBoundArray[6] = { 0 };
	double LowerBoundArray[6] = { 0 };	
	
	double zeroValuedCounters[6] = { 0 };
	for (int currentLogSize = 0; currentLogSize < 6; ++currentLogSize)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			bool allZeroes = true;
			for (int currentOffset = 0; currentOffset < (1 << currentLogSize); ++currentOffset)
			{
				if (read_counter(current_counter_index + currentOffset, row) > 0)
				{
					allZeroes = false;
				}
			}
			zeroValuedCounters[currentLogSize] += allZeroes;
			current_counter_index += (1 << currentLogSize);
		}
		UpperBoundArray[currentLogSize] = zeroValuedCounters[currentLogSize] + mergedZeroValuedCountersMatrix[row][currentLogSize];
		LowerBoundArray[currentLogSize] = zeroValuedCounters[currentLogSize];
	}

	assert(UpperBoundArray[log_counters_per_word - 1] == LowerBoundArray[log_counters_per_word - 1] && "Bounds missmatch!");

	double currentWidth = width / (32.0 / counterSize);
	double p = UpperBoundArray[log_counters_per_word - 1] / currentWidth;
	double baselineEstimate = logl(p) / logl(1.0 - 1.0 / currentWidth);

	double currentEstimate = baselineEstimate;
	for (int i = log_counters_per_word - 2; i >= 0; --i)
	{

		//cout << "row " << row << " bounds:\t" << UpperBoundArray[i] << " " << LowerBoundArray[i] << " " << i << endl;

		double currentWidth = width >> i;
		double p_upper = UpperBoundArray[i] / currentWidth;
		double p_lower = LowerBoundArray[i] / currentWidth;
		double ithUpperEstimate = logl(p_lower) / logl(1.0 - 1.0 / currentWidth);
		double ithLowerEstimate = logl(p_upper) / logl(1.0 - 1.0 / currentWidth);

		//cout << "row " << row << " estimations:\t" << currentEstimate << "\t" << ithLowerEstimate << " " << ithUpperEstimate << endl;
		//cout << "row " << row << " currentWidth:\t" << currentWidth << endl;

		if (currentEstimate < ithLowerEstimate)
		{
			currentEstimate = ithLowerEstimate;
		}

		if (currentEstimate > ithUpperEstimate)
		{
			currentEstimate = ithUpperEstimate;
		}
	}
	//cout << "row " << row << " result: " << currentEstimate << endl << endl;
	return currentEstimate;

}

double FineGrainedSalsaCMSBaselineSanity::countDistinctRowEstimateAverage(uint row)
{
	double UpperBoundArray[6] = { 0 };
	double LowerBoundArray[6] = { 0 };

	double zeroValuedCounters[6] = { 0 };
	for (int currentLogSize = 0; currentLogSize < 6; ++currentLogSize)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			bool allZeroes = true;
			for (int currentOffset = 0; currentOffset < (1 << currentLogSize); ++currentOffset)
			{
				if (read_counter(current_counter_index + currentOffset, row) > 0)
				{
					allZeroes = false;
				}
			}
			zeroValuedCounters[currentLogSize] += allZeroes;
			current_counter_index += (1 << currentLogSize);
		}
		UpperBoundArray[currentLogSize] = zeroValuedCounters[currentLogSize] + mergedZeroValuedCountersMatrix[row][currentLogSize];
		LowerBoundArray[currentLogSize] = zeroValuedCounters[currentLogSize];
	}

	assert(UpperBoundArray[log_counters_per_word - 1] == LowerBoundArray[log_counters_per_word - 1] && "Bounds missmatch!");

	double currentWidth = width / (32.0 / counterSize);
	double p = UpperBoundArray[log_counters_per_word - 1] / currentWidth;
	double baselineEstimate = logl(p) / logl(1.0 - 1.0 / currentWidth);

	if (counterSize == 32)
	{
		return baselineEstimate;
	}

	currentWidth = width >> log_counters_per_word - 2;
	double p_upper = UpperBoundArray[log_counters_per_word - 2] / currentWidth;
	double p_lower = LowerBoundArray[log_counters_per_word - 2] / currentWidth;
	double ithUpperEstimate = logl(p_lower) / logl(1.0 - 1.0 / currentWidth);
	double ithLowerEstimate = logl(p_upper) / logl(1.0 - 1.0 / currentWidth);

	return sqrt(ithUpperEstimate*ithLowerEstimate);
	
	/*
	double currentEstimate = baselineEstimate;
	for (int i = log_counters_per_word - 2; i >= 0; --i)
	{
		double currentWidth = width >> i;
		double p_upper = UpperBoundArray[i] / currentWidth;
		double p_lower = LowerBoundArray[i] / currentWidth;
		double ithUpperEstimate = logl(p_lower) / logl(1.0 - 1.0 / currentWidth);
		double ithLowerEstimate = logl(p_upper) / logl(1.0 - 1.0 / currentWidth);

		cout << currentEstimate << "\t" << ithLowerEstimate << " " << ithUpperEstimate << endl;

		if (currentEstimate < ithLowerEstimate)
		{
			currentEstimate = ithLowerEstimate;
		}

		if (currentEstimate > ithUpperEstimate)
		{
			currentEstimate = ithUpperEstimate;
		}
	}

	return currentEstimate;
	*/
}

double FineGrainedSalsaCMSBaselineSanity::provable_query_count_distinct()
{

	double* estimationsArray = new double[height];

	for (int row = 0; row < height; ++row)
	{
		estimationsArray[row] = countDistinctRowEstimate(row);
	}

	sort(estimationsArray, estimationsArray + height);

	//for (int row = 0; row < height; ++row)
	//{
	//	cout << estimationsArray[row] << "\t";
	//}
	//cout << endl;

	double finalEstimate;

	if (height % 2 == 0)
	{
		finalEstimate = 0.5*(estimationsArray[height >> 1] + estimationsArray[(height >> 1) - 1]);
	}
	else
	{
		finalEstimate = estimationsArray[height >> 1];
	}
	delete[] estimationsArray;

	return finalEstimate;

}

double FineGrainedSalsaCMSBaselineSanity::bound_average_query_count_distinct()
{
	double* estimationsArray = new double[height];

	for (int row = 0; row < height; ++row)
	{
		estimationsArray[row] = countDistinctRowEstimateAverage(row);
	}

	sort(estimationsArray, estimationsArray + height);

	//for (int row = 0; row < height; ++row)
	//{
	//	cout << estimationsArray[row] << "\t";
	//}
	//cout << endl;

	double finalEstimate;

	if (height % 2 == 0)
	{
		finalEstimate = 0.5*(estimationsArray[height >> 1] + estimationsArray[(height >> 1) - 1]);
	}
	else
	{
		finalEstimate = estimationsArray[height >> 1];
	}
	delete[] estimationsArray;

	return finalEstimate;
}

double FineGrainedSalsaCMSBaselineSanity::global_query_count_distinct()
{

	double* estimationsArray = new double[height];

	double UpperBoundArray[6] = { 0 };
	double LowerBoundArray[6] = { 0 };

	double zeroValuedCounters[6] = { 0 };
	for (int currentLogSize = 0; currentLogSize < 6; ++currentLogSize)
	{
		for (int row = 0; row < height; ++row)
		{
			int current_counter_index = 0;
			while (current_counter_index < width)
			{
				bool allZeroes = true;
				for (int currentOffset = 0; currentOffset < (1 << currentLogSize); ++currentOffset)
				{
					if (read_counter(current_counter_index + currentOffset, row) > 0)
					{
						allZeroes = false;
					}
				}
				zeroValuedCounters[currentLogSize] += allZeroes;
				current_counter_index += (1 << currentLogSize);
			}
		}
		UpperBoundArray[currentLogSize] = zeroValuedCounters[currentLogSize] + mergedZeroValuedCounters[currentLogSize];
		LowerBoundArray[currentLogSize] = zeroValuedCounters[currentLogSize];
	}

	assert(UpperBoundArray[log_counters_per_word - 1] == LowerBoundArray[log_counters_per_word - 1] && "Bounds missmatch!");

	double currentWidth = (height * width) / (32.0 / counterSize);
	double p = UpperBoundArray[log_counters_per_word - 1] / currentWidth;
	double baselineEstimate = logl(p) / logl(1.0 - 1.0 / currentWidth);

	double currentEstimate = baselineEstimate;
	for (int i = log_counters_per_word - 2; i >= 0; --i)
	{
		currentWidth = (height * width) >> i;
		double p_upper = UpperBoundArray[i] / currentWidth;
		double p_lower = LowerBoundArray[i] / currentWidth;
		double ithUpperEstimate = logl(p_lower) / logl(1.0 - 1.0 / currentWidth);
		double ithLowerEstimate = logl(p_upper) / logl(1.0 - 1.0 / currentWidth);

		//cout << currentEstimate << "\t" << ithLowerEstimate << " " << ithUpperEstimate << endl;

		if (currentEstimate < ithLowerEstimate)
		{
			currentEstimate = ithLowerEstimate;
		}

		if (currentEstimate > ithUpperEstimate)
		{
			currentEstimate = ithUpperEstimate;
		}
	}

	return currentEstimate / height;

}

double FineGrainedSalsaCMSBaselineSanity::best_guess_query_count_distinct()
{
	
	double* estimationsArray = new double[height];

	for (int row = 0; row < height; ++row)
	{


		double TotalMinSizeBuckets = 0;
		double NonZeroMinSizeBuckets = 0;

		int current_counter_index = 0;
		while (current_counter_index < width)
		{

			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			uint64_t counter_value = read_counter(current_counter_index, row);

			if (log_counter_size_value == 0)
			{
				++TotalMinSizeBuckets;
				if (counter_value > 0)
				{
					++NonZeroMinSizeBuckets;
				}
			}

			current_counter_index += (1 << log_counter_size_value);
		}

		double hitMinSizeFraction = NonZeroMinSizeBuckets / TotalMinSizeBuckets;
		double hitWidth = NonZeroMinSizeBuckets;

		current_counter_index = 0;
		while (current_counter_index < width)
		{

			uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
			uint64_t counter_value = read_counter(current_counter_index, row);

			if (log_counter_size_value > 0)
			{
				hitWidth += 1 + hitMinSizeFraction * ((1 << log_counter_size_value) - 1);
			}

			current_counter_index += (1 << log_counter_size_value);
		}

		double p = 1.0 - hitWidth / width;
		estimationsArray[row] = logl(p) / logl(1.0 - 1.0 / width);
	}
		
	sort(estimationsArray, estimationsArray + height);

	//for (int row = 0; row < height; ++row)
	//{
	//	cout << estimationsArray[row] << "\t";
	//}
	//cout << endl;
		
	double finalEstimate;

	if (height % 2 == 0)
	{
		finalEstimate = 0.5*(estimationsArray[height >> 1] + estimationsArray[(height >> 1) - 1]);
	}
	else
	{
		finalEstimate = estimationsArray[height >> 1];
	}
	delete[] estimationsArray;

	return finalEstimate;
	
}

double FineGrainedSalsaCMSBaselineSanity::global_best_guess_query_count_distinct()
{

	double UpperBoundArray[6] = { 0 };
	double LowerBoundArray[6] = { 0 };

	double zeroValuedCounters[6] = { 0 };
	for (int currentLogSize = 0; currentLogSize < 6; ++currentLogSize)
	{
		for (int row = 0; row < height; ++row)
		{
			int current_counter_index = 0;
			while (current_counter_index < width)
			{
				bool allZeroes = true;
				for (int currentOffset = 0; currentOffset < (1 << currentLogSize); ++currentOffset)
				{
					if (read_counter(current_counter_index + currentOffset, row) > 0)
					{
						allZeroes = false;
					}
				}
				zeroValuedCounters[currentLogSize] += allZeroes;
				current_counter_index += (1 << currentLogSize);
			}
		}
		UpperBoundArray[currentLogSize] = zeroValuedCounters[currentLogSize] + mergedZeroValuedCounters[currentLogSize];
		LowerBoundArray[currentLogSize] = zeroValuedCounters[currentLogSize];
	}

	assert(UpperBoundArray[log_counters_per_word - 1] == LowerBoundArray[log_counters_per_word - 1] && "Bounds missmatch!");

	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

	double estimationsArray[6] = { 0 };

	double TotalLogSizeEqualsIBuckets[6] = { 0 };
	double NonZeroLogSizeEqualsIBuckets[6] = { 0 };

	for (int currentLogSize = 0; currentLogSize < 6; ++currentLogSize)
	{
		for (int row = 0; row < height; ++row)
		{
			int current_counter_index = 0;
			while (current_counter_index < width)
			{
				bool allZeroes = true;
				for (int currentOffset = 0; currentOffset < (1 << currentLogSize); ++currentOffset)
				{
					if (read_counter(current_counter_index + currentOffset, row) > 0)
					{
						allZeroes = false;
					}
				}				
				if (log_counter_size(current_counter_index, row) <= currentLogSize)
				{
					NonZeroLogSizeEqualsIBuckets[currentLogSize] += (1 - allZeroes);
					++TotalLogSizeEqualsIBuckets[currentLogSize];
				}
				current_counter_index += (1 << currentLogSize);
			}
		}

		double hitLogSizeEqualsIFraction = NonZeroLogSizeEqualsIBuckets[currentLogSize] / TotalLogSizeEqualsIBuckets[currentLogSize];
		double hitLogSizeEqualsIWidth = NonZeroLogSizeEqualsIBuckets[currentLogSize];

		for (int row = 0; row < height; ++row)
		{
			int current_counter_index = 0;
			while (current_counter_index < width)
			{

				uint8_t log_counter_size_value = log_counter_size(current_counter_index, row);
				uint64_t counter_value = read_counter(current_counter_index, row);

				if (log_counter_size_value > currentLogSize)
				{
					hitLogSizeEqualsIWidth += 1 + hitLogSizeEqualsIFraction * ((1 << (log_counter_size_value - currentLogSize)) - 1);
				}

				current_counter_index += (1 << log_counter_size_value);
			}
		}
		double currentWidth = (height * width) >> currentLogSize;

		double p = 1.0 - hitLogSizeEqualsIWidth / currentWidth;
		estimationsArray[currentLogSize] = logl(p) / logl(1.0 - 1.0 / currentWidth);

	}

	////////////////////////////////////////////////////
	////////////////////////////////////////////////////

	double currentWidth = (height * width) / (32.0 / counterSize);
	double p = UpperBoundArray[log_counters_per_word - 1] / currentWidth;
	double baselineEstimate = logl(p) / logl(1.0 - 1.0 / currentWidth);

	double currentEstimate = baselineEstimate;
	for (int i = log_counters_per_word - 2; i >= 0; --i)
	{
		currentWidth = (height * width) >> i;
		double p_upper = UpperBoundArray[i] / currentWidth;
		double p_lower = LowerBoundArray[i] / currentWidth;
		double ithUpperEstimate = logl(p_lower) / logl(1.0 - 1.0 / currentWidth);
		double ithLowerEstimate = logl(p_upper) / logl(1.0 - 1.0 / currentWidth);

		//cout << currentEstimate << "\t" << ithLowerEstimate << " " << ithUpperEstimate << endl;

		if (currentEstimate < ithLowerEstimate)
		{
			currentEstimate = estimationsArray[i];
		}

		if (currentEstimate > ithUpperEstimate)
		{
			currentEstimate = estimationsArray[i];
		}
	}

	return currentEstimate / height;

}

void FineGrainedSalsaCMSBaselineSanity::increment(const char * str)
{
	for (int row = 0; row < height; ++row) 
	{

		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);
		uint8_t log_counter_size_value = log_counter_size(index, row);
		uint64_t counter_value = read_counter(index, row);
		int bitSize = counterSize * (1 << log_counter_size_value);
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

inline uint64_t FineGrainedSalsaCMSBaselineSanity::read_counter(uint index, uint row)
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

uint64_t FineGrainedSalsaCMSBaselineSanity::query(const char * str)
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
	return min;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

MaximumSalsaCUSBaseline::MaximumSalsaCUSBaseline()
{
}

MaximumSalsaCUSBaseline::~MaximumSalsaCUSBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cus_counters[i];
		delete[] baseline_cus_merges[i];
	}

	delete[] bobhash;
	delete[] baseline_cus_counters;

	delete[] baseline_cus_counters_16;
	delete[] baseline_cus_counters_32;
	delete[] baseline_cus_counters_64;

	delete[] baseline_cus_merges;
}

void MaximumSalsaCUSBaseline::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 16 == 0 && "We assume that (w % 16 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cus_counters = new uint8_t*[height];

	baseline_cus_counters_16 = new uint16_t*[height];
	baseline_cus_counters_32 = new uint32_t*[height];
	baseline_cus_counters_64 = new uint64_t*[height];

	baseline_cus_merges = new uint16_t*[height];
	bobhash = new BOBHash[height];

	incremet_indices = new uint[height];
	incremet_values = new uint64_t[height];

	for (int row = 0; row < height; ++row)
	{
		baseline_cus_counters[row] = new uint8_t[width]();

		baseline_cus_counters_16[row] = (uint16_t*)baseline_cus_counters[row];
		baseline_cus_counters_32[row] = (uint32_t*)baseline_cus_counters[row];
		baseline_cus_counters_64[row] = (uint64_t*)baseline_cus_counters[row];

		baseline_cus_merges[row] = new uint16_t[width >> 4]();
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

}

inline uint8_t MaximumSalsaCUSBaseline::log_counter_size(uint index, uint row)
{

	// 16 counters in each index  
	uint32_t byte_index = index >> 4;

	// index % 16
	uint16_t index_mod_16 = index & 0b1111;

	uint16_t mod_2_bit_mask = merges_lookup_table_8[index_mod_16]; // 1 << (index_mod_16 & 0b1110);
	uint16_t mod_4_bit_mask = merges_lookup_table_16[index_mod_16]; // 1 << ((index_mod_16 & 0b1100) | 0b1);
	uint16_t mod_8_bit_mask = merges_lookup_table_32[index_mod_16]; // 1 << ((index_mod_16 & 0b1000) | 0b11);

	uint16_t merges_word = baseline_cus_merges[row][byte_index];

	//uint16_t mm2 = merges_word & mod_2_bit_mask;
	//uint16_t mm4 = merges_word & mod_4_bit_mask;
	//uint16_t mm8 = merges_word & mod_8_bit_mask;

	return (merges_word & mod_2_bit_mask) != 0 ? ((merges_word & mod_4_bit_mask) == 0 ? 1 : (merges_word & mod_8_bit_mask) == 0 ? 2 : 3) : 0;

	// we assume that counter sizes cannot exceed 8X the original
	// if the original is 8 bits then we have a 64-bit limit

}

void MaximumSalsaCUSBaseline::merge_8_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_2_bit_mask = 1 << (index & 0b1110);

	if ((baseline_cus_merges[row][byte_index] & mod_2_bit_mask))
	{
		// counters already merged
	}
	else
	{
		baseline_cus_merges[row][byte_index] |= mod_2_bit_mask;

		// happens in increment
		//uint16_t* p = (uint16_t*)&baseline_cus_counters[row][index & 0xFFFFFFFE];
		//*p = max(baseline_cus_counters[row][index], baseline_cus_counters[row][index ^ 0b1]);
	}
}

void MaximumSalsaCUSBaseline::merge_16_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_4_bit_mask = 1 << ((index & 0b1100) | 0b1);

	if ((baseline_cus_merges[row][byte_index] & mod_4_bit_mask))
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

		baseline_cus_merges[row][byte_index] |= mod_4_bit_mask;

		// happens in increment
		//uint32_t* p = (uint32_t*)&baseline_cus_counters_16[row][index_16 & 0xFFFFFFFE];		
		//*p = max(baseline_cus_counters_16[row][index_16], baseline_cus_counters_16[row][index_16 ^ 0b1]);
	}
}

void MaximumSalsaCUSBaseline::merge_32_bit_counters(uint index, uint row)
{
	int byte_index = index >> 4;
	uint16_t mod_8_bit_mask = 1 << ((index & 0b1000) | 0b11);

	if ((baseline_cus_merges[row][byte_index] & mod_8_bit_mask))
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

		baseline_cus_merges[row][byte_index] |= mod_8_bit_mask;

		// happens in increment
		//uint64_t* p = (uint64_t*)&baseline_cus_counters_16[row][index_32 & 0xFFFFFFFE];		
		//*p = max(baseline_cus_counters_32[row][index_32], baseline_cus_counters_32[row][index_32 ^ 0b1]);
	}
}

void MaximumSalsaCUSBaseline::increment(const char * str)
{
	uint64_t minimum = -1;

	for (int row = 0; row < height; ++row) {

		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;

		incremet_values[row] = read_counter(index, row);
		incremet_indices[row] = index;
		minimum = (minimum < incremet_values[row]) ? minimum : incremet_values[row];
	}

	if (minimum < (1 << 8))
	{
		for (int row = 0; row < height; ++row)
		{
			if (minimum != incremet_values[row])
			{

			}
			else
			{
				uint index = incremet_indices[row];

				if (++baseline_cus_counters[row][index])
				{
					// no overflow
				}
				else
				{
					merge_8_bit_counters(index, row);

					// the current counter became 0 while its real value is 256...
					uint index_16 = index >> 1;
					baseline_cus_counters_16[row][index_16] = 1 << 8;
				}
			}
		}
	}
	else if (minimum < (1 << 16))
	{
		for (int row = 0; row < height; ++row)
		{
			if (minimum != incremet_values[row])
			{

			}
			else
			{
				uint index = incremet_indices[row];

				uint index_16 = index >> 1;
				if (++baseline_cus_counters_16[row][index_16])
				{
					// no overflow
				}
				else
				{
					merge_16_bit_counters(index, row);

					// the current counter became 0 while its real value is 2^16...
					uint index_32 = index >> 2;
					baseline_cus_counters_32[row][index_32] = 1 << 16;
				}
			}
		}
	}
	else if (minimum < ((uint64_t)1 << 32))
	{
		for (int row = 0; row < height; ++row)
		{
			if (minimum != incremet_values[row])
			{

			}
			else
			{
				uint index = incremet_indices[row];

				uint index_32 = index >> 2;
				if (++baseline_cus_counters_32[row][index_32])
				{
					// no overflow
				}
				else
				{
					merge_32_bit_counters(index, row);

					// the current counter became 0 while its real value is 2^32...
					uint index_64 = index >> 3;
					baseline_cus_counters_64[row][index_64] = ((uint64_t)1) << 32;
				}
			}
		}
	}
	else
	{
		for (int row = 0; row < height; ++row)
		{
			if (minimum != incremet_values[row])
			{

			}
			else
			{
				uint index = incremet_indices[row];

				uint index_64 = index >> 3;
				++baseline_cus_counters_64[row][index_64];
			}
		}
	}

}

inline uint64_t MaximumSalsaCUSBaseline::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);

	switch (log_counter_size_value) {
	case 0:
		return baseline_cus_counters[row][index];
	case 1:
		return baseline_cus_counters_16[row][index >> 1];
	case 2:
		return baseline_cus_counters_32[row][index >> 2];
	case 3:
		return baseline_cus_counters_64[row][index >> 3];
	case 4:
		// not implemented
		return -1;
	default:
		// BOB!
		return -1;
	}
}

uint64_t MaximumSalsaCUSBaseline::query(const char * str)
{
	int row = 0;

	uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
	uint64_t min = read_counter(index, row++);

	while (row < height) {
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}


