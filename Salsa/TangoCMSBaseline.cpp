#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "TangoCMSBaseline.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

bool isBitSet(uint64_t merges_word, uint8_t bitIndex)
{
	return !!(merges_word & (((uint64_t)1) << bitIndex));
}

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

TangoBaseline::TangoBaseline()
{
}

TangoBaseline::~TangoBaseline()
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

void TangoBaseline::initialize(int counterSize, int width, int height, int seed)
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
	//assert((counterSize <= 8) && "We assume that counterSize <= 8");

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
}

inline pair<uint8_t, uint8_t> TangoBaseline::counter_size_and_offset(uint index, uint row)
{

	// 64 counters in each index  
	uint64_t byte_index = index >> 6;

	// in-word index
	uint8_t index_mod_64 = index & 0b111111;

	uint64_t merges_word = baseline_cms_merges[row][byte_index];

	int size = 1;
	uint8_t bitIndex = index_mod_64;

	while ((bitIndex < 64) && (isBitSet(merges_word, bitIndex++)))
	{
		++size;
	}

	uint8_t offset = index_mod_64 - 1;

	while ((offset >= 0) && (isBitSet(merges_word, offset)))
	{
		++size;
		--offset;
	}

	return make_pair(size, (offset + 1) & (counters_per_word - 1));

}

void TangoBaseline::increment(const char * str)
{
	for (int row = 0; row < height; ++row) {
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;

		pair<uint8_t, uint8_t> size_and_offset = counter_size_and_offset(index, row);

		int size = size_and_offset.first;
		int offset = size_and_offset.second;

		uint64_t counter_value = read_counter(index, row);

		if (counter_value < (((uint64_t)1 << (size*counterSize))) - 1)
		{
			// no overflow
			int bitSize = size * counterSize;
			int bitOffset = offset * counterSize;
			int negativeOffset = 64 - bitSize - bitOffset;
			int wordOffset = index >> log_counters_per_word;

			uint64_t& cnt = baseline_cms_counters[row][wordOffset];

			uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
			uint64_t cnt_masked = cnt & cnt_mask;
			uint64_t cnt_shifted = (counter_value + 1) << negativeOffset;

			cnt = cnt_shifted | cnt_masked;

		}
		else
		{
			/*
			if (size == 1) // merging two size-1 counters
			{
				int bitSize = counterSize << 1;
				int index_pairity = index & 0b1;
				int even_index = index - index_pairity;
				int bitOffset = (offset  - index_pairity) * counterSize;
				int negativeOffset = 64 - bitSize - bitOffset;
				int wordOffset = even_index >> log_counters_per_word;

				// write to counter
				uint64_t& cnt = baseline_cms_counters[row][wordOffset];

				uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
				uint64_t cnt_masked = cnt & cnt_mask;
				uint64_t cnt_shifted = ((uint64_t)1 << counterSize) << negativeOffset;

				cnt = cnt_shifted | cnt_masked;

				// turn on bit
				uint64_t byte_index = even_index >> 6;

				uint8_t bitIndex = even_index & 0b111111;

				uint64_t& merges_word = baseline_cms_merges[row][byte_index];

				merges_word = merges_word | ((uint64_t)1 << bitIndex);
			}
			else
			{
			*/
			int log_size = 0;
			int i = size;
			while (i >>= 1) ++log_size;

			int in_word_offset = index & ((1 << log_counters_per_word) - 1);

			////

			bool merge_right;
			bool found_available_counter = false;

			int merge_bounries[6];

			// what the fuck? evil bitwise level hacking
			merge_bounries[0] = (index & 0b1) == 0 ? index | 0b1 : index & 0xFFFFFFFE;
			merge_bounries[1] = (index & 0b10) == 0 ? index | 0b11 : index & 0xFFFFFFFC;
			merge_bounries[2] = (index & 0b100) == 0 ? index | 0b111 : index & 0xFFFFFFF8;
			merge_bounries[3] = (index & 0b1000) == 0 ? index | 0b1111 : index & 0xFFFFFFF0;
			merge_bounries[4] = (index & 0b10000) == 0 ? index | 0b11111 : index & 0xFFFFFFE0;
			merge_bounries[5] = (index & 0b100000) == 0 ? index | 0b111111 : index & 0xFFFFFFC0;

			uint64_t& curr_merges_word = baseline_cms_merges[row][index >> 6];

			for (int j = 0; (j < 6) && (!(found_available_counter)); ++j)
			{
				if (index > merge_bounries[j])
				{
					for (int current_index = index - 1; current_index >= merge_bounries[j]; --current_index)
					{
						uint8_t bitIndex = current_index & 0b111111;
						if (!(isBitSet(curr_merges_word, bitIndex)))
						{
							found_available_counter = true;
							merge_right = false;
							break;
						}
					}
				}
				else
				{
					for (int current_index = index; current_index < merge_bounries[j]; ++current_index)
					{
						uint8_t bitIndex = current_index & 0b111111;
						if (!(isBitSet(curr_merges_word, bitIndex)))
						{
							found_available_counter = true;
							merge_right = true;
							break;
						}
					}
				}
			}

			////

			if (merge_right)
			{
				int paired_counter_index = (index - in_word_offset) + (size + offset);

				pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

				int paired_size = paired_size_and_offset.first;
				int paired_offset = paired_size_and_offset.second;

				uint64_t paired_counter_value = read_counter(paired_counter_index, row);

				uint64_t max_value = paired_counter_value > (counter_value + 1) ? paired_counter_value : counter_value + 1;

				////
				int bitSize = counterSize * (size + paired_size);
				//int index_pairity = index & 0b1;
				//int even_index = index - index_pairity;
				int bitOffset = offset * counterSize;
				int negativeOffset = 64 - bitSize - bitOffset;
				int wordOffset = index >> log_counters_per_word;

				// write to counter
				uint64_t& cnt = baseline_cms_counters[row][wordOffset];

				uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
				uint64_t cnt_masked = cnt & cnt_mask;
				uint64_t cnt_shifted = max_value << negativeOffset;

				cnt = cnt_shifted | cnt_masked;

				// turn on bit
				uint64_t byte_index = index >> 6;

				uint8_t bitIndex = (paired_counter_index - 1) & 0b111111;

				uint64_t& merges_word = baseline_cms_merges[row][byte_index];

				merges_word = merges_word | ((uint64_t)1 << bitIndex);
				////
			}
			else // merge left
			{
				int paired_counter_index = (index - in_word_offset) + (offset - 1);

				pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

				int paired_size = paired_size_and_offset.first;
				int paired_offset = paired_size_and_offset.second;

				uint64_t paired_counter_value = read_counter(paired_counter_index, row);

				uint64_t max_value = paired_counter_value > (counter_value + 1) ? paired_counter_value : counter_value + 1;

				////
				int bitSize = counterSize * (size + paired_size);
				//int index_pairity = index & 0b1;
				//int even_index = index - index_pairity;
				int bitOffset = paired_offset * counterSize;
				int negativeOffset = 64 - bitSize - bitOffset;
				int wordOffset = index >> log_counters_per_word;

				// write to counter
				uint64_t& cnt = baseline_cms_counters[row][wordOffset];

				uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
				uint64_t cnt_masked = cnt & cnt_mask;
				uint64_t cnt_shifted = max_value << negativeOffset;

				cnt = cnt_shifted | cnt_masked;

				// turn on bit
				uint64_t byte_index = index >> 6;

				uint8_t bitIndex = paired_counter_index & 0b111111;

				uint64_t& merges_word = baseline_cms_merges[row][byte_index];

				merges_word = merges_word | ((uint64_t)1 << bitIndex);
				////

			}

			// overflow
			//cout << "No more bob! " << counter_value << endl;
		//}

		}


	}
}

inline uint64_t TangoBaseline::read_counter(uint index, uint row)
{

	pair<uint8_t, uint8_t> size_and_offset = counter_size_and_offset(index, row);

	int size = size_and_offset.first;
	int offset = size_and_offset.second;

	int bitSize = size * counterSize;
	int bitOffset = offset * counterSize;
	int negativeOffset = 64 - bitSize - bitOffset;
	int wordOffset = index >> log_counters_per_word;

	return  (baseline_cms_counters[row][wordOffset] << bitOffset) >> (negativeOffset + bitOffset);

}

uint64_t TangoBaseline::query(const char * str)
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

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

smarTangoBaselineSanity::smarTangoBaselineSanity()
{
}

smarTangoBaselineSanity::~smarTangoBaselineSanity()
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

void smarTangoBaselineSanity::initialize(int counterSize, int width, int height, int seed)
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
}

inline pair<uint8_t, uint8_t> smarTangoBaselineSanity::counter_size_and_offset(uint index, uint row)
{

	// 64 counters in each index  
	uint64_t byte_index = index >> 6;

	// in-word index
	uint8_t index_mod_64 = index & 0b111111;

	uint64_t merges_word = baseline_cms_merges[row][byte_index];

	int size = 1;
	uint8_t bitIndex = index_mod_64;

	while ((bitIndex < 64) && (isBitSet(merges_word, bitIndex++)))
	{
		++size;
	}

	uint8_t offset = index_mod_64 - 1;

	while ((offset >= 0) && (isBitSet(merges_word, offset)))
	{
		++size;
		--offset;
	}

	return make_pair(size, (offset + 1) & (counters_per_word - 1));

}

void smarTangoBaselineSanity::increment(const char * str)
{
	for (int row = 0; row < height; ++row) {
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);

		pair<uint8_t, uint8_t> size_and_offset = counter_size_and_offset(index, row);

		int size = size_and_offset.first;
		int offset = size_and_offset.second;

		uint64_t counter_value = read_counter(index, row);

		int bitSize = size * counterSize;
		int bitOffset = offset * counterSize;
		int negativeOffset = 64 - bitSize - bitOffset;

		if (counter_value < (((uint64_t)1 << (size*counterSize))) - 1)
		{
			// no overflow

			int wordOffset = index >> log_counters_per_word;

			uint64_t& cnt = baseline_cms_counters[row][wordOffset];

			uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
			uint64_t cnt_masked = cnt & cnt_mask;
			uint64_t cnt_shifted = (counter_value + 1) << negativeOffset;

			cnt = cnt_shifted | cnt_masked;

		}
		else
		{

			int log_size = 0;
			int i = size;
			while (i >>= 1) ++log_size;

			int in_word_offset = index & ((1 << log_counters_per_word) - 1);

			bool merge_right;
			
			if ((offset == 0) || (offset == 32))
			{
				merge_right = true;
			}
			else if ((negativeOffset == 0) || (negativeOffset == 32))
			{
				merge_right = false;
			}
			else
			{
				int left_index = (index - in_word_offset) + (offset - 1);
				int right_index = (index - in_word_offset) + (size + offset);

				uint64_t left_value = read_counter(left_index, row);
				uint64_t right_value = read_counter(right_index, row);
				
				/*
				if (left_value > right_value)
				{
					merge_right = true;
				}
				else
				{
					merge_right = false;
				}
				*/
								
				pair<uint8_t, uint8_t> left_size_and_offset = counter_size_and_offset(left_index, row);
				pair<uint8_t, uint8_t> right_size_and_offset = counter_size_and_offset(right_index, row);

				int left_size = left_size_and_offset.first;
				int left_offset = left_size_and_offset.second;

				int right_size = right_size_and_offset.first;
				int right_offset = right_size_and_offset.second;
				
				if (left_size == right_size)
				{
					if (abs((long long)((counter_value + 1) - left_value)) > abs((long long)((counter_value + 1) - right_value)))
					{
						merge_right = true;
					}
					else
					{
						merge_right = false;
					}
				}
				else
				{
					if (left_size > right_size)
					{
						merge_right = true;
					}
					else
					{
						merge_right = false;
					}
				}
				

			}


			/*
			bool found_available_counter = false;
	
			int merge_bounries[6];

			merge_bounries[0] = (index & 0b1) == 0 ? index | 0b1 : index & 0xFFFFFFFE;
			merge_bounries[1] = (index & 0b10) == 0 ? index | 0b11 : index & 0xFFFFFFFC;
			merge_bounries[2] = (index & 0b100) == 0 ? index | 0b111 : index & 0xFFFFFFF8;
			merge_bounries[3] = (index & 0b1000) == 0 ? index | 0b1111 : index & 0xFFFFFFF0;
			merge_bounries[4] = (index & 0b10000) == 0 ? index | 0b11111 : index & 0xFFFFFFE0;
			merge_bounries[5] = (index & 0b100000) == 0 ? index | 0b111111 : index & 0xFFFFFFC0;

			uint64_t& curr_merges_word = baseline_cms_merges[row][index >> 6];

			for (int j = 0; (j < 6) && (!(found_available_counter)); ++j)
			{
				if (index > merge_bounries[j])
				{
					for (int current_index = index - 1; current_index >= merge_bounries[j]; --current_index)
					{
						uint8_t bitIndex = current_index & 0b111111;
						if (!(isBitSet(curr_merges_word, bitIndex)))
						{
							found_available_counter = true;
							merge_right = false;
							break;
						}
					}
				}
				else
				{
					for (int current_index = index; current_index < merge_bounries[j]; ++current_index)
					{
						uint8_t bitIndex = current_index & 0b111111;
						if (!(isBitSet(curr_merges_word, bitIndex)))
						{
							found_available_counter = true;
							merge_right = true;
							break;
						}
					}
				}
									
			}
			*/

			if (merge_right)
			{
				int paired_counter_index = (index - in_word_offset) + (size + offset);

				pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

				int paired_size = paired_size_and_offset.first;
				int paired_offset = paired_size_and_offset.second;

				uint64_t paired_counter_value = read_counter(paired_counter_index, row);

				uint64_t max_value = paired_counter_value > (counter_value + 1) ? paired_counter_value : counter_value + 1;

				int bitSize = counterSize * (size + paired_size);
				int bitOffset = offset * counterSize;
				int negativeOffset = 64 - bitSize - bitOffset;
				int wordOffset = index >> log_counters_per_word;

				// write to counter
				uint64_t& cnt = baseline_cms_counters[row][wordOffset];

				uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
				uint64_t cnt_masked = cnt & cnt_mask;
				uint64_t cnt_shifted = max_value << negativeOffset;

				cnt = cnt_shifted | cnt_masked;

				// turn on bit
				uint64_t byte_index = index >> 6;

				uint8_t bitIndex = (paired_counter_index - 1) & 0b111111;

				uint64_t& merges_word = baseline_cms_merges[row][byte_index];

				merges_word = merges_word | ((uint64_t)1 << bitIndex);

			}

			else // merge left
			{
				int paired_counter_index = (index - in_word_offset) + (offset - 1);

				pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

				int paired_size = paired_size_and_offset.first;
				int paired_offset = paired_size_and_offset.second;

				uint64_t paired_counter_value = read_counter(paired_counter_index, row);

				uint64_t max_value = paired_counter_value > (counter_value + 1) ? paired_counter_value : counter_value + 1;

				int bitSize = counterSize * (size + paired_size);
				int bitOffset = paired_offset * counterSize;
				int negativeOffset = 64 - bitSize - bitOffset;
				int wordOffset = index >> log_counters_per_word;

				// write to counter
				uint64_t& cnt = baseline_cms_counters[row][wordOffset];

				uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
				uint64_t cnt_masked = cnt & cnt_mask;
				uint64_t cnt_shifted = max_value << negativeOffset;

				cnt = cnt_shifted | cnt_masked;

				// turn on bit
				uint64_t byte_index = index >> 6;

				uint8_t bitIndex = paired_counter_index & 0b111111;

				uint64_t& merges_word = baseline_cms_merges[row][byte_index];

				merges_word = merges_word | ((uint64_t)1 << bitIndex);

			}

		}


	}
}

inline uint64_t smarTangoBaselineSanity::read_counter(uint index, uint row)
{

	pair<uint8_t, uint8_t> size_and_offset = counter_size_and_offset(index, row);

	int size = size_and_offset.first;
	int offset = size_and_offset.second;

	int bitSize = size * counterSize;
	int bitOffset = offset * counterSize;
	int negativeOffset = 64 - bitSize - bitOffset;
	int wordOffset = index >> log_counters_per_word;

	return  (baseline_cms_counters[row][wordOffset] << bitOffset) >> (negativeOffset + bitOffset);

}

uint64_t smarTangoBaselineSanity::query(const char * str)
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

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

TangoBaselineSanity::TangoBaselineSanity()
{
}

TangoBaselineSanity::~TangoBaselineSanity()
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

void TangoBaselineSanity::initialize(int counterSize, int width, int height, int seed)
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
}

inline pair<uint8_t, uint8_t> TangoBaselineSanity::counter_size_and_offset(uint index, uint row)
{

	// 64 counters in each index  
	uint64_t byte_index = index >> 6;

	// in-word index
	uint8_t index_mod_64 = index & 0b111111;

	uint64_t merges_word = baseline_cms_merges[row][byte_index];

	int size = 1;
	uint8_t bitIndex = index_mod_64;

	while ((bitIndex < 64) && (isBitSet(merges_word, bitIndex++)))
	{
		++size;
	}

	uint8_t offset = index_mod_64 - 1;

	while ((offset >= 0) && (isBitSet(merges_word, offset)))
	{
		++size;
		--offset;
	}

	return make_pair(size, (offset + 1) & (counters_per_word - 1));

}

void TangoBaselineSanity::increment(const char * str)
{
	for (int row = 0; row < height; ++row) {
		uint index = ((bobhash[row].run(str, FT_SIZE)) & width_mask) >> (6 - log_counters_per_word);

		pair<uint8_t, uint8_t> size_and_offset = counter_size_and_offset(index, row);

		int size = size_and_offset.first;
		int offset = size_and_offset.second;

		uint64_t counter_value = read_counter(index, row);

		if (counter_value < (((uint64_t)1 << (size*counterSize))) - 1)
		{
			// no overflow
			int bitSize = size * counterSize;
			int bitOffset = offset * counterSize;
			int negativeOffset = 64 - bitSize - bitOffset;
			int wordOffset = index >> log_counters_per_word;

			uint64_t& cnt = baseline_cms_counters[row][wordOffset];

			uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
			uint64_t cnt_masked = cnt & cnt_mask;
			uint64_t cnt_shifted = (counter_value + 1) << negativeOffset;

			cnt = cnt_shifted | cnt_masked;

		}
		else
		{

			int log_size = 0;
			int i = size;
			while (i >>= 1) ++log_size;

			int in_word_offset = index & ((1 << log_counters_per_word) - 1);

			bool merge_right;
			bool found_available_counter = false;

			int merge_bounries[6];

			merge_bounries[0] = (index & 0b1) == 0 ? index | 0b1 : index & 0xFFFFFFFE;
			merge_bounries[1] = (index & 0b10) == 0 ? index | 0b11 : index & 0xFFFFFFFC;
			merge_bounries[2] = (index & 0b100) == 0 ? index | 0b111 : index & 0xFFFFFFF8;
			merge_bounries[3] = (index & 0b1000) == 0 ? index | 0b1111 : index & 0xFFFFFFF0;
			merge_bounries[4] = (index & 0b10000) == 0 ? index | 0b11111 : index & 0xFFFFFFE0;
			merge_bounries[5] = (index & 0b100000) == 0 ? index | 0b111111 : index & 0xFFFFFFC0;

			uint64_t& curr_merges_word = baseline_cms_merges[row][index >> 6];

			for (int j = 0; (j < 6) && (!(found_available_counter)); ++j)
			{
				if (index > merge_bounries[j])
				{
					for (int current_index = index - 1; current_index >= merge_bounries[j]; --current_index)
					{
						uint8_t bitIndex = current_index & 0b111111;
						if (!(isBitSet(curr_merges_word, bitIndex)))
						{
							found_available_counter = true;
							merge_right = false;
							break;
						}
					}
				}
				else
				{
					for (int current_index = index; current_index < merge_bounries[j]; ++current_index)
					{
						uint8_t bitIndex = current_index & 0b111111;
						if (!(isBitSet(curr_merges_word, bitIndex)))
						{
							found_available_counter = true;
							merge_right = true;
							break;
						}
					}
				}
			}

			if (merge_right)
			{
				int paired_counter_index = (index - in_word_offset) + (size + offset);

				pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

				int paired_size = paired_size_and_offset.first;
				int paired_offset = paired_size_and_offset.second;

				uint64_t paired_counter_value = read_counter(paired_counter_index, row);

				uint64_t max_value = paired_counter_value > (counter_value + 1) ? paired_counter_value : counter_value + 1;

				int bitSize = counterSize * (size + paired_size);
				int bitOffset = offset * counterSize;
				int negativeOffset = 64 - bitSize - bitOffset;
				int wordOffset = index >> log_counters_per_word;

				// write to counter
				uint64_t& cnt = baseline_cms_counters[row][wordOffset];

				uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
				uint64_t cnt_masked = cnt & cnt_mask;
				uint64_t cnt_shifted = max_value << negativeOffset;

				cnt = cnt_shifted | cnt_masked;

				// turn on bit
				uint64_t byte_index = index >> 6;

				uint8_t bitIndex = (paired_counter_index - 1) & 0b111111;

				uint64_t& merges_word = baseline_cms_merges[row][byte_index];

				merges_word = merges_word | ((uint64_t)1 << bitIndex);

			}
			else // merge left
			{
				int paired_counter_index = (index - in_word_offset) + (offset - 1);

				pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

				int paired_size = paired_size_and_offset.first;
				int paired_offset = paired_size_and_offset.second;

				uint64_t paired_counter_value = read_counter(paired_counter_index, row);

				uint64_t max_value = paired_counter_value > (counter_value + 1) ? paired_counter_value : counter_value + 1;

				int bitSize = counterSize * (size + paired_size);
				int bitOffset = paired_offset * counterSize;
				int negativeOffset = 64 - bitSize - bitOffset;
				int wordOffset = index >> log_counters_per_word;

				// write to counter
				uint64_t& cnt = baseline_cms_counters[row][wordOffset];

				uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
				uint64_t cnt_masked = cnt & cnt_mask;
				uint64_t cnt_shifted = max_value << negativeOffset;

				cnt = cnt_shifted | cnt_masked;

				// turn on bit
				uint64_t byte_index = index >> 6;

				uint8_t bitIndex = paired_counter_index & 0b111111;

				uint64_t& merges_word = baseline_cms_merges[row][byte_index];

				merges_word = merges_word | ((uint64_t)1 << bitIndex);

			}

		}

	}
}

inline uint64_t TangoBaselineSanity::read_counter(uint index, uint row)
{

	pair<uint8_t, uint8_t> size_and_offset = counter_size_and_offset(index, row);

	int size = size_and_offset.first;
	int offset = size_and_offset.second;

	int bitSize = size * counterSize;
	int bitOffset = offset * counterSize;
	int negativeOffset = 64 - bitSize - bitOffset;
	int wordOffset = index >> log_counters_per_word;

	return  (baseline_cms_counters[row][wordOffset] << bitOffset) >> (negativeOffset + bitOffset);

}

uint64_t TangoBaselineSanity::query(const char * str)
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
