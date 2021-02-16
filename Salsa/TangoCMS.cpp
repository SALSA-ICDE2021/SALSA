#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "TangoCMS.hpp"

using namespace std;

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

smarTango::smarTango(int counterSize, int width, int height, int seed, int downsamplings_per_merge) :
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

	////
	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	counter_indices = new int[height];
	counter_sizes = new uint8_t[height];
	counter_offsets = new uint8_t[height];
	counter_values = new uint64_t[height];

	largest_counter = 1;
	this->downsamplings_per_merge = downsamplings_per_merge;
	remaining_downsamplings_per_merge = downsamplings_per_merge;
	////
}

smarTango::~smarTango()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms_counters[i];
		delete[] baseline_cms_merges[i];
	}

	delete[] bobhash;
	delete[] baseline_cms_counters;

	delete[] baseline_cms_merges;

	delete[] counter_sizes;
	delete[] counter_offsets;

}

inline pair<uint8_t, uint8_t> smarTango::counter_size_and_offset(uint index, uint row)
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

void smarTango::increment_h(const char * str)
{

	// remember all counter indices and the sizes
	bool is_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		counter_indices[row] = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		counter_values[row] = read_counter(counter_indices[row], row);

		pair<uint8_t, uint8_t> current_counter_size_and_offset = counter_size_and_offset(counter_indices[row], row);

		counter_sizes[row] = current_counter_size_and_offset.first;
		counter_offsets[row] = current_counter_size_and_offset.second;

		int bitSize = counterSize * counter_sizes[row];

		if (counter_values[row] == (((uint64_t)1 << bitSize) - 1))
		{
			is_overflow = true;
		}
	}

	// safe to increment counters
	if (!is_overflow)
	{
		for (int row = 0; row < height; ++row) 
		{
			uint index = counter_indices[row];
			int size = counter_sizes[row];
			int offset = counter_offsets[row];
			uint64_t counter_value = counter_values[row];

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
	}
	// NOT safe to increment counters
	else
	{
		handle_overflow_and_increment();
	}

}

void smarTango::increment(const char * str)
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

inline uint64_t smarTango::read_counter(uint index, uint row)
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

uint64_t smarTango::query(const char * str)
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
	return ((uint64_t)min) << minus_log_p;
}

////
void smarTango::handle_overflow_and_increment()
{

	bool is_max_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		int bitSize = counterSize * counter_sizes[row];
		if ((counter_values[row] == (((uint64_t)1 << bitSize) - 1)) && (what_would_be_the_merged_counter_size(row) > largest_counter))
		{
			is_max_overflow = true;
		}
	}

	if (!is_max_overflow) // merge
	{
		increment_and_potentially_merge();
	}
	else if (remaining_downsamplings_per_merge > 0) // do downsampling
	{
		--remaining_downsamplings_per_merge;

		++minus_log_p;
		LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
		divide_counters();

		rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
	}
	else
	{
		// update of largest_counter and remaining_downsamplings_per_merge moved to "increment_and_potentially_merge"
		increment_and_potentially_merge();
	}
	
	
}

void smarTango::divide_counters()
{	
	cout << "divide_counters" << endl;
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			pair<uint8_t, uint8_t> current_counter_size_and_offset = counter_size_and_offset(current_counter_index, row);

			int size = current_counter_size_and_offset.first;
			int offset = current_counter_size_and_offset.second;
			uint64_t counter_value = read_counter(current_counter_index, row);

			if (current_counter_index == counter_indices[row])
			{
				++counter_value;
			}

			int bitSize = size * counterSize;
			int bitOffset = offset * counterSize;
			int negativeOffset = 64 - bitSize - bitOffset;

			int wordOffset = current_counter_index >> log_counters_per_word;

			uint64_t& cnt = baseline_cms_counters[row][wordOffset];

			uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
			uint64_t cnt_masked = cnt & cnt_mask;

			uint64_t new_val = counter_value >> 1;
			if (counter_value % 2)
			{
				new_val += rand() % 2;
			}
			uint64_t cnt_shifted = new_val << negativeOffset;

			cnt = cnt_shifted | cnt_masked;

			current_counter_index += size;
		}
	}	
}

void smarTango::increment_and_potentially_merge()
{
	for (int row = 0; row < height; ++row) {

		uint index = counter_indices[row];
		int size = counter_sizes[row];
		int offset = counter_offsets[row];
		uint64_t counter_value = counter_values[row];

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
			bool found_available_counter = false;

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

				////
				if (largest_counter < size + paired_size)
				{
					remaining_downsamplings_per_merge += downsamplings_per_merge*(size + paired_size - largest_counter);
					largest_counter = size + paired_size;
				}
				
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

				////
				if (largest_counter < size + paired_size)
				{
					remaining_downsamplings_per_merge += downsamplings_per_merge * (size + paired_size - largest_counter);
					largest_counter = size + paired_size;
				}


			}

		}


	}
}

uint8_t smarTango::what_would_be_the_merged_counter_size(uint row)
{

	bool merge_right;
	bool found_available_counter = false;

	uint index = counter_indices[row];
	int size = counter_sizes[row];
	int offset = counter_offsets[row];
	uint64_t counter_value = counter_values[row];

	int in_word_offset = index & ((1 << log_counters_per_word) - 1);

	int bitSize = size * counterSize;
	int bitOffset = offset * counterSize;
	int negativeOffset = 64 - bitSize - bitOffset;

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

	int paired_size;

	if (merge_right)
	{
		int paired_counter_index = (index - in_word_offset) + (size + offset);

		pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

		paired_size = paired_size_and_offset.first;
	}
	else
	{
		int paired_counter_index = (index - in_word_offset) + (offset - 1);

		pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

		paired_size = paired_size_and_offset.first;
	}
	return size + paired_size;
}
////

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////

Tango::Tango(int counterSize, int width, int height, int seed, int downsamplings_per_merge) :
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

	////
	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	counter_indices = new int[height];
	counter_sizes = new uint8_t[height];
	counter_offsets = new uint8_t[height];
	counter_values = new uint64_t[height];

	largest_counter = 1;
	this->downsamplings_per_merge = downsamplings_per_merge;
	remaining_downsamplings_per_merge = downsamplings_per_merge;
	////
}

Tango::~Tango()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms_counters[i];
		delete[] baseline_cms_merges[i];
	}

	delete[] bobhash;
	delete[] baseline_cms_counters;

	delete[] baseline_cms_merges;

	delete[] counter_sizes;
	delete[] counter_offsets;

}

inline pair<uint8_t, uint8_t> Tango::counter_size_and_offset(uint index, uint row)
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

void Tango::increment_h(const char * str)
{

	// remember all counter indices and the sizes
	bool is_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		counter_indices[row] = (bobhash[row].run(str, FT_SIZE)) & width_mask;
		counter_values[row] = read_counter(counter_indices[row], row);

		pair<uint8_t, uint8_t> current_counter_size_and_offset = counter_size_and_offset(counter_indices[row], row);

		counter_sizes[row] = current_counter_size_and_offset.first;
		counter_offsets[row] = current_counter_size_and_offset.second;

		int bitSize = counterSize * counter_sizes[row];

		if (counter_values[row] == (((uint64_t)1 << bitSize) - 1))
		{
			is_overflow = true;
		}
	}

	// safe to increment counters
	if (!is_overflow)
	{
		for (int row = 0; row < height; ++row)
		{
			uint index = counter_indices[row];
			int size = counter_sizes[row];
			int offset = counter_offsets[row];
			uint64_t counter_value = counter_values[row];

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
	}
	// NOT safe to increment counters
	else
	{
		handle_overflow_and_increment();
	}

}

void Tango::increment(const char * str)
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

inline uint64_t Tango::read_counter(uint index, uint row)
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

uint64_t Tango::query(const char * str)
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
	return ((uint64_t)min) << minus_log_p;
}

////
void Tango::handle_overflow_and_increment()
{

	bool is_max_overflow = false;
	for (int row = 0; row < height; ++row)
	{
		int bitSize = counterSize * counter_sizes[row];
		if ((counter_values[row] == (((uint64_t)1 << bitSize) - 1)) && (what_would_be_the_merged_counter_size(row) > largest_counter))
		{
			is_max_overflow = true;
		}
	}

	if (!is_max_overflow) // merge
	{
		increment_and_potentially_merge();
	}
	else if (remaining_downsamplings_per_merge > 0) // do downsampling
	{
		--remaining_downsamplings_per_merge;

		++minus_log_p;
		LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
		divide_counters();

		rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
	}
	else
	{
		// update of largest_counter and remaining_downsamplings_per_merge moved to "increment_and_potentially_merge"
		increment_and_potentially_merge();
	}


}

void Tango::divide_counters()
{
	cout << "divide_counters" << endl;
	for (int row = 0; row < height; ++row)
	{
		int current_counter_index = 0;
		while (current_counter_index < width)
		{
			pair<uint8_t, uint8_t> current_counter_size_and_offset = counter_size_and_offset(current_counter_index, row);

			int size = current_counter_size_and_offset.first;
			int offset = current_counter_size_and_offset.second;
			uint64_t counter_value = read_counter(current_counter_index, row);

			if (current_counter_index == counter_indices[row])
			{
				++counter_value;
			}

			int bitSize = size * counterSize;
			int bitOffset = offset * counterSize;
			int negativeOffset = 64 - bitSize - bitOffset;

			int wordOffset = current_counter_index >> log_counters_per_word;

			uint64_t& cnt = baseline_cms_counters[row][wordOffset];

			uint64_t cnt_mask = ~((((uint64_t)1 << bitSize) - 1) << negativeOffset);
			uint64_t cnt_masked = cnt & cnt_mask;

			uint64_t new_val = counter_value >> 1;
			if (counter_value % 2)
			{
				new_val += rand() % 2;
			}
			uint64_t cnt_shifted = new_val << negativeOffset;

			cnt = cnt_shifted | cnt_masked;

			current_counter_index += size;
		}
	}
}

void Tango::increment_and_potentially_merge()
{
	for (int row = 0; row < height; ++row) {

		uint index = counter_indices[row];
		int size = counter_sizes[row];
		int offset = counter_offsets[row];
		uint64_t counter_value = counter_values[row];

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

			////
			int log_size = 0;
			int i = size;
			while (i >>= 1) ++log_size;

			int in_word_offset = index & ((1 << log_counters_per_word) - 1);

			////

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
			////

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

				////
				if (largest_counter < size + paired_size)
				{
					remaining_downsamplings_per_merge += downsamplings_per_merge * (size + paired_size - largest_counter);
					largest_counter = size + paired_size;
				}

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

				////
				if (largest_counter < size + paired_size)
				{
					remaining_downsamplings_per_merge += downsamplings_per_merge * (size + paired_size - largest_counter);
					largest_counter = size + paired_size;
				}


			}

		}


	}
}

uint8_t Tango::what_would_be_the_merged_counter_size(uint row)
{

	bool merge_right;
	bool found_available_counter = false;

	uint index = counter_indices[row];
	int size = counter_sizes[row];
	int offset = counter_offsets[row];
	uint64_t counter_value = counter_values[row];

	int in_word_offset = index & ((1 << log_counters_per_word) - 1);

	int bitSize = size * counterSize;
	int bitOffset = offset * counterSize;
	int negativeOffset = 64 - bitSize - bitOffset;

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

	int paired_size;

	if (merge_right)
	{
		int paired_counter_index = (index - in_word_offset) + (size + offset);

		pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

		paired_size = paired_size_and_offset.first;
	}
	else
	{
		int paired_counter_index = (index - in_word_offset) + (offset - 1);

		pair<uint8_t, uint8_t> paired_size_and_offset = counter_size_and_offset(paired_counter_index, row);

		paired_size = paired_size_and_offset.first;
	}
	return size + paired_size;
}
////

////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////