#include "smartSALSA.hpp"

#define uint unsigned int
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

SmartEncodingSALSA::SmartEncodingSALSA()
{
}

SmartEncodingSALSA::~SmartEncodingSALSA()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] basic_counters[i];
		delete[] twelve_counter_encodings[i];
	}

#ifdef USE_BOBHASH
	delete[] bobhash;
#endif // USE_BOBHASH

	delete[] basic_counters;

	delete[] size_16_counter_alias;
	delete[] size_32_counter_alias;

	delete[] twelve_counter_encodings;
}

void SmartEncodingSALSA::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	int index = width_mask;
	log_width = 0;
	while (index >>= 1) ++log_width;

	assert(width > 0 && "We assume too much!");
	assert(width % 16 == 0 && "We assume that (w % 16 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	for (int i = 0; i < 512; ++i) {
		mod12_lookup[i] = i % 12;
	}

#ifdef USE_ONE_XXHASH
	if (height > 1) {
		assert(width <= ((uint64_t)1 << (64 / height)) && "One XXHASH is not enough!");
	}
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	height_plus_1_over_2 = (height + 1) / 2;
	assert(width < ((uint64_t)1 << (64 / height_plus_1_over_2)) && "Two XXHASHes are not enough!");
#endif // USE_TWO_XXHASHES

	basic_counters = new uint8_t * [height];

	size_16_counter_alias = new uint16_t * [height];
	size_32_counter_alias = new uint32_t * [height];

	twelve_counter_encodings = new uint8_t * [height];

#ifdef USE_BOBHASH
	bobhash = new BOBHash[height];
#endif // USE_BOBHASH

	for (int row = 0; row < height; ++row)
	{
		basic_counters[row] = new uint8_t[width]();

		size_16_counter_alias[row] = (uint16_t*)basic_counters[row];
		size_32_counter_alias[row] = (uint32_t*)basic_counters[row];

		twelve_counter_encodings[row] = new uint8_t[(width / 12) + 1]();

#ifdef USE_BOBHASH
		bobhash[row].initialize(seed * (7 + row) + row + 100);
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
		seeds[row] = seed * (7 + row) + row + 100;
#endif // USE_XXHASH
	}

	int logSizes[5][4];
	for (int four_encoding = 0; four_encoding < 5; ++four_encoding) {
		for (int counter_number = 0; counter_number < 4; ++counter_number) {
			if (four_encoding == 4) { // all are merged 
				logSizes[four_encoding][counter_number] = 2;
			}
			else if (four_encoding == 3) { // merged in pairs
				logSizes[four_encoding][counter_number] = 1;
			}
			else if (four_encoding == 2) { // left pair is merged
				if (counter_number < 2) {
					logSizes[four_encoding][counter_number] = 1;
				}
				else {
					logSizes[four_encoding][counter_number] = 0;
				}
			}
			else if (four_encoding == 1) {
				if (counter_number < 2) {
					logSizes[four_encoding][counter_number] = 0;
				}
				else {
					logSizes[four_encoding][counter_number] = 1;
				}
			}
			else if (four_encoding == 0) { // none are merged
				logSizes[four_encoding][counter_number] = 0;
			}

		}
	}

	//for (int encoding = 0; encoding < 125 * 12; ++encoding)
	for (int encoding = 0; encoding < 256 * 12; ++encoding)
	{
		//int counter_number = encoding % 12;
		//int residual_encoding = encoding / 12;
		int counter_number = encoding >> 8;
		int residual_encoding = encoding & 0xFF;

		int first_four_encoding = residual_encoding / 25;
		int second_four_encoding = (residual_encoding / 5) % 5;
		int third_four_encoding = residual_encoding % 5;

		if (counter_number < 4) { //first four
			smart_encoding_lookup_table[encoding] = logSizes[first_four_encoding][counter_number % 4];
		}
		else if (counter_number < 8) {//second four
			smart_encoding_lookup_table[encoding] = logSizes[second_four_encoding][counter_number % 4];
		}
		else {//third four
			smart_encoding_lookup_table[encoding] = logSizes[third_four_encoding][counter_number % 4];
		}
		int merges_byte = residual_encoding;
		int smarter_0_encoding = residual_encoding;
		int smarter_1_encoding = residual_encoding;
		if (smart_encoding_lookup_table[encoding] == 0) { // for merging 8b counters
			if (counter_number < 4) { //counters 0-3
				int first_four_encoding = merges_byte / 25;
				if (first_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 2) { //counters 0-1
						merges_byte += 25 * 2;
						//first_four_encoding = 2;
					}
					else {//counters 2-3
						merges_byte += 25;
						//first_four_encoding = 1;
					}
					//merges_byte = (first_four_encoding * 25) + (merges_byte % 25);
				}
				else { // The other counters were merged
					merges_byte = (25 * 3) + (merges_byte % 25);
				}
			}
			else if (counter_number < 8) { //counters 4-7
				int second_four_encoding = (merges_byte / 5) % 5;
				if (second_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 6) { //counters 4-5
						merges_byte += 5 * 2;
						//second_four_encoding = 2;
					}
					else {//counters 6-7
						merges_byte += 5;
						//second_four_encoding = 1;
					}
				}
				else { // The other counters were merged
					int first_four_encoding = merges_byte / 25;
					int third_four_encoding = merges_byte % 5;
					merges_byte = (first_four_encoding * 25) + (3 * 5) + third_four_encoding;
				}
			}
			else { //counters 8-11
				int third_four_encoding = merges_byte % 5;
				if (third_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 10) { //counters 8-9
						merges_byte += 2;
						//third_four_encoding = 2;
					}
					else {//counters 10-11
						merges_byte += 1;
						//third_four_encoding = 1;
					}
				}
				else { // The other counters were merged
					merges_byte = merges_byte - third_four_encoding + 3;
				}
			}
			smarter_0_encoding = merges_byte;
		}
		else { // for merging 16b counters
			if (counter_number < 4) { //counters 0-3
		//int first_four_encoding = merges_byte / 25;
		//merges_byte += (4 - first_four_encoding) * 25;
				merges_byte = 4 * 25 + (merges_byte % 25);
			}
			else if (counter_number < 8) { //counters 4-7
				int first_four_encoding = merges_byte / 25;
				int third_four_encoding = merges_byte % 5;
				merges_byte = (first_four_encoding * 25) + (4 * 5) + third_four_encoding;
			}
			else { //counters 8-11
				merges_byte = (merges_byte / 5) * 5 + 4;
			}
			smarter_1_encoding = merges_byte;
		}
		//smart_merging_lookup_table[encoding] = merges_byte;
		smarter_merging_lookup_table_0[encoding] = smarter_0_encoding;
		smarter_merging_lookup_table_1[encoding] = smarter_1_encoding;
	}

}

inline uint8_t SmartEncodingSALSA::log_counter_size(uint index, uint row)
{
	uint index_over_12 = index / 12;
	uint index_mod_12 = index - 12 * index_over_12; //index % 12;

	uint16_t merges_byte = (uint16_t)twelve_counter_encodings[row][index_over_12];


	return smart_encoding_lookup_table[(index_mod_12 << 8) + merges_byte];
	//return smart_encoding_lookup_table[12 * merges_byte + index_mod_12];
	

	// we assume that counter sizes cannot exceed 4X the original
	// if the original is 8 bits then we have a 32-bit limit

}

#if false
void SmartEncodingSALSA::merge_8_bit_counters(uint index, uint row)
{
	//uint index_mod_12 = index - 12 * index_over_12; //index % 12;
	uint index_mod_12 = index % 12;

	uint index_over_12 = index / 12;
	uint8_t& merges_byte = twelve_counter_encodings[row][index_over_12];
	//merges_byte = smarter_merging_lookup_table_0[12 * merges_byte + index_mod_12];
	merges_byte = smarter_merging_lookup_table_0[(index_mod_12 << 8) + merges_byte];
	
	/*
	return;

	uint8_t log_size = smart_encoding_lookup_table[12 * merges_byte + index_mod_12];
	if (log_size > 0)
	{
		// counters already merged
	}
	else
	{
		merges_byte = smart_merging_lookup_table[merges_byte*12 + index_mod_12];
		return;
		if (index_mod_12 < 4) { //counters 0-3
			int first_four_encoding = merges_byte / 25;
			if (first_four_encoding == 0) { // all counters were 8-bits
				if (index_mod_12 < 2) { //counters 0-1
					merges_byte += 25 * 2;
					//first_four_encoding = 2;
				}
				else {//counters 2-3
					merges_byte += 25;
					//first_four_encoding = 1;
				}
				//merges_byte = (first_four_encoding * 25) + (merges_byte % 25);
			}
			else { // The other counters were merged
				merges_byte = (25 * 3) + (merges_byte % 25);
			}
		}
		else if (index_mod_12 < 8) { //counters 4-7
			int second_four_encoding = (merges_byte / 5) % 5;
			if (second_four_encoding == 0) { // all counters were 8-bits
				if (index_mod_12 < 6) { //counters 4-5
					merges_byte += 5 * 2;
					//second_four_encoding = 2;
				}
				else {//counters 6-7
					merges_byte += 5;
					//second_four_encoding = 1;
				}
			}
			else { // The other counters were merged
				int first_four_encoding = merges_byte / 25;
				int third_four_encoding = merges_byte % 5;
				merges_byte = (first_four_encoding * 25) + (3 * 5) + third_four_encoding;
			}
		}
		else { //counters 8-11
			int third_four_encoding = merges_byte % 5;
			if (third_four_encoding == 0) { // all counters were 8-bits
				if (index_mod_12 < 10) { //counters 8-9
					merges_byte += 2;
					//third_four_encoding = 2;
				}
				else {//counters 10-11
					merges_byte += 1;
					//third_four_encoding = 1;
				}
			}
			else { // The other counters were merged
				merges_byte = merges_byte - third_four_encoding + 3;
			}
		}

		// happens in increment
		//uint16_t* p = (uint16_t*)&baseline_cms_counters[row][index & 0xFFFFFFFE];
		//*p = max(baseline_cms_counters[row][index], baseline_cms_counters[row][index ^ 0b1]);
	}*/
}

//void SmartEncodingSALSA::merge_16_bit_counters(uint index, uint row)
{
	uint index_over_12 = index / 12;
	uint index_mod_12 = index - 12 * index_over_12; //index % 12;

	uint8_t& merges_byte = twelve_counter_encodings[row][index_over_12];
	//merges_byte = smarter_merging_lookup_table_1[12 * merges_byte + index_mod_12];
	merges_byte = smarter_merging_lookup_table_1[(index_mod_12 << 8) + merges_byte];
	return;
	/*
	merges_byte = smart_merging_lookup_table[merges_byte * 12 + index_mod_12];
	return;

	//uint8_t log_size = smart_merges_lookup_table[12 * ((uint16_t)merges_byte) + index_mod_12];
	if (index_mod_12 < 4) { //counters 0-3
		//int first_four_encoding = merges_byte / 25;
		//merges_byte += (4 - first_four_encoding) * 25;
		merges_byte = 4 * 25 + (merges_byte % 25);
	}
	else if (index_mod_12 < 8) { //counters 4-7
		int first_four_encoding = merges_byte / 25;
		int third_four_encoding = merges_byte % 5;
		merges_byte = (first_four_encoding * 25) + (4 * 5) + third_four_encoding;
	}
	else { //counters 8-11
		merges_byte = (merges_byte/5)*5 + 4;
	}
	*/
}
#endif
void SmartEncodingSALSA::increment(const char* str)
{

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
#endif // USE_TWO_XXHASHES

	for (int row = 0; row < height; ++row) {

#ifdef USE_BOBHASH
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[row])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
		uint index;
		if (row <= height_plus_1_over_2)
		{
			index = hashes0 & width_mask;
			hashes0 >>= log_width;
		}
		else
		{
			index = hashes1 & width_mask;
			hashes1 >>= log_width;
		}
#endif // USE_TWO_XXHASHES


		//uint8_t log_counter_size_value = log_counter_size(index, row);
		uint index_over_12 = index / 12;
		//uint index_mod_12 = index - 12 * index_over_12; //index % 12;
		//uint index_mod_12 = mod12_lookup[((mod12_lookup[index >> 16] + mod12_lookup[(index >> 8) & 0xFF]) << 2) + mod12_lookup[index & 0xFF]];
		uint index_mod_12 = mod12_lookup[(mod12_lookup[(index >> 16) + ((index >> 8) & 0xFF)] << 2) + mod12_lookup[index & 0xFF]];

		uint8_t &merges_byte = twelve_counter_encodings[row][index_over_12];
		uint16_t lookup_key = (index_mod_12 << 8) | merges_byte;


		//uint8_t log_counter_size_value = smart_encoding_lookup_table[12 * merges_byte + index_mod_12];
		uint8_t log_counter_size_value = smart_encoding_lookup_table[lookup_key];

		

		/*
		uint inc_index = index & ~((1 << log_counter_size_value) - 1);
		
		uint32_t& cnt = *((uint32_t*)(basic_counters[row] + inc_index));
		//++cnt;
		if (++cnt & 0xFF) {
			//continue;
		}	
		else {
			if (log_counter_size_value == 0) {
				if ((inc_index & 0b1) || (((inc_index & 0b1) == 0) && ((cnt & 0xFFFF) == 0))) {
					--cnt;
				}
				uint index_16 = index >> 1;
				size_16_counter_alias[row][index_16] = 1 << 8;
				merges_byte = smarter_merging_lookup_table_0[lookup_key];
			}
			else if (log_counter_size_value == 1) {
				if ((cnt & 0xFFFF) == 0) {
					if (inc_index & 0b10) {
						--cnt;
					}

					uint index_32 = index >> 2;
					size_32_counter_alias[row][index_32] = 1 << 16;
					merges_byte = smarter_merging_lookup_table_1[lookup_key];
				}
			}
			continue;
		}*/
		uint inc_index = (log_counter_size_value == 0) ? index : (log_counter_size_value == 1) ? index & 0xFFFFFFFE : index & 0xFFFFFFFC;
		//uint8_t* &basic_counters_row = basic_counters[row];
		//if (++basic_counters_row[inc_index]) {
		if (++basic_counters[row][inc_index]) {
			//continue;
		}
		else {
			//uint8_t& merges_byte = eight_counter_encodings[row][index_over_8];
			if (log_counter_size_value == 1) {
				if (++basic_counters[row][inc_index + 1]) {
					//if (++basic_counters_row[inc_index + 1] & 0xFF) {
				}
				else {
					uint index_32 = index >> 2;
					size_32_counter_alias[row][index_32] = 1 << 16;
					merges_byte = smarter_merging_lookup_table_1[lookup_key];
				}

			}
			else if (log_counter_size_value == 0) {
				//uint index_16 = index >> 1;
				//size_16_counter_alias[row][index_16] = 1 << 8;
				basic_counters[row][inc_index | 1] = 1;
				basic_counters[row][inc_index & 0xFFFFFFFE] = 0;
				merges_byte = smarter_merging_lookup_table_0[lookup_key];
			}
			else {
				if (++basic_counters[row][inc_index + 1]) {
				}
				else {
					++basic_counters[row][inc_index + 2];
				}
			}

			//continue;
		}
#if false
		// we are here... 
		if (log_counter_size_value == 0)
		{
			if (++basic_counters[row][index])
			{
				// no overflow
			}
			else
			{
				//merge_8_bit_counters(index, row);
				merges_byte = smarter_merging_lookup_table_0[lookup_key];

				// the current counter became 0 while its real value is 256...
				uint index_16 = index >> 1;
				size_16_counter_alias[row][index_16] = 1 << 8;
			}
		}
		else if (log_counter_size_value == 1)
		{
			uint index_16 = index >> 1;
			if (++size_16_counter_alias[row][index_16])
			{
				// no overflow
			}
			else
			{
				//merge_16_bit_counters(index, row);
				merges_byte = smarter_merging_lookup_table_1[lookup_key];

				// the current counter became 0 while its real value is 2^16...
				uint index_32 = index >> 2;
				size_32_counter_alias[row][index_32] = 1 << 16;
			}
		}
		else if (log_counter_size_value == 2)
		{
			uint index_32 = index >> 2;
			if (++size_32_counter_alias[row][index_32])
			{
				// no overflow
			}
			else
			{
				assert(false);
			}
		}
		else // if (log_counter_size_value == 3)
		{
			assert(false);
		}
#endif
	}
}

inline uint64_t SmartEncodingSALSA::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);

	switch (log_counter_size_value) {
	case 0:
		return basic_counters[row][index];
	case 1:
		return size_16_counter_alias[row][index >> 1];
	case 2:
		return size_32_counter_alias[row][index >> 2];
	case 3:
		// not implemented
		return -1;
	default:
		// BOB!
		return -1;
	}
}

uint64_t SmartEncodingSALSA::query(const char* str)
{

#ifdef USE_BOBHASH
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
	uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[0])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint index = hashes & width_mask;
	hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
	uint index = hashes0 & width_mask;
	hashes0 >>= log_width;
#endif // USE_TWO_XXHASHES

	int row = 0;

	uint64_t min = read_counter(index, row++);

	while (row < height) {

#ifdef USE_BOBHASH
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[row])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
		uint index;
		if (row <= height_plus_1_over_2)
		{
			index = hashes0 & width_mask;
			hashes0 >>= log_width;
		}
		else
		{
			index = hashes1 & width_mask;
			hashes1 >>= log_width;
		}
#endif 
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

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

FastEncodingSALSA::FastEncodingSALSA()
{
}

FastEncodingSALSA::~FastEncodingSALSA()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] basic_counters[i];
		delete[] eight_counter_encodings[i];
	}

#ifdef USE_BOBHASH
	delete[] bobhash;
#endif // USE_BOBHASH

	delete[] basic_counters;

	delete[] size_16_counter_alias;
	delete[] size_32_counter_alias;

	delete[] eight_counter_encodings;
}

void FastEncodingSALSA::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	int index = width_mask;
	log_width = 0;
	while (index >>= 1) ++log_width;

	assert(width > 0 && "We assume too much!");
	assert(width % 16 == 0 && "We assume that (w % 16 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

#ifdef USE_ONE_XXHASH
	if (height > 1) {
		assert(width <= ((uint64_t)1 << (64 / height)) && "One XXHASH is not enough!");
	}
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	height_plus_1_over_2 = (height + 1) / 2;
	assert(width < ((uint64_t)1 << (64 / height_plus_1_over_2)) && "Two XXHASHes are not enough!");
#endif // USE_TWO_XXHASHES

	basic_counters = new uint8_t * [height];

	size_16_counter_alias = new uint16_t * [height];
	size_32_counter_alias = new uint32_t * [height];

	eight_counter_encodings = new uint8_t * [height];

#ifdef USE_BOBHASH
	bobhash = new BOBHash[height];
#endif // USE_BOBHASH

	for (int row = 0; row < height; ++row)
	{
		basic_counters[row] = new uint8_t[width]();

		size_16_counter_alias[row] = (uint16_t*)basic_counters[row];
		size_32_counter_alias[row] = (uint32_t*)basic_counters[row];

		eight_counter_encodings[row] = new uint8_t[width >> 3]();

#ifdef USE_BOBHASH
		bobhash[row].initialize(seed * (7 + row) + row + 100);
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
		seeds[row] = seed * (7 + row) + row + 100;
#endif // USE_XXHASH
	}



	int logSizes[5][4];
	for (int four_encoding = 0; four_encoding < 5; ++four_encoding) {
		for (int counter_number = 0; counter_number < 4; ++counter_number) {
			if (four_encoding == 4) { // all are merged 
				logSizes[four_encoding][counter_number] = 2;
			}
			else if (four_encoding == 3) { // merged in pairs
				logSizes[four_encoding][counter_number] = 1;
			}
			else if (four_encoding == 2) { // left pair is merged
				if (counter_number < 2) {
					logSizes[four_encoding][counter_number] = 1;
				}
				else {
					logSizes[four_encoding][counter_number] = 0;
				}
			}
			else if (four_encoding == 1) {
				if (counter_number < 2) {
					logSizes[four_encoding][counter_number] = 0;
				}
				else {
					logSizes[four_encoding][counter_number] = 1;
				}
			}
			else if (four_encoding == 0) { // none are merged
				logSizes[four_encoding][counter_number] = 0;
			}

		}
	}

	for (int encoding = 0; encoding < 64 * 5; ++encoding)
	{
		int counter_number = encoding & 0b111;
		int residual_encoding = encoding >> 3;
		int first_four_encoding = (residual_encoding & 0b111000) >> 3;
		int second_four_encoding = residual_encoding & 0b111;

		if (counter_number < 4) { //first four
			smart_encoding_lookup_table[encoding] = logSizes[first_four_encoding][counter_number % 4];
		}
		else {//second four
			smart_encoding_lookup_table[encoding] = logSizes[second_four_encoding][counter_number % 4];
		}
		int merges_byte = residual_encoding;
		int smarter_0_encoding = residual_encoding;
		int smarter_1_encoding = residual_encoding;
		if (smart_encoding_lookup_table[encoding] == 0) { // for merging 8b counters
			if (counter_number < 4) { //counters 0-3
				if (first_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 2) { //counters 0-1
						merges_byte += 16;
						//first_four_encoding = 2;
					}
					else {//counters 2-3
						merges_byte += 8;
						//first_four_encoding = 1;
					}
				}
				else { // The other counters were merged
					merges_byte = 24 + (merges_byte & 0b111);
				}
			}
			else { //counters 4-7
				if (second_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 6) { //counters 4-5
						merges_byte += 2;
						//third_four_encoding = 2;
					}
					else {//counters 6-7
						merges_byte += 1;
						//third_four_encoding = 1;
					}
				}
				else { // The other counters were merged
					merges_byte = merges_byte | 0b11;
				}
			}
			smarter_0_encoding = merges_byte;
		}
		else { // for merging 16b counters
			if (counter_number < 4) { //counters 0-3
				merges_byte = 32 | (merges_byte & 0b111);
			}
			else { //counters 4-7
				merges_byte = (merges_byte & 0b111000) + 4;
			}
			smarter_1_encoding = merges_byte;
		}
		//smart_merging_lookup_table[encoding] = merges_byte;
		smarter_merging_lookup_table_0[encoding] = smarter_0_encoding;
		smarter_merging_lookup_table_1[encoding] = smarter_1_encoding;
	}

}

inline uint8_t FastEncodingSALSA::log_counter_size(uint index, uint row)
{
	uint index_over_8 = index >> 3;
	uint index_mod_8 = index & 0b111;

	uint16_t merges_byte = (uint16_t)eight_counter_encodings[row][index_over_8];

	return smart_encoding_lookup_table[merges_byte << 3 | index_mod_8];

	// we assume that counter sizes cannot exceed 4X the original
	// if the original is 8 bits then we have a 32-bit limit

}

#if false
void FastEncodingSALSA::merge_8_bit_counters(uint index, uint row)
{
	uint index_over_8 = index >> 3;
	uint index_mod_8 = index & 0b111;

	uint8_t& merges_byte = eight_counter_encodings[row][index_over_8];
	merges_byte = smarter_merging_lookup_table_0[merges_byte << 3 | index_mod_8];
}

void FastEncodingSALSA::merge_16_bit_counters(uint index, uint row)
{
	uint index_over_8 = index >> 3;
	uint index_mod_8 = index & 0b111;

	uint8_t& merges_byte = eight_counter_encodings[row][index_over_8];
	merges_byte = smarter_merging_lookup_table_1[merges_byte << 3 | index_mod_8];
}
#endif

void FastEncodingSALSA::increment(const char* str)
{

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
#endif // USE_TWO_XXHASHES

	for (int row = 0; row < height; ++row) {

#ifdef USE_BOBHASH
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[row])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
		uint index;
		if (row <= height_plus_1_over_2)
		{
			index = hashes0 & width_mask;
			hashes0 >>= log_width;
		}
		else
		{
			index = hashes1 & width_mask;
			hashes1 >>= log_width;
		}
#endif // USE_TWO_XXHASHES


		//uint8_t log_counter_size_value = log_counter_size(index, row);
		/*
		uint index_over_8 = index >> 3;
		uint index_mod_8 = index & 0b111;

		uint8_t &merges_byte = eight_counter_encodings[row][index_over_8];

		uint8_t log_counter_size_value = 0;
		if (index_mod_8 < 4) {
			uint quadruple = merges_byte >> 3;
			if (quadruple == 0) {

			}
			else if (quadruple == 3) {
				log_counter_size_value = 1;
			}
			else if (quadruple == 4) {
				log_counter_size_value = 2;
			}
			else if (quadruple == 1){
				log_counter_size_value = ((index_mod_8 & 0b11)  > 1);
			}
			else {
				log_counter_size_value = ((index_mod_8 & 0b11) < 2);
			}
		}
		else {
			uint quadruple = merges_byte & 0b111;
			if (quadruple == 0) {

			}
			else if (quadruple == 3) {
				log_counter_size_value = 1;
			}
			else if (quadruple == 4) {
				log_counter_size_value = 2;
			}
			else if (quadruple == 1) {
				log_counter_size_value = ((index_mod_8 & 0b11) > 1);
			}
			else {
				log_counter_size_value = ((index_mod_8 & 0b11) < 2);
			}
		}*/


		uint index_over_8 = index >> 3;
		uint index_mod_8 = index & 0b111;

		uint8_t& merges_byte = eight_counter_encodings[row][index_over_8];
		uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;

		uint8_t log_counter_size_value =  smart_encoding_lookup_table[lookup_key];
#if true
		uint inc_index = index & ~((1 << log_counter_size_value) - 1);

		uint32_t& cnt = *((uint32_t*)(basic_counters[row] + inc_index));
		//++cnt;
		if (++cnt & 0xFF) {
			//continue;
		}
		else {
			if (log_counter_size_value == 0) {
					--cnt;
				uint index_16 = index >> 1;
				size_16_counter_alias[row][index_16] = 1 << 8;
				merges_byte = smarter_merging_lookup_table_0[lookup_key];
			}
			else if ((log_counter_size_value == 1) && (cnt & 0xFFFF) == 0) {
				--cnt;

				uint index_32 = index >> 2;
				size_32_counter_alias[row][index_32] = 1 << 16;
				merges_byte = smarter_merging_lookup_table_1[lookup_key];
			}
			//continue;
		}
#else
		// we are here... 
		if (log_counter_size_value == 0)
		{
			if (++basic_counters[row][index])
			{
				// no overflow
			}
			else
			{
				//merge_8_bit_counters(index, row);
				//uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;
				merges_byte = smarter_merging_lookup_table_0[lookup_key];

				// the current counter became 0 while its real value is 256...
				uint index_16 = index >> 1;
				size_16_counter_alias[row][index_16] = 1 << 8;
			}
		}
		else if (log_counter_size_value == 1)
		{
			uint index_16 = index >> 1;
			if (++size_16_counter_alias[row][index_16])
			{
				// no overflow
			}
			else
			{
				//merge_16_bit_counters(index, row);
				//uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;
				merges_byte = smarter_merging_lookup_table_1[lookup_key];

				// the current counter became 0 while its real value is 2^16...
				uint index_32 = index >> 2;
				size_32_counter_alias[row][index_32] = 1 << 16;
			}
		}
		else if (log_counter_size_value == 2)
		{
			uint index_32 = index >> 2;
			++size_32_counter_alias[row][index_32];
			/*if (++size_32_counter_alias[row][index_32])
			{
				// no overflow
			}
			else
			{
				assert(false);
			}*/
		}
		/*
		else // if (log_counter_size_value == 3)
		{
			assert(false);
		}*/
#endif
	}
}

inline uint64_t FastEncodingSALSA::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);

	switch (log_counter_size_value) {
	case 0:
		return basic_counters[row][index];
	case 1:
		return size_16_counter_alias[row][index >> 1];
	case 2:
		return size_32_counter_alias[row][index >> 2];
	case 3:
		// not implemented
		return -1;
	default:
		// BOB!
		return -1;
	}
}

uint64_t FastEncodingSALSA::query(const char* str)
{

#ifdef USE_BOBHASH
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
	uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[0])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint index = hashes & width_mask;
	hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
	uint index = hashes0 & width_mask;
	hashes0 >>= log_width;
#endif // USE_TWO_XXHASHES

	int row = 0;

	uint64_t min = read_counter(index, row++);

	while (row < height) {

#ifdef USE_BOBHASH
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[row])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
		uint index;
		if (row <= height_plus_1_over_2)
		{
			index = hashes0 & width_mask;
			hashes0 >>= log_width;
		}
		else
		{
			index = hashes1 & width_mask;
			hashes1 >>= log_width;
		}
#endif 
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

StupidEncodingSALSA::StupidEncodingSALSA()
{
}

StupidEncodingSALSA::~StupidEncodingSALSA()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] basic_counters[i];
		delete[] eight_counter_encodings[i];
	}

#ifdef USE_BOBHASH
	delete[] bobhash;
#endif // USE_BOBHASH

	delete[] basic_counters;

	delete[] size_16_counter_alias;
	delete[] size_32_counter_alias;

	delete[] eight_counter_encodings;
}

void StupidEncodingSALSA::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	int index = width_mask;
	log_width = 0;
	while (index >>= 1) ++log_width;

	assert(width > 0 && "We assume too much!");
	assert(width % 16 == 0 && "We assume that (w % 16 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

#ifdef USE_ONE_XXHASH
	if (height > 1) {
		assert(width <= ((uint64_t)1 << (64 / height)) && "One XXHASH is not enough!");
	}
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	height_plus_1_over_2 = (height + 1) / 2;
	assert(width < ((uint64_t)1 << (64 / height_plus_1_over_2)) && "Two XXHASHes are not enough!");
#endif // USE_TWO_XXHASHES

	basic_counters = new uint8_t * [height];

	size_16_counter_alias = new uint16_t * [height];
	size_32_counter_alias = new uint32_t * [height];

	eight_counter_encodings = new uint8_t * [height];

#ifdef USE_BOBHASH
	bobhash = new BOBHash[height];
#endif // USE_BOBHASH

	for (int row = 0; row < height; ++row)
	{
		basic_counters[row] = new uint8_t[width]();

		size_16_counter_alias[row] = (uint16_t*)basic_counters[row];
		size_32_counter_alias[row] = (uint32_t*)basic_counters[row];

		eight_counter_encodings[row] = new uint8_t[width >> 3]();

#ifdef USE_BOBHASH
		bobhash[row].initialize(seed * (7 + row) + row + 100);
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
		seeds[row] = seed * (7 + row) + row + 100;
#endif // USE_XXHASH
	}


	/*
	int logSizes[5][4];
	for (int four_encoding = 0; four_encoding < 5; ++four_encoding) {
		for (int counter_number = 0; counter_number < 4; ++counter_number) {
			if (four_encoding == 4) { // all are merged 
				logSizes[four_encoding][counter_number] = 2;
			}
			else if (four_encoding == 3) { // merged in pairs
				logSizes[four_encoding][counter_number] = 1;
			}
			else if (four_encoding == 2) { // left pair is merged
				if (counter_number < 2) {
					logSizes[four_encoding][counter_number] = 1;
				}
				else {
					logSizes[four_encoding][counter_number] = 0;
				}
			}
			else if (four_encoding == 1) {
				if (counter_number < 2) {
					logSizes[four_encoding][counter_number] = 0;
				}
				else {
					logSizes[four_encoding][counter_number] = 1;
				}
			}
			else if (four_encoding == 0) { // none are merged
				logSizes[four_encoding][counter_number] = 0;
			}

		}
	}

	for (int encoding = 0; encoding < 64 * 5; ++encoding)
	{
		int counter_number = encoding & 0b111;
		int residual_encoding = encoding >> 3;
		int first_four_encoding = (residual_encoding & 0b111000) >> 3;
		int second_four_encoding = residual_encoding & 0b111;

		if (counter_number < 4) { //first four
			smart_encoding_lookup_table[encoding] = logSizes[first_four_encoding][counter_number % 4];
		}
		else {//second four
			smart_encoding_lookup_table[encoding] = logSizes[second_four_encoding][counter_number % 4];
		}
		int merges_byte = residual_encoding;
		int smarter_0_encoding = residual_encoding;
		int smarter_1_encoding = residual_encoding;
		if (smart_encoding_lookup_table[encoding] == 0) { // for merging 8b counters
			if (counter_number < 4) { //counters 0-3
				if (first_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 2) { //counters 0-1
						merges_byte += 16;
						//first_four_encoding = 2;
					}
					else {//counters 2-3
						merges_byte += 8;
						//first_four_encoding = 1;
					}
				}
				else { // The other counters were merged
					merges_byte = 24 + (merges_byte & 0b111);
				}
			}
			else { //counters 4-7
				if (second_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 6) { //counters 4-5
						merges_byte += 2;
						//third_four_encoding = 2;
					}
					else {//counters 6-7
						merges_byte += 1;
						//third_four_encoding = 1;
					}
				}
				else { // The other counters were merged
					merges_byte = merges_byte | 0b11;
				}
			}
			smarter_0_encoding = merges_byte;
		}
		else { // for merging 16b counters
			if (counter_number < 4) { //counters 0-3
				merges_byte = 32 | (merges_byte & 0b111);
			}
			else { //counters 4-7
				merges_byte = (merges_byte & 0b111000) + 4;
			}
			smarter_1_encoding = merges_byte;
		}
		//smart_merging_lookup_table[encoding] = merges_byte;
		smarter_merging_lookup_table_0[encoding] = smarter_0_encoding;
		smarter_merging_lookup_table_1[encoding] = smarter_1_encoding;
	}
	*/
}

inline uint8_t StupidEncodingSALSA::log_counter_size(uint index, uint row)
{
	uint index_over_8 = index >> 3;
	//uint index_mod_8 = index & 0b111;

	uint index_mod_8_over_2 = (index & 0b111) >> 1;

	uint8_t merges_byte = (uint16_t)eight_counter_encodings[row][index_over_8];

	return (merges_byte & (0b11 << (index_mod_8_over_2 << 1))) >> (index_mod_8_over_2 << 1);

	// we assume that counter sizes cannot exceed 4X the original
	// if the original is 8 bits then we have a 32-bit limit

}

#if false
void StupidEncodingSALSA::merge_8_bit_counters(uint index, uint row)
{
	uint index_over_8 = index >> 3;
	uint index_mod_8 = index & 0b111;

	uint8_t& merges_byte = eight_counter_encodings[row][index_over_8];
	merges_byte = smarter_merging_lookup_table_0[merges_byte << 3 | index_mod_8];
}

void StupidEncodingSALSA::merge_16_bit_counters(uint index, uint row)
{
	uint index_over_8 = index >> 3;
	uint index_mod_8 = index & 0b111;

	uint8_t& merges_byte = eight_counter_encodings[row][index_over_8];
	merges_byte = smarter_merging_lookup_table_1[merges_byte << 3 | index_mod_8];
}
#endif

void StupidEncodingSALSA::increment(const char* str)
{

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
#endif // USE_TWO_XXHASHES

	for (int row = 0; row < height; ++row) {

#ifdef USE_BOBHASH
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[row])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
		uint index;
		if (row <= height_plus_1_over_2)
		{
			index = hashes0 & width_mask;
			hashes0 >>= log_width;
		}
		else
		{
			index = hashes1 & width_mask;
			hashes1 >>= log_width;
		}
#endif // USE_TWO_XXHASHES


		//uint8_t log_counter_size_value = log_counter_size(index, row);
		/*
		uint index_over_8 = index >> 3;
		uint index_mod_8 = index & 0b111;

		uint8_t &merges_byte = eight_counter_encodings[row][index_over_8];

		uint8_t log_counter_size_value = 0;
		if (index_mod_8 < 4) {
			uint quadruple = merges_byte >> 3;
			if (quadruple == 0) {

			}
			else if (quadruple == 3) {
				log_counter_size_value = 1;
			}
			else if (quadruple == 4) {
				log_counter_size_value = 2;
			}
			else if (quadruple == 1){
				log_counter_size_value = ((index_mod_8 & 0b11)  > 1);
			}
			else {
				log_counter_size_value = ((index_mod_8 & 0b11) < 2);
			}
		}
		else {
			uint quadruple = merges_byte & 0b111;
			if (quadruple == 0) {

			}
			else if (quadruple == 3) {
				log_counter_size_value = 1;
			}
			else if (quadruple == 4) {
				log_counter_size_value = 2;
			}
			else if (quadruple == 1) {
				log_counter_size_value = ((index_mod_8 & 0b11) > 1);
			}
			else {
				log_counter_size_value = ((index_mod_8 & 0b11) < 2);
			}
		}*/


		uint index_over_8 = index >> 3;
		//uint index_mod_8 = index & 0b111;
		//uint index_mod_8_over_2 = (index & 0b111) >> 1;
		uint index_mod_8_mod_2 = index & 0b110;

		uint8_t& merges_byte = eight_counter_encodings[row][index_over_8];
		//uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;

		uint64_t log_counter_size_value = (merges_byte >> index_mod_8_mod_2) & 0b11;//(merges_byte & (0b11 << index_mod_8_mod_2)) >> index_mod_8_mod_2; //smart_encoding_lookup_table[lookup_key];
#if true
		//uint inc_index = index & ~((1 << log_counter_size_value) - 1);
		//uint inc_index = (log_counter_size_value == 1) ? index & 0xFFFE : (log_counter_size_value == 0) ? index : index & 0xFFFC;
		uint inc_index = (log_counter_size_value == 0) ? index  : (log_counter_size_value == 1) ? index & 0xFFFFFFFE : index & 0xFFFFFFFC;
		//uint8_t* &basic_counters_row = basic_counters[row];
		//if (++basic_counters_row[inc_index]) {
		if (++basic_counters[row][inc_index]) {
			//continue;
		}
		else {
			//uint8_t& merges_byte = eight_counter_encodings[row][index_over_8];
			if (log_counter_size_value == 1) {
				if (++basic_counters[row][inc_index+1]) {
				//if (++basic_counters_row[inc_index + 1] & 0xFF) {
				}
				else {
					uint index_32 = index >> 2;
					size_32_counter_alias[row][index_32] = 1 << 16;
					merges_byte &= ~(0b11 << index_mod_8_mod_2);
					merges_byte |= 2 << index_mod_8_mod_2;

					int paired_index_mod_8_mod_2 = index_mod_8_mod_2 ^ 0b10;

					merges_byte &= ~(0b11 << paired_index_mod_8_mod_2);
					merges_byte |= 2 << paired_index_mod_8_mod_2;
				}

			}
			else if (log_counter_size_value == 0) {
				//uint index_16 = index >> 1;
				//size_16_counter_alias[row][index_16] = 1 << 8;
				basic_counters[row][inc_index | 1] = 1;
				basic_counters[row][inc_index & 0xFFFFFFFE] = 0;
				merges_byte = (merges_byte | (1 << index_mod_8_mod_2));
			}
			else {
				if (++basic_counters[row][inc_index + 1]) {
				}
				else {
					++basic_counters[row][inc_index + 2];
				}
			}

			//continue;
		}
#else
		// we are here... 
		if (log_counter_size_value == 0)
		{
			if (++basic_counters[row][index])
			{
				// no overflow
			}
			else
			{
				//merge_8_bit_counters(index, row);
				//uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;
				merges_byte = (merges_byte | (1 << index_mod_8_mod_2));

				// the current counter became 0 while its real value is 256...
				uint index_16 = index >> 1;
				size_16_counter_alias[row][index_16] = 1 << 8;
			}
		}
		else if (log_counter_size_value == 1)
		{
			uint index_16 = index >> 1;
			if (++size_16_counter_alias[row][index_16])
			{
				// no overflow
			}
			else
			{
				//merge_16_bit_counters(index, row);
				//uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;
				merges_byte &= ~(0b11 << (index_mod_8_over_2 << 1));
				merges_byte |= 2 << (index_mod_8_over_2 << 1);

				int paired_index_mod_8_over_2 = index_mod_8_over_2 ^ 1;

				merges_byte &= ~(0b11 << (paired_index_mod_8_over_2 << 1));
				merges_byte |= 2 << (paired_index_mod_8_over_2 << 1);

				// the current counter became 0 while its real value is 2^16...
				uint index_32 = index >> 2;
				size_32_counter_alias[row][index_32] = 1 << 16;
			}
		}
		else if (log_counter_size_value == 2)
		{
			uint index_32 = index >> 2;
			++size_32_counter_alias[row][index_32];
			/*if (++size_32_counter_alias[row][index_32])
			{
				// no overflow
			}
			else
			{
				assert(false);
			}*/
		}
		/*
		else // if (log_counter_size_value == 3)
		{
			assert(false);
		}*/
#endif
	}
}

inline uint64_t StupidEncodingSALSA::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);

	switch (log_counter_size_value) {
	case 0:
		return basic_counters[row][index];
	case 1:
		return size_16_counter_alias[row][index >> 1];
	case 2:
		return size_32_counter_alias[row][index >> 2];
	case 3:
		// not implemented
		return -1;
	default:
		// BOB!
		return -1;
	}
}

uint64_t StupidEncodingSALSA::query(const char* str)
{

#ifdef USE_BOBHASH
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
	uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[0])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint index = hashes & width_mask;
	hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
	uint index = hashes0 & width_mask;
	hashes0 >>= log_width;
#endif // USE_TWO_XXHASHES

	int row = 0;

	uint64_t min = read_counter(index, row++);

	while (row < height) {

#ifdef USE_BOBHASH
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[row])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
		uint index;
		if (row <= height_plus_1_over_2)
		{
			index = hashes0 & width_mask;
			hashes0 >>= log_width;
		}
		else
		{
			index = hashes1 & width_mask;
			hashes1 >>= log_width;
		}
#endif 
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

StupiderEncodingSALSA::StupiderEncodingSALSA()
{
}

StupiderEncodingSALSA::~StupiderEncodingSALSA()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] basic_counters[i];
		delete[] sixty_four_counter_encodings[i];
	}

#ifdef USE_BOBHASH
	delete[] bobhash;
#endif // USE_BOBHASH

	delete[] basic_counters;

	delete[] size_16_counter_alias;
	delete[] size_32_counter_alias;

	delete[] sixty_four_counter_encodings;
}

void StupiderEncodingSALSA::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	int index = width_mask;
	log_width = 0;
	while (index >>= 1) ++log_width;

	assert(width > 0 && "We assume too much!");
	assert(width % 16 == 0 && "We assume that (w % 16 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

#ifdef USE_ONE_XXHASH
	if (height > 1) {
		assert(width <= ((uint64_t)1 << (64 / height)) && "One XXHASH is not enough!");
	}
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	height_plus_1_over_2 = (height + 1) / 2;
	assert(width < ((uint64_t)1 << (64 / height_plus_1_over_2)) && "Two XXHASHes are not enough!");
#endif // USE_TWO_XXHASHES

	basic_counters = new uint8_t * [height];

	size_16_counter_alias = new uint16_t * [height];
	size_32_counter_alias = new uint32_t * [height];

	sixty_four_counter_encodings = new uint64_t * [height];

#ifdef USE_BOBHASH
	bobhash = new BOBHash[height];
#endif // USE_BOBHASH

	for (int row = 0; row < height; ++row)
	{
		basic_counters[row] = new uint8_t[width]();

		size_16_counter_alias[row] = (uint16_t*)basic_counters[row];
		size_32_counter_alias[row] = (uint32_t*)basic_counters[row];

		sixty_four_counter_encodings[row] = new uint64_t[width >> 6]();

#ifdef USE_BOBHASH
		bobhash[row].initialize(seed * (7 + row) + row + 100);
#endif // USE_BOBHASH

#if defined USE_XXHASH || defined USE_ONE_XXHASH || defined USE_TWO_XXHASHES
		seeds[row] = seed * (7 + row) + row + 100;
#endif // USE_XXHASH
	}


	/*
	int logSizes[5][4];
	for (int four_encoding = 0; four_encoding < 5; ++four_encoding) {
		for (int counter_number = 0; counter_number < 4; ++counter_number) {
			if (four_encoding == 4) { // all are merged
				logSizes[four_encoding][counter_number] = 2;
			}
			else if (four_encoding == 3) { // merged in pairs
				logSizes[four_encoding][counter_number] = 1;
			}
			else if (four_encoding == 2) { // left pair is merged
				if (counter_number < 2) {
					logSizes[four_encoding][counter_number] = 1;
				}
				else {
					logSizes[four_encoding][counter_number] = 0;
				}
			}
			else if (four_encoding == 1) {
				if (counter_number < 2) {
					logSizes[four_encoding][counter_number] = 0;
				}
				else {
					logSizes[four_encoding][counter_number] = 1;
				}
			}
			else if (four_encoding == 0) { // none are merged
				logSizes[four_encoding][counter_number] = 0;
			}

		}
	}

	for (int encoding = 0; encoding < 64 * 5; ++encoding)
	{
		int counter_number = encoding & 0b111;
		int residual_encoding = encoding >> 3;
		int first_four_encoding = (residual_encoding & 0b111000) >> 3;
		int second_four_encoding = residual_encoding & 0b111;

		if (counter_number < 4) { //first four
			smart_encoding_lookup_table[encoding] = logSizes[first_four_encoding][counter_number % 4];
		}
		else {//second four
			smart_encoding_lookup_table[encoding] = logSizes[second_four_encoding][counter_number % 4];
		}
		int merges_byte = residual_encoding;
		int smarter_0_encoding = residual_encoding;
		int smarter_1_encoding = residual_encoding;
		if (smart_encoding_lookup_table[encoding] == 0) { // for merging 8b counters
			if (counter_number < 4) { //counters 0-3
				if (first_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 2) { //counters 0-1
						merges_byte += 16;
						//first_four_encoding = 2;
					}
					else {//counters 2-3
						merges_byte += 8;
						//first_four_encoding = 1;
					}
				}
				else { // The other counters were merged
					merges_byte = 24 + (merges_byte & 0b111);
				}
			}
			else { //counters 4-7
				if (second_four_encoding == 0) { // all counters were 8-bits
					if (counter_number < 6) { //counters 4-5
						merges_byte += 2;
						//third_four_encoding = 2;
					}
					else {//counters 6-7
						merges_byte += 1;
						//third_four_encoding = 1;
					}
				}
				else { // The other counters were merged
					merges_byte = merges_byte | 0b11;
				}
			}
			smarter_0_encoding = merges_byte;
		}
		else { // for merging 16b counters
			if (counter_number < 4) { //counters 0-3
				merges_byte = 32 | (merges_byte & 0b111);
			}
			else { //counters 4-7
				merges_byte = (merges_byte & 0b111000) + 4;
			}
			smarter_1_encoding = merges_byte;
		}
		//smart_merging_lookup_table[encoding] = merges_byte;
		smarter_merging_lookup_table_0[encoding] = smarter_0_encoding;
		smarter_merging_lookup_table_1[encoding] = smarter_1_encoding;
	}
	*/
}

inline uint8_t StupiderEncodingSALSA::log_counter_size(uint index, uint row)
{
	uint index_over_64 = index >> 6;
	//uint index_mod_8 = index & 0b111;

	//uint index_mod_8_over_2 = (index & 0b111) >> 1;
	uint index_mod_64_over_2 = (index & 0b111111) >> 1;

	//uint8_t merges_byte = (uint16_t)eight_counter_encodings[row][index_over_8];
	uint64_t merges_word = sixty_four_counter_encodings[row][index_over_64];

	return (merges_word & ((uint64_t)0b11 << (index_mod_64_over_2 << 1))) >> (index_mod_64_over_2 << 1);

	// we assume that counter sizes cannot exceed 4X the original
	// if the original is 8 bits then we have a 32-bit limit

}

void StupiderEncodingSALSA::increment(const char* str)
{

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
#endif // USE_TWO_XXHASHES

	for (int row = 0; row < height; ++row) {

#ifdef USE_BOBHASH
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[row])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
		uint index;
		if (row <= height_plus_1_over_2)
		{
			index = hashes0 & width_mask;
			hashes0 >>= log_width;
		}
		else
		{
			index = hashes1 & width_mask;
			hashes1 >>= log_width;
		}
#endif // USE_TWO_XXHASHES


		//uint8_t log_counter_size_value = log_counter_size(index, row);
		/*
		uint index_over_8 = index >> 3;
		uint index_mod_8 = index & 0b111;

		uint8_t &merges_byte = eight_counter_encodings[row][index_over_8];

		uint8_t log_counter_size_value = 0;
		if (index_mod_8 < 4) {
			uint quadruple = merges_byte >> 3;
			if (quadruple == 0) {

			}
			else if (quadruple == 3) {
				log_counter_size_value = 1;
			}
			else if (quadruple == 4) {
				log_counter_size_value = 2;
			}
			else if (quadruple == 1){
				log_counter_size_value = ((index_mod_8 & 0b11)  > 1);
			}
			else {
				log_counter_size_value = ((index_mod_8 & 0b11) < 2);
			}
		}
		else {
			uint quadruple = merges_byte & 0b111;
			if (quadruple == 0) {

			}
			else if (quadruple == 3) {
				log_counter_size_value = 1;
			}
			else if (quadruple == 4) {
				log_counter_size_value = 2;
			}
			else if (quadruple == 1) {
				log_counter_size_value = ((index_mod_8 & 0b11) > 1);
			}
			else {
				log_counter_size_value = ((index_mod_8 & 0b11) < 2);
			}
		}*/


		//uint index_over_8 = index >> 3;
		//uint index_mod_8 = index & 0b111;
		//uint index_mod_8_over_2 = (index & 0b111) >> 1;
		//uint index_mod_8_mod_2 = index & 0b110;

		//uint64_t merges_byte = eight_counter_encodings[row][index_over_8];
		//uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;

		uint index_over_64 = index >> 6;
		//uint index_mod_8 = index & 0b111;

		//uint index_mod_8_over_2 = (index & 0b111) >> 1;
		uint index_mod_8_mod_2 = index & 0b111110;

		//uint8_t merges_byte = (uint16_t)eight_counter_encodings[row][index_over_8];
		uint64_t merges_word = sixty_four_counter_encodings[row][index_over_64];

		uint64_t log_counter_size_value = (merges_word >> index_mod_8_mod_2) & 0b11;//(merges_byte & (0b11 << index_mod_8_mod_2)) >> index_mod_8_mod_2; //smart_encoding_lookup_table[lookup_key];
#if true
		//uint inc_index = index & ~((1 << log_counter_size_value) - 1);
		//uint inc_index = (log_counter_size_value == 1) ? index & 0xFFFE : (log_counter_size_value == 0) ? index : index & 0xFFFC;
		uint inc_index = (log_counter_size_value == 0) ? index : (log_counter_size_value == 1) ? index & 0xFFFFFFFE : index & 0xFFFFFFFC;
		//uint8_t* &basic_counters_row = basic_counters[row];
		//if (++basic_counters_row[inc_index]) {
		if (++basic_counters[row][inc_index]) {
			//continue;
		}
		else {
			uint64_t &merges_word = sixty_four_counter_encodings[row][index_over_64];
			if (log_counter_size_value == 1) {
				if (++basic_counters[row][inc_index + 1]) {
					//if (++basic_counters_row[inc_index + 1] & 0xFF) {
				}
				else {
					uint index_32 = index >> 2;
					size_32_counter_alias[row][index_32] = 1 << 16;
					merges_word &= ~((uint64_t) 0b11 << index_mod_8_mod_2);
					merges_word |= (uint64_t)2 << index_mod_8_mod_2;

					int paired_index_mod_64_mod_2 = index_mod_8_mod_2 ^ 0b10;

					merges_word &= ~((uint64_t)0b11 << paired_index_mod_64_mod_2);
					merges_word |= (uint64_t)2 << paired_index_mod_64_mod_2;
				}

			}
			else if (log_counter_size_value == 0) {
				//uint index_16 = index >> 1;
				//size_16_counter_alias[row][index_16] = 1 << 8;
				basic_counters[row][inc_index | 1] = 1;
				basic_counters[row][inc_index & 0xFFFFFFFE] = 0;
				merges_word = (merges_word | ((uint64_t)1 << index_mod_8_mod_2));
			}
			else {
				if (++basic_counters[row][inc_index + 1]) {
				}
				else {
					++basic_counters[row][inc_index + 2];
				}
			}

			//continue;
		}
#else
		// we are here... 
		if (log_counter_size_value == 0)
		{
			if (++basic_counters[row][index])
			{
				// no overflow
			}
			else
			{
				//merge_8_bit_counters(index, row);
				//uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;
				merges_byte = (merges_byte | (1 << index_mod_8_mod_2));

				// the current counter became 0 while its real value is 256...
				uint index_16 = index >> 1;
				size_16_counter_alias[row][index_16] = 1 << 8;
			}
		}
		else if (log_counter_size_value == 1)
		{
			uint index_16 = index >> 1;
			if (++size_16_counter_alias[row][index_16])
			{
				// no overflow
			}
			else
			{
				//merge_16_bit_counters(index, row);
				//uint16_t lookup_key = ((uint16_t)merges_byte << 3) | index_mod_8;
				merges_byte &= ~(0b11 << (index_mod_8_over_2 << 1));
				merges_byte |= 2 << (index_mod_8_over_2 << 1);

				int paired_index_mod_8_over_2 = index_mod_8_over_2 ^ 1;

				merges_byte &= ~(0b11 << (paired_index_mod_8_over_2 << 1));
				merges_byte |= 2 << (paired_index_mod_8_over_2 << 1);

				// the current counter became 0 while its real value is 2^16...
				uint index_32 = index >> 2;
				size_32_counter_alias[row][index_32] = 1 << 16;
			}
		}
		else if (log_counter_size_value == 2)
		{
			uint index_32 = index >> 2;
			++size_32_counter_alias[row][index_32];
			/*if (++size_32_counter_alias[row][index_32])
			{
				// no overflow
			}
			else
			{
				assert(false);
			}*/
		}
		/*
		else // if (log_counter_size_value == 3)
		{
			assert(false);
		}*/
#endif
	}
}

inline uint64_t StupiderEncodingSALSA::read_counter(uint index, uint row)
{
	uint8_t log_counter_size_value = log_counter_size(index, row);

	switch (log_counter_size_value) {
	case 0:
		return basic_counters[row][index];
	case 1:
		return size_16_counter_alias[row][index >> 1];
	case 2:
		return size_32_counter_alias[row][index >> 2];
	case 3:
		// not implemented
		return -1;
	default:
		// BOB!
		return -1;
	}
}

uint64_t StupiderEncodingSALSA::query(const char* str)
{

#ifdef USE_BOBHASH
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
	uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[0])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
	uint64_t hashes = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint index = hashes & width_mask;
	hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
	uint64_t hashes0 = xxh::xxhash3<64>(str, FT_SIZE, seeds[0]);
	uint64_t hashes1 = xxh::xxhash3<64>(str, FT_SIZE, seeds[1]);
	uint index = hashes0 & width_mask;
	hashes0 >>= log_width;
#endif // USE_TWO_XXHASHES

	int row = 0;

	uint64_t min = read_counter(index, row++);

	while (row < height) {

#ifdef USE_BOBHASH
		uint index = (bobhash[row].run(str, FT_SIZE)) & width_mask;
#endif // USE_BOBHASH

#ifdef USE_XXHASH
		uint index = (xxh::xxhash3<64>(str, FT_SIZE, seeds[row])) & width_mask;
#endif // USE_XXHASH

#ifdef USE_ONE_XXHASH
		uint index = hashes & width_mask;
		hashes >>= log_width;
#endif // USE_ONE_XXHASH

#ifdef USE_TWO_XXHASHES
		uint index;
		if (row <= height_plus_1_over_2)
		{
			index = hashes0 & width_mask;
			hashes0 >>= log_width;
		}
		else
		{
			index = hashes1 & width_mask;
			hashes1 >>= log_width;
		}
#endif 
		uint64_t temp = read_counter(index, row++);
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}