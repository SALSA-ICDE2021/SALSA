#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "CountMinAEE.hpp"

using namespace std;


///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

inline unsigned CountMinCC::accelerated_geometric(double LogOneMinusP) {
	// Return a randomly drawn geometric variable with parameter p, such that LogOneMinusP = log(1-p)
	return log(dis_arr(generator_arr)) / LogOneMinusP + 1;
}

CountMinCC::CountMinCC(int width, int height, int seed) :

	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX)),
	RandByteArray((char *)&rand_bits),
	weighted_rand_bits(_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr())),
	WeightedRandByteArray((uint16_t *)&weighted_rand_bits)
{

	this->seed = seed;
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	cc_cms = new cc_counter_t*[height];
	for (int i = 0; i < height; ++i)
	{
		cc_cms[i] = new cc_counter_t[width]();
	}

	bobhash = new BOBHash[height];
	for (int i = 0; i < height; ++i)
	{
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}

	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");

	width_over_four = width >> 2;

	indices = new int[height];
}

CountMinCC::~CountMinCC()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] cc_cms[i];
	}

	delete[] cc_cms;
	delete[] bobhash;
	delete[] indices;
}

void CountMinCC::divide_array_counters()
{
	for (int i = 0; i < height; ++i)
	{
		uint64_t* four_counters_array = (uint64_t*)cc_cms[i];
		for (unsigned j = 0; j < width_over_four; ++j)
		{
			uint64_t &four_counters = four_counters_array[j];
			four_counters = (four_counters & division_mask) >> 1;
		}
	}
}

void CountMinCC::increment(const char * str)
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

void CountMinCC::increment_h(const char * str)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	for (int i = 0; i < height; ++i) 
	{
		cc_counter_t &comp_counter = cc_cms[i][bobhash[i].run(str, FT_SIZE) % width];
		if (comp_counter == 65535)
		{
			++minus_log_p;
			LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
			divide_array_counters();

			// if p decreases we need to reduce sampling rate accordingly.
			rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
		}
		comp_counter += 1; 
	}
}

void CountMinCC::add(const char * str, int c)
{
	if (minus_log_p == 0)
	{
		add_h(str, c);
		return;
	}

	else if (minus_log_p < minus_log_p_bound)
	{

		if (weighted_rand_pos < 16)
		{
		}
		else
		{
			weighted_rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
			WeightedRandByteArray = (uint16_t *)&weighted_rand_bits;
			weighted_rand_pos = 0;
		}

		w1 = c >> minus_log_p;
		w2 = c - (w1 << minus_log_p);

		uint32_t increment = w1;
		
		if (WeightedRandByteArray[weighted_rand_pos] <= w2 << (16 - minus_log_p))
		{
			increment += 1;
		}
		++weighted_rand_pos;

		add_h(str, increment);
	}

	else
	{		
		if (minus_log_p >= 15) 
		{
			uint32_t increment = 0;
			while (w >= interval)
			{
				++increment;
				w -= interval;
				interval = accelerated_geometric(LogOneMinusP);
			}
			interval -= w;

			if (increment == 0)
			{
			}
			else
			{
				add_h(str, increment);
			}	
		}		
		else
		{
			w1 = c >> minus_log_p;
			w2 = c - (w1 << minus_log_p);

			uint32_t increment = w1;

			while (w2 >= interval)
			{
				++increment;
				w2 -= interval;
				interval = accelerated_geometric(LogOneMinusP);
			}
			interval -= w2;

			if (increment == 0)
			{
			}
			else
			{
				add_h(str, increment);
			}
		}		
	}
}

void CountMinCC::add_h(const char * str, int c)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	indices[0] = bobhash[0].run(str, FT_SIZE) % width;
	cc_counter_t max = cc_cms[0][indices[0]];
	for (int i = 1; i < height; ++i)
	{
		indices[i] = bobhash[i].run(str, FT_SIZE) % width;
		if (max < cc_cms[i][indices[i]])
		{
			max = cc_cms[i][indices[i]];
		}
	}

	if (max > 65535 - c)
	{
		// TODO: devide p several times?
		++minus_log_p;
		LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
		interval = accelerated_geometric(LogOneMinusP);

		for (int i = 0; i < height; ++i)
		{
			for (int j = 0; j < width; ++j)
			{
				if (j != indices[i]) 
				{
					cc_cms[i][j] /= 2;
				}
				else
				{
					curr_sum = cc_cms[i][j] + c;
					cc_cms[i][j] = curr_sum / 2;
				}
			}
		}
		return;
	}
	else
	{
		for (int i = 0; i < height; ++i)
		{
			cc_cms[i][indices[i]] += c;
		}
	}
}

uint64_t CountMinCC::query(const char * str)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	cc_counter_t min = cc_cms[0][bobhash[0].run(str, FT_SIZE) % width];
	for (int i = 1; i < height; ++i)
	{
		cc_counter_t &comp_counter = cc_cms[i][bobhash[i].run(str, FT_SIZE) % width];
		if (min > comp_counter)
		{
			min = comp_counter;
		}
	}
	return ((uint64_t)min) << minus_log_p;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

inline unsigned CountMinCCSanity::accelerated_geometric(double LogOneMinusP) {
	// Return a randomly drawn geometric variable with parameter p, such that LogOneMinusP = log(1-p)
	return log(dis_arr(generator_arr)) / LogOneMinusP + 1;
}

CountMinCCSanity::CountMinCCSanity(int width, int height, int seed) :

	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX)),
	RandByteArray((char *)&rand_bits),
	weighted_rand_bits(_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr())),
	WeightedRandByteArray((uint16_t *)&weighted_rand_bits)
{

	numPackets = 0;

	this->seed = seed;
	this->width = width;
	this->height = height;

	width_mask = (width - 1) << 4;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	cc_cms = new cc_counter_t*[height];
	for (int i = 0; i < height; ++i)
	{
		cc_cms[i] = new cc_counter_t[width]();
	}

	bobhash = new BOBHash[height];
	for (int i = 0; i < height; ++i)
	{
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}

	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");

	width_over_four = width >> 2;

	indices = new int[height];
}

CountMinCCSanity::~CountMinCCSanity()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] cc_cms[i];
	}

	delete[] cc_cms;
	delete[] bobhash;
	delete[] indices;
}

void CountMinCCSanity::divide_array_counters()
{
	for (int i = 0; i < height; ++i)
	{
		uint64_t* four_counters_array = (uint64_t*)cc_cms[i];
		for (unsigned j = 0; j < width_over_four; ++j)
		{
			uint64_t &four_counters = four_counters_array[j];
			four_counters = (four_counters & division_mask) >> 1;
		}
	}
}

void CountMinCCSanity::increment(const char * str)
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

void CountMinCCSanity::increment_h(const char * str)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	for (int i = 0; i < height; ++i)
	{
		cc_counter_t &comp_counter = cc_cms[i][(bobhash[i].run(str, FT_SIZE) & width_mask) >> 4];
		if (comp_counter == 65535)
		{
			++minus_log_p;
			LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
			divide_array_counters();

			// if p decreases we need to reduce sampling rate accordingly.
			rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
		}
		comp_counter += 1;
	}
}

void CountMinCCSanity::add(const char * str, int c)
{
	if (minus_log_p == 0)
	{
		add_h(str, c);
		return;
	}

	else if (minus_log_p < minus_log_p_bound)
	{

		if (weighted_rand_pos < 16)
		{
		}
		else
		{
			weighted_rand_bits = _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr());
			WeightedRandByteArray = (uint16_t *)&weighted_rand_bits;
			weighted_rand_pos = 0;
		}

		w1 = c >> minus_log_p;
		w2 = c - (w1 << minus_log_p);

		uint32_t increment = w1;

		if (WeightedRandByteArray[weighted_rand_pos] <= w2 << (16 - minus_log_p))
		{
			increment += 1;
		}
		++weighted_rand_pos;

		add_h(str, increment);
	}

	else
	{
		if (minus_log_p >= 15)
		{
			uint32_t increment = 0;
			while (w >= interval)
			{
				++increment;
				w -= interval;
				interval = accelerated_geometric(LogOneMinusP);
			}
			interval -= w;

			if (increment == 0)
			{
			}
			else
			{
				add_h(str, increment);
			}
		}
		else
		{
			w1 = c >> minus_log_p;
			w2 = c - (w1 << minus_log_p);

			uint32_t increment = w1;

			while (w2 >= interval)
			{
				++increment;
				w2 -= interval;
				interval = accelerated_geometric(LogOneMinusP);
			}
			interval -= w2;

			if (increment == 0)
			{
			}
			else
			{
				add_h(str, increment);
			}
		}
	}
}

void CountMinCCSanity::add_h(const char * str, int c)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	indices[0] = bobhash[0].run(str, FT_SIZE) % width;
	cc_counter_t max = cc_cms[0][indices[0]];
	for (int i = 1; i < height; ++i)
	{
		indices[i] = bobhash[i].run(str, FT_SIZE) % width;
		if (max < cc_cms[i][indices[i]])
		{
			max = cc_cms[i][indices[i]];
		}
	}

	if (max > 65535 - c)
	{
		// TODO: devide p several times?
		++minus_log_p;
		LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
		interval = accelerated_geometric(LogOneMinusP);

		for (int i = 0; i < height; ++i)
		{
			for (int j = 0; j < width; ++j)
			{
				if (j != indices[i])
				{
					cc_cms[i][j] /= 2;
				}
				else
				{
					curr_sum = cc_cms[i][j] + c;
					cc_cms[i][j] = curr_sum / 2;
				}
			}
		}
		return;
	}
	else
	{
		for (int i = 0; i < height; ++i)
		{
			cc_cms[i][indices[i]] += c;
		}
	}
}

uint64_t CountMinCCSanity::query(const char * str)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	cc_counter_t min = cc_cms[0][(bobhash[0].run(str, FT_SIZE) & width_mask) >> 4];
	for (int i = 1; i < height; ++i)
	{
		cc_counter_t &comp_counter = cc_cms[i][(bobhash[i].run(str, FT_SIZE) & width_mask) >> 4];
		if (min > comp_counter)
		{
			min = comp_counter;
		}
	}
	return ((uint64_t)min) << minus_log_p;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

inline unsigned ConservativeUpdateCC::accelerated_geometric(double LogOneMinusP) {
	// Return a randomly drawn geometric variable with parameter p, such that LogOneMinusP = log(1-p)
	return log(dis_arr(generator_arr)) / LogOneMinusP + 1;
}

ConservativeUpdateCC::ConservativeUpdateCC(int width, int height, int seed) :

	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr())),
	RandByteArray((char *)&rand_bits)

{

	this->seed = seed;
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	cc_cus = new cc_counter_t*[height];
	for (int i = 0; i < height; ++i)
	{
		cc_cus[i] = new cc_counter_t[width]();
	}

	bobhash = new BOBHash[height];
	for (int i = 0; i < height; ++i)
	{
		bobhash[i].initialize(i + 1000);
	}

	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");

	width_over_four = width >> 2;

	indices = new int[height];
}

ConservativeUpdateCC::~ConservativeUpdateCC()
{
	delete[] indices;
	delete[] bobhash;
	
	for (int i = 0; i < height; ++i)
	{
		delete[] cc_cus[i];
	}
	delete[] cc_cus;	
}

void ConservativeUpdateCC::divide_array_counters()
{
	for (int i = 0; i < height; ++i)
	{
		uint64_t* four_counters_array = (uint64_t*)cc_cus[i];
		for (unsigned j = 0; j < width_over_four; ++j)
		{
			uint64_t &four_counters = four_counters_array[j];
			four_counters = (four_counters & division_mask) >> 1;
		}
	}
}

void ConservativeUpdateCC::increment(const char * str)
{
	if (minus_log_p == 0)
	{
		increment_h(str);
	}
	else if (minus_log_p < minus_log_p_bound)
	{
		if (RandByteArray[rand_pos >> 3] & (1 << ((rand_pos & 0x7) - 1))) // TODO: lookup table 
		{
			increment_h(str);
		}
		else
		{
		}

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

void ConservativeUpdateCC::increment_h(const char * str)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	indices[0] = bobhash[0].run(str, FT_SIZE) % width;
	cc_counter_t min = cc_cus[0][indices[0]];
	for (int i = 1; i < height; ++i)
	{
		indices[i] = bobhash[i].run(str, FT_SIZE) % width;
		if (min > cc_cus[i][indices[i]])
		{
			min = cc_cus[i][indices[i]];
		}
	}

	if (min == 65535)
	{
		++minus_log_p;
		LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
		divide_array_counters();
		min = 32767;
	}

	for (int i = 0; i < height; ++i)
	{
		cc_counter_t &comp_counter = cc_cus[i][indices[i]];
		if (comp_counter == min)
		{
			++comp_counter;
		}		
	}
}

uint64_t ConservativeUpdateCC::query(const char * str)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	cc_counter_t min = cc_cus[0][bobhash[0].run(str, FT_SIZE) % width];
	for (int i = 1; i < height; ++i)
	{
		cc_counter_t &comp_counter = cc_cus[i][bobhash[i].run(str, FT_SIZE) % width];
		if (min > comp_counter)
		{
			min = comp_counter;
		}
	}
	return ((uint64_t)min) << minus_log_p;
}

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

inline unsigned CountMinCC_MaxSpeed::accelerated_geometric(double LogOneMinusP) {
	// Return a randomly drawn geometric variable with parameter p, such that LogOneMinusP = log(1-p)
	return log(dis_arr(generator_arr)) / LogOneMinusP + 1;
}

CountMinCC_MaxSpeed::CountMinCC_MaxSpeed(int width, int height, int seed, double epsilon_cc, double delta_cc) :

	gen_arr(seed),
	dis_arr(0, 1),
	rand_bits(_mm256_set_epi64x(UINT64_MAX, UINT64_MAX, UINT64_MAX, UINT64_MAX)),
	RandByteArray((char *)&rand_bits),
	weighted_rand_bits(_mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr())),
	WeightedRandByteArray((uint16_t *)&weighted_rand_bits)
{

	this->seed = seed;
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	cc_cms = new cc_counter_t*[height];
	for (int i = 0; i < height; ++i)
	{
		cc_cms[i] = new cc_counter_t[width]();
	}

	bobhash = new BOBHash[height];
	for (int i = 0; i < height; ++i)
	{
		bobhash[i].initialize(i + 1000);
	}

	generator_arr.seed(seed);
	LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p_bound));
	interval = accelerated_geometric(LogOneMinusP);

	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");

	width_over_four = width >> 2;

	indices = new int[height];

	N_current = 0;
	N_prime = 2 * (1 + epsilon_cc / 3)*(1 / (epsilon_cc*epsilon_cc))*log(2 / delta_cc) + 1;
	Shifted_N_prime = N_prime;

	assert((N_prime <= 1 << 16) && "N_prime is too large! For this (eplison, delta) we need more than 16 bits!");
}

CountMinCC_MaxSpeed::~CountMinCC_MaxSpeed()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] cc_cms[i];
	}

	delete[] cc_cms;
	delete[] bobhash;
	delete[] indices;
}

void CountMinCC_MaxSpeed::divide_array_counters()
{
	for (int i = 0; i < height; ++i)
	{
		uint64_t* four_counters_array = (uint64_t*)cc_cms[i];
		for (unsigned j = 0; j < width_over_four; ++j)
		{
			uint64_t &four_counters = four_counters_array[j];
			four_counters = (four_counters & division_mask) >> 1;
		}
	}
}

void CountMinCC_MaxSpeed::increment(const char * str)
{
	++N_current;

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

void CountMinCC_MaxSpeed::increment_h(const char * str)
{

	if (N_current > Shifted_N_prime)
	{
		Shifted_N_prime *= 2;
		++minus_log_p;
		LogOneMinusP = log(1 - 1.0 / (1 << minus_log_p));
		divide_array_counters();

		// if p decreases we need to reduce sampling rate accordingly.
		rand_bits = _mm256_and_si256(rand_bits, _mm256_set_epi64x(gen_arr(), gen_arr(), gen_arr(), gen_arr()));
	}

	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	for (int i = 0; i < height; ++i)
	{
		cc_counter_t &comp_counter = cc_cms[i][bobhash[i].run(str, FT_SIZE) % width];
		++comp_counter;

		assert((comp_counter <= 65535) && "This cannot happen! MaxSpeed is Bob!");
	}
}

uint64_t CountMinCC_MaxSpeed::query(const char * str)
{
	// approximation of the correct logic - as we may divide after going only over a fraction of the counters.
	cc_counter_t min = cc_cms[0][bobhash[0].run(str, FT_SIZE) % width];
	for (int i = 1; i < height; ++i)
	{
		cc_counter_t &comp_counter = cc_cms[i][bobhash[i].run(str, FT_SIZE) % width];
		if (min > comp_counter)
		{
			min = comp_counter;
		}
	}
	return ((uint64_t)min) << minus_log_p;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
