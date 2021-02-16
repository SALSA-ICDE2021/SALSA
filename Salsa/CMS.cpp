#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string.h>

#include "CMS.hpp"

using namespace std;

///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////

CountMinBaseline::CountMinBaseline()
{
}

CountMinBaseline::~CountMinBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms[i];
	}
	delete[] bobhash;
	delete[] baseline_cms;
}

void CountMinBaseline::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cms = new uint32_t*[height];
	bobhash = new BOBHash[height];

	for (int i = 0; i < height; ++i)
	{
		baseline_cms[i] = new uint32_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}
}

void CountMinBaseline::increment(const char * str)
{
	for (int i = 0; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		++baseline_cms[i][index];
	}
}

uint64_t CountMinBaseline::query(const char * str)
{
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
	uint64_t min = baseline_cms[0][index];
	for (int i = 1; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = baseline_cms[i][index];
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

CountMinBaselineSanity::CountMinBaselineSanity()
{
}

CountMinBaselineSanity::~CountMinBaselineSanity()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms[i];
	}
	delete[] bobhash;
	delete[] baseline_cms;
}

void CountMinBaselineSanity::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = (width - 1) << 5;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cms = new uint32_t*[height];
	bobhash = new BOBHash[height];

	for (int i = 0; i < height; ++i)
	{
		baseline_cms[i] = new uint32_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}
}

void CountMinBaselineSanity::increment(const char * str)
{
	for (int i = 0; i < height; ++i) {
		uint index = ((bobhash[i].run(str, FT_SIZE)) & width_mask) >> 5;
		++baseline_cms[i][index];
	}
}

uint64_t CountMinBaselineSanity::query(const char * str)
{
	uint index = ((bobhash[0].run(str, FT_SIZE)) & width_mask) >> 5;
	uint64_t min = baseline_cms[0][index];
	for (int i = 1; i < height; ++i) {
		uint index = ((bobhash[i].run(str, FT_SIZE)) & width_mask) >> 5;
		uint64_t temp = baseline_cms[i][index];
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

WeightedCountMinBaseline::WeightedCountMinBaseline()
{
}

WeightedCountMinBaseline::~WeightedCountMinBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cms[i];
	}
	delete[] bobhash;
	delete[] baseline_cms;
}

void WeightedCountMinBaseline::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cms = new uint64_t*[height];
	bobhash = new BOBHash[height];

	for (int i = 0; i < height; ++i)
	{
		baseline_cms[i] = new uint64_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}
}

void WeightedCountMinBaseline::add(const char * str, int c)
{
	for (int i = 0; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		baseline_cms[i][index] += c;
	}
}

uint64_t WeightedCountMinBaseline::query(const char * str)
{
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
	uint64_t min = baseline_cms[0][index];
	for (int i = 1; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = baseline_cms[i][index];
		if (min > temp)
		{
			min = temp;
		}
	}
	return min;
}

/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////

ConservativeUpdateBaseline::ConservativeUpdateBaseline()
{
}

ConservativeUpdateBaseline::~ConservativeUpdateBaseline()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cu[i];
	}
	delete[] bobhash;
	delete[] baseline_cu;

	delete[] counters;
}

void ConservativeUpdateBaseline::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = width - 1;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cu = new uint32_t*[height];
	bobhash = new BOBHash[height];

	for (int i = 0; i < height; ++i)
	{
		baseline_cu[i] = new uint32_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}

	counters = new uint32_t*[height];
}

void ConservativeUpdateBaseline::increment(const char * str)
{
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
	counters[0] = baseline_cu[0] + index;
	uint64_t min = *counters[0];
	
	for (int i = 1; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		counters[i] = baseline_cu[i] + index;
		uint64_t temp = *counters[i];
		if (min > temp)
		{
			min = temp;
		}
	}

	for (int i = 0; i < height; ++i) {
		if (*counters[i] == min)
		{
			++(*counters[i]);
		}	
	}
}

void ConservativeUpdateBaseline::insert(const char * str, int weight)
{
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
	counters[0] = baseline_cu[0] + index;
	uint64_t min = *counters[0];

	for (int i = 1; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		counters[i] = baseline_cu[i] + index;
		uint64_t temp = *counters[i];
		if (min > temp)
		{
			min = temp;
		}
	}

	uint64_t new_min = min + weight;

	for (int i = 0; i < height; ++i) {
		if (*counters[i] < new_min)
		{
			*counters[i] = new_min;
		}
	}
}

uint64_t ConservativeUpdateBaseline::query(const char * str)
{
	uint index = (bobhash[0].run(str, FT_SIZE)) & width_mask;
	uint64_t min = baseline_cu[0][index];
	for (int i = 1; i < height; ++i) {
		uint index = (bobhash[i].run(str, FT_SIZE)) & width_mask;
		uint64_t temp = baseline_cu[i][index];
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

ConservativeUpdateBaselineSanity::ConservativeUpdateBaselineSanity()
{
}

ConservativeUpdateBaselineSanity::~ConservativeUpdateBaselineSanity()
{
	for (int i = 0; i < height; ++i)
	{
		delete[] baseline_cu[i];
	}
	delete[] bobhash;
	delete[] baseline_cu;

	delete[] counters;
}

void ConservativeUpdateBaselineSanity::initialize(int width, int height, int seed)
{
	this->width = width;
	this->height = height;

	width_mask = (width - 1) << 5;

	assert(width > 0 && "We assume too much!");
	assert(width % 4 == 0 && "We assume that (w % 4 == 0)!");
	assert((width & (width - 1)) == 0 && "We assume that width is a power of 2!");

	baseline_cu = new uint32_t*[height];
	bobhash = new BOBHash[height];

	for (int i = 0; i < height; ++i)
	{
		baseline_cu[i] = new uint32_t[width]();
		bobhash[i].initialize(seed*(7 + i) + i + 100);
	}

	counters = new uint32_t*[height];
}

void ConservativeUpdateBaselineSanity::increment(const char * str)
{ 
	uint index = ((bobhash[0].run(str, FT_SIZE)) & width_mask) >> 5;
	counters[0] = baseline_cu[0] + index;
	uint64_t min = *counters[0];

	for (int i = 1; i < height; ++i) {
		uint index = ((bobhash[i].run(str, FT_SIZE)) & width_mask) >> 5;
		counters[i] = baseline_cu[i] + index;
		uint64_t temp = *counters[i];
		if (min > temp)
		{
			min = temp;
		}
	}

	for (int i = 0; i < height; ++i) {
		if (*counters[i] == min)
		{
			++(*counters[i]);
		}
	}
}

uint64_t ConservativeUpdateBaselineSanity::query(const char * str)
{
	uint index = ((bobhash[0].run(str, FT_SIZE)) & width_mask) >> 5;
	uint64_t min = baseline_cu[0][index];
	for (int i = 1; i < height; ++i) {
		uint index = ((bobhash[i].run(str, FT_SIZE)) & width_mask) >> 5;
		uint64_t temp = baseline_cu[i][index];
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
