#ifndef UnivMon_H
#define UnivMon_H

#include "Defs.hpp"
#include "CountSketch.hpp"
#include "BobHash.hpp"

class UnivSketchFgSalsa 
{
private:
	FineGrainedSalsaCountSketchHH sketches[UnivMon_CS_LVLS];
	BOBHash LastLevelBob;
	BOBHash BobArray[UnivMon_CS_LVLS];

	long double num_increments;

public:
	UnivSketchFgSalsa(int counterSize, int width, int height, int seed, int downsamplings_per_merge, int HH_num, bool maximum_merge);
	~UnivSketchFgSalsa();

	void increment(const char * str);

	long double query_count_distint();
	long double query_Fp(double p);
	long double query_entropy();
};

#endif
