#include "UnivMon.hpp"

UnivSketchFgSalsa::UnivSketchFgSalsa(int counterSize, int width, int height, int seed, int downsamplings_per_merge, int HH_num, bool maximum_merge) :
	sketches{
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed, downsamplings_per_merge, HH_num, maximum_merge), // 1
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 1, downsamplings_per_merge, HH_num, maximum_merge), // 2
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 2, downsamplings_per_merge, HH_num, maximum_merge), // 3
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 3, downsamplings_per_merge, HH_num, maximum_merge), // 4
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 4, downsamplings_per_merge, HH_num, maximum_merge), // 5
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 5, downsamplings_per_merge, HH_num, maximum_merge), // 6
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 6, downsamplings_per_merge, HH_num, maximum_merge), // 7
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 7, downsamplings_per_merge, HH_num, maximum_merge), // 8
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 8, downsamplings_per_merge, HH_num, maximum_merge), // 9
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 9, downsamplings_per_merge, HH_num, maximum_merge), // 10
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 10, downsamplings_per_merge, HH_num, maximum_merge), // 11
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 11, downsamplings_per_merge, HH_num, maximum_merge), // 12
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 12, downsamplings_per_merge, HH_num, maximum_merge), // 13
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 13, downsamplings_per_merge, HH_num, maximum_merge), // 14
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 14, downsamplings_per_merge, HH_num, maximum_merge), // 15
		FineGrainedSalsaCountSketchHH(counterSize, width, height, seed + 15, downsamplings_per_merge, HH_num, maximum_merge)  // 16
}
{
	LastLevelBob.initialize(seed + 15);
	for (int lvl = 0; lvl < UnivMon_CS_LVLS; lvl++)
	{
		BobArray[lvl].initialize(seed + 42 + 10 * lvl);
	}
	num_increments = 0;
}

UnivSketchFgSalsa::~UnivSketchFgSalsa()
{
}

void UnivSketchFgSalsa::increment(const char * str) {

	uint16_t levelHash = LastLevelBob.run(str, FT_SIZE) & 0xFFFF;
	int lastLevel = __lzcnt16(levelHash | 0b1);
	//cout << "lastLevel" << lastLevel << endl;
	for (int lvl = 0; lvl <= lastLevel; lvl++)
	{
		sketches[lvl].increment(str);
	}
	++num_increments;

}

long double UnivSketchFgSalsa::query_count_distint() {

	long double y_bottom = (long double)sketches[UnivMon_CS_LVLS - 1].HH().size();
	long double y_1 = 0.0, y_2 = 0.0;

	//cout << y_bottom << endl;

	y_2 = y_bottom;

	for (int lvl = UnivMon_CS_LVLS - 2; lvl >= 0; lvl--)
	{
		long double indSum = 0.0;
		vector<pair<string, uint64_t>> HH_lvl = sketches[lvl].HH();

		for (int j = 0; j < HH_lvl.size(); j++)
		{

			uint16_t levelHash = LastLevelBob.run(HH_lvl[j].first.c_str(), FT_SIZE) & 0xFFFF;
			int lastLevel = __lzcnt16(levelHash | 0b1);

			long double hash = (lastLevel == lvl) ? 0 : 1;

			indSum += (1.0 - 2.0 * hash);
		}
		
		y_1 = 2.0 * y_2 + indSum;
		y_2 = y_1;
		//cout << y_1 << " " << lvl << " " << HH_lvl.size() << endl;
	}
	return y_1;

}

long double UnivSketchFgSalsa::query_Fp(double p)
{

	assert(p <= 2 && "We assume too much!");
	assert(p >= 0 && "We assume too much!");

	long double y_bottom = 0.0;
	vector<pair<string, uint64_t>> HH_lvl = sketches[UnivMon_CS_LVLS - 1].HH();
	for (int i = 0; i < HH_lvl.size(); ++i)
	{
		long double w = (long double)HH_lvl[i].second;
		w = w > 1.0 ? w : 1.0;
		y_bottom += powl(w, p);
	}

	long double y_1 = 0.0, y_2 = 0.0;

	y_2 = y_bottom;

	for (int lvl = UnivMon_CS_LVLS - 2; lvl >= 0; lvl--)
	{
		long double indSum = 0.0;
		vector<pair<string, uint64_t>> HH_lvl = sketches[lvl].HH();

		for (int j = 0; j < HH_lvl.size(); j++)
		{
			long double w = (long double)HH_lvl[j].second;
			w = w > 1.0 ? w : 1.0;
			uint16_t levelHash = LastLevelBob.run(HH_lvl[j].first.c_str(), FT_SIZE) & 0xFFFF;
			int lastLevel = __lzcnt16(levelHash | 0b1);
			long double hash = (lastLevel == lvl) ? 0 : 1;
			indSum += (1.0 - 2.0 * hash)*powl(w, p);
		}
		
		y_1 = 2.0 * y_2 + indSum;
		y_2 = y_1;
	}
	return y_1;

}

long double UnivSketchFgSalsa::query_entropy()
{
	long double y_bottom = 0.0;
	vector<pair<string, uint64_t>> HH_lvl = sketches[UnivMon_CS_LVLS - 1].HH();
	for (int i = 0; i < HH_lvl.size(); ++i)
	{
		long double w = (long double)HH_lvl[i].second;
		if (w > 0.0)
		{
			y_bottom += w * log2l(w);
		}
	}

	long double y_1 = 0.0, y_2 = 0.0;

	y_2 = y_bottom;

	for (int lvl = UnivMon_CS_LVLS - 2; lvl >= 0; lvl--)
	{
		long double indSum = 0.0;
		vector<pair<string, uint64_t>> HH_lvl = sketches[lvl].HH();

		for (int j = 0; j < HH_lvl.size(); j++)
		{
			long double w = (long double)HH_lvl[j].second;
			if (w > 0.0)
			{
				uint16_t levelHash = LastLevelBob.run(HH_lvl[j].first.c_str(), FT_SIZE) & 0xFFFF;
				int lastLevel = __lzcnt16(levelHash | 0b1);
				long double hash = (lastLevel == lvl) ? 0 : 1;

				indSum += (1.0 - 2.0 * hash)*w*log2l(w);
			}
		}
		y_1 = 2.0 * y_2 + indSum;
		y_2 = y_1;
	}
	long double entropy = log2l((long double)num_increments) - (y_1 / (long double)num_increments);
	return entropy;
}

