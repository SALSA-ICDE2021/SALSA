#ifndef UnivMon_TESTS
#define UnivMon_TESTS

#include "UnivMon.hpp"

void test_univmon_final_count_distinct(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge);
void test_univmon_final_entropy(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge);
void test_univmon_final_Fp(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge, double p);
void test_univmon_final_all(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge);

#endif