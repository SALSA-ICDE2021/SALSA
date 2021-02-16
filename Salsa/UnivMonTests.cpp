#include <string>
#include "UnivMonTests.hpp"

void test_univmon_final_count_distinct(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge)
{

	UnivSketchFgSalsa baseline(32, width, height, seed, downsamplings_per_merge, HH_num, maximum_merge);
	UnivSketchFgSalsa not_baseline(salsa_bits, width*(32 / salsa_bits), height, seed, downsamplings_per_merge, HH_num, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	int mycounter = 0;
	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		if (++mycounter % 100000 == 0)
		{
			cout << mycounter << endl;
		}

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		true_sizes[ft_key] += 1;
		baseline.increment(data + i);
		not_baseline.increment(data + i);
	}


	cout << "True Count Distinct:" << true_sizes.size() << endl;
	cout << "Baseline Count Distinct:" << baseline.query_count_distint() << endl;
	cout << "Salsa Count Distinct:" << not_baseline.query_count_distint() << endl;

}

void test_univmon_final_entropy(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge)
{

	UnivSketchFgSalsa baseline(32, width, height, seed, downsamplings_per_merge, HH_num, maximum_merge);
	UnivSketchFgSalsa not_baseline(salsa_bits, width*(32 / salsa_bits), height, seed, downsamplings_per_merge, HH_num, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	int mycounter = 0;
	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		if (++mycounter % 100000 == 0)
		{
			cout << mycounter << endl;
		}

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		true_sizes[ft_key] += 1;
		baseline.increment(data + i);
		not_baseline.increment(data + i);
	}

	long double true_entopy = 0;
	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		true_entopy += it->second*log2l(it->second);
	}
	true_entopy = log2l((long double)N) - (true_entopy / (long double)N);

	cout << "True entropy:" << true_entopy << endl;
	cout << "Baseline entropy:" << baseline.query_entropy() << endl;
	cout << "Salsa entropy:" << not_baseline.query_entropy() << endl;

}

void test_univmon_final_Fp(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge, double p)
{

	UnivSketchFgSalsa baseline(32, width, height, seed, downsamplings_per_merge, HH_num, maximum_merge);
	UnivSketchFgSalsa not_baseline(salsa_bits, width*(32 / salsa_bits), height, seed, downsamplings_per_merge, HH_num, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	int mycounter = 0;
	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		if (++mycounter % 100000 == 0)
		{
			cout << mycounter << endl;
		}

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		true_sizes[ft_key] += 1;
		baseline.increment(data + i);
		not_baseline.increment(data + i);
	}

	long double Fp = 0;
	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		Fp += powl(it->second, p);
	}

	cout << "True Fp:" << Fp << endl;
	cout << "Baseline Fp:" << baseline.query_Fp(p) << endl;
	cout << "Salsa Fp:" << not_baseline.query_Fp(p) << endl;

}

void test_univmon_final_all(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge)
{

	UnivSketchFgSalsa not_baseline(salsa_bits, width*(32 / salsa_bits), height, seed, downsamplings_per_merge, HH_num, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	int mycounter = 0;
	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		//if (++mycounter % 100000 == 0)
		//{
		//	cout << mycounter << endl;
		//}

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		true_sizes[ft_key] += 1;
		not_baseline.increment(data + i);
	}

	long double true_entopy = 0;
	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		true_entopy += it->second*log2l(it->second);
	}
	true_entopy = log2l((long double)N) - (true_entopy / (long double)N);

	long double true_count_distinct = true_sizes.size();

	long double true_Fp[201] = { 0 };
	for (int i = 0; i <= 200; ++i)
	{
		double p = (double)i / 100;
		for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
		{
			true_Fp[i] += powl(it->second, p);
		}
	}

	long double univ_entopy = not_baseline.query_entropy();
	long double univ_count_distinct = not_baseline.query_count_distint();

	long double univ_Fp[201] = { 0 };
	for (int i = 0; i <= 200; ++i)
	{
		double p = (double)i / 100;
		univ_Fp[i] = not_baseline.query_Fp(p);
	}

	ofstream results_file_error_on_arrival;
	string file_name_error_on_arrival = "test_fg_salsa_univmon_fp_";
	file_name_error_on_arrival.append(to_string(salsa_bits));
	file_name_error_on_arrival.append("_seed_");
	file_name_error_on_arrival.append(to_string(seed));
	file_name_error_on_arrival.append("_downsamplings_per_merge_");
	file_name_error_on_arrival.append(to_string(downsamplings_per_merge));
	file_name_error_on_arrival.append(".txt");
	results_file_error_on_arrival.open(file_name_error_on_arrival, ofstream::out | ofstream::app);
	results_file_error_on_arrival << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height;

	for (int i = 0; i <= 200; ++i)
	{
		double p = (double)i / 100;
		results_file_error_on_arrival << "\tp\t" << p << "\ttruefp\t" << true_Fp[i] << "\tsalsafp\t" << univ_Fp[i];
	}
	results_file_error_on_arrival << endl;


	//cout << "True count distinct:" << true_count_distinct << endl;
	//cout << "Salsa " << salsa_bits << " bits count distinct:" << univ_count_distinct << endl;

	ofstream results_file_entropy;
	string file_name_entropy = "test_fg_salsa_univmon_entropy_";
	file_name_entropy.append(to_string(salsa_bits));
	file_name_entropy.append("_seed_");
	file_name_entropy.append(to_string(seed));
	file_name_entropy.append("_downsamplings_per_merge_");
	file_name_entropy.append(to_string(downsamplings_per_merge));
	file_name_entropy.append(".txt");
	results_file_entropy.open(file_name_entropy, ofstream::out | ofstream::app);
	results_file_entropy << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\ttrueentropy\t" << true_entopy << "\tsalsaentropy\t" << univ_entopy << endl;

}
