#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <string>
#include <unordered_map>

#include "CountSketch.hpp"
#include "CountSketchTests.hpp"

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_baseline_count_sketch_error_on_arrival(int N, int width, int height, int seed, const char* data)
{

	CountSketchBaseline ds(width, height, seed);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		long double point_error = (int64_t)ds.query(data + i) - (int64_t)true_sizes[ft_key];

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		ds.increment(data + i);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_baseline_count_sketch_error_on_arrival_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_baseline_count_sketch_speed(int N, int width, int height, int seed, const char* data)
{

	CountSketchBaseline ds(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		ds.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_baseline_count_sketch_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	results_file.open("test_baseline_count_sketch_speed.txt", ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_baseline_count_sketch_five_rows_error_on_arrival(int N, int width, int height, int seed, const char* data)
{
	assert(height == 5 && "We assume that height is 5");

	CountSketchBaselineFiveRows ds(width, seed);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		long double point_error = (int64_t)ds.query(data + i) - (int64_t)true_sizes[ft_key];

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		ds.increment(data + i);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_baseline_count_sketch_five_rows_error_on_arrival_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_baseline_count_sketch_five_rows_speed(int N, int width, int height, int seed, const char* data)
{
	assert(height == 5 && "We assume that height is 5");

	CountSketchBaselineFiveRows ds(width, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		ds.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_baseline_count_five_rows_sketch_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	results_file.open("test_baseline_count_five_rows_sketch_speed.txt", ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_fg_salsa_count_sketch_error_on_arrival(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, bool maximum_merge)
{

	FineGrainedSalsaCountSketch cms_baseline(salsa_bits, width, height, seed, downsamplings_per_merge, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		long double point_error = (int64_t)cms_baseline.query(data + i) - (int64_t)true_sizes[ft_key];

		if (false)

		{
			int ggg = 4953;
			//cout << "we are not BOB! BOB is guilty! " << 559 << endl;
			//cout << "we are not BOB! BOB is guilty! " << ft_key << endl;
			uint64_t ftkey = (((uint64_t)ft_to_bobkey_1.run(data + ggg, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + ggg, FT_SIZE);
			cout << (int64_t)cms_baseline.query(data + ggg) << " ";
			cout << (int64_t)true_sizes[ftkey] << " " << i << " ";

			cout << cms_baseline.read_row_value_ft(data + ggg, 0) << " ";
			cout << cms_baseline.read_row_value_ft(data + ggg, 1) << " ";
			cout << cms_baseline.read_row_value_ft(data + ggg, 2) << " ";
			cout << cms_baseline.read_row_value_ft(data + ggg, 3) << " " << point_error << endl;
		}

		//cout << point_error << " " << ft_key << " " << (int64_t)cms_baseline.query(data + i) << " " << (int64_t)true_sizes[ft_key] << endl;

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cms_baseline.increment(data + i);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string file_name = "no_max_test_fg_salsa_count_sketch_error_on_arrival_tb_";
	file_name.append(to_string(salsa_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append("_downsamplings_per_merge_");
	file_name.append(to_string(downsamplings_per_merge));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_fg_salsa_split_counters_sanity_count_sketch(int N, int width, int height, int seed, const char* data, bool maximum_merge)
{

	CountSketchBaseline baseline(width, height, seed);
	FineGrainedSalsaCountSketch fg_salsa(32, width, height, seed, 0, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	/*
	int last_4216 = 0;
	int last_4217 = 0;
	int last_4218 = 0;
	int last_4219 = 0;

	int last_2108 = 0;
	int last_2109 = 0;
	*/

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		long double point_error_1 = (int64_t)baseline.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_2 = (int64_t)fg_salsa.query(data + i) - (int64_t)true_sizes[ft_key];

		if (abs(point_error_2) != abs(point_error_1))

		{
			cout << "we are not BOB! BOB is guilty! " << i << endl;
			cout << "we are not BOB! BOB is guilty! " << ft_key << endl;
			cout << (int64_t)baseline.query(data + i) << endl;
			cout << (int64_t)fg_salsa.query(data + i) << endl;
			cout << (int64_t)true_sizes[ft_key];
			system("pause");
			exit(1);
		}

		//cout << i << "\t" << tango.read_row_value(4105, 0)
		//	<< "\t" << tango.read_row_value(5937, 1)
		//	<< "\t" << tango.read_row_value(4108, 2)
		//	<< "\t" << tango.read_row_value(5119, 3) << endl;

		//if (tango.query(data + 1137630) > maximum_salsa_cms_baseline_sanity.query(data + 1137630))
		//{
		//	cout << i << "\t" << tango.query(data + 1137630) << endl; // 
		//	system("pause");
		//}
		/*
		int query_i = 4489225;

		if ((last_4216 != tango_2_bit.read_row_value(4216, 0)) ||
			(last_4217 != tango_2_bit.read_row_value(4217, 0)) ||
			(last_4218 != tango_2_bit.read_row_value(4218, 0)) ||
			(last_4219 != tango_2_bit.read_row_value(4219, 0)) ||
			(last_2108 != tango_4_bit.read_row_value(2108, 0)) ||
			(last_2109 != tango_4_bit.read_row_value(2109, 0)))//(tango_2_bit.read_row_value_ft(data + query_i, 0) > tango_4_bit.read_row_value_ft(data + query_i, 0))
		{
			cout << i << "\t" << 0 << "\t 2bit \t" << tango_2_bit.read_row_value(4216, 0) << " "
												   << tango_2_bit.read_row_value(4217, 0) << " "
												   << tango_2_bit.read_row_value(4218, 0) << " "
												   << tango_2_bit.read_row_value(4219, 0) << " "
																						  << endl;

			cout << i << "\t" << 0 << "\t 4bit \t" << tango_4_bit.read_row_value(2108, 0) << " "
												   << tango_4_bit.read_row_value(2109, 0) << " "
												   << endl;

			last_4216 = tango_2_bit.read_row_value(4216, 0);
			last_4217 = tango_2_bit.read_row_value(4217, 0);
			last_4218 = tango_2_bit.read_row_value(4218, 0);
			last_4219 = tango_2_bit.read_row_value(4219, 0);

			last_2108 = tango_4_bit.read_row_value(2108, 0);
			last_2109 = tango_4_bit.read_row_value(2109, 0);

			//cout << i << "\t" << 0 << "\t 2bit \t" << tango_2_bit.read_row_value_ft(data + query_i, 0) << endl; //
			//cout << i << "\t" << 0 << "\t 4bit \t" << tango_4_bit.read_row_value_ft(data + query_i, 0) << endl; //
			//system("pause");
		}
		*/
		/*
		int p = 546;
		
		if (i >= 0)
		{
			

			cout << "baseline: " << baseline.read_row_value_ft(data + p, 0) << " ";
			cout << baseline.read_row_value_ft(data + p, 1) << " ";
			cout << baseline.read_row_value_ft(data + p, 2) << " ";
			cout << baseline.read_row_value_ft(data + p, 3) << endl;

			cout << "salsa fg: " << fg_salsa.read_row_value_ft(data + p, 0) << " ";
			cout << fg_salsa.read_row_value_ft(data + p, 1) << " ";
			cout << fg_salsa.read_row_value_ft(data + p, 2) << " ";
			cout << fg_salsa.read_row_value_ft(data + p, 3) << " ";
			cout << (int)fg_salsa.log_counter_size(58, 3) << "\t\t" << i << endl;
		}
		*/

		/*
		if (tango.read_row_value_ft(data + query_i, 1) > maximum_salsa_cms_baseline_sanity.read_row_value_ft(data + query_i, 1))
		{
			cout << i << "\t" << 1 << "\t" << tango.read_row_value_ft(data + query_i, 1) << endl; //
			cout << i << "\t" << 1 << "\t" << maximum_salsa_cms_baseline_sanity.read_row_value_ft(data + query_i, 1) << endl; //
			system("pause");
		}
		if (tango.read_row_value_ft(data + query_i, 2) > maximum_salsa_cms_baseline_sanity.read_row_value_ft(data + query_i, 2))
		{
			cout << i << "\t" << 2 << "\t" << tango.read_row_value_ft(data + query_i, 2) << endl; //
			cout << i << "\t" << 2 << "\t" << maximum_salsa_cms_baseline_sanity.read_row_value_ft(data + query_i, 2) << endl; //
			system("pause");
		}
		if (tango.read_row_value_ft(data + query_i, 3) > maximum_salsa_cms_baseline_sanity.read_row_value_ft(data + query_i, 3))
		{
			cout << i << "\t" << 3 << "\t" << tango.read_row_value_ft(data + query_i, 3) << endl; //
			cout << i << "\t" << 3 << "\t" << maximum_salsa_cms_baseline_sanity.read_row_value_ft(data + query_i, 3) << endl; //
			system("pause");
		}
		*/

		//if (i == 4473313)
		//{
		//	cout << i << endl;
		//}

		true_sizes[ft_key] += 1;
		baseline.increment(data + i);
		fg_salsa.increment(data + i);

	}
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_fg_salsa_count_sketch_hh(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge)
{

	FineGrainedSalsaCountSketchHH baseline(32, width, height, seed, downsamplings_per_merge, HH_num, maximum_merge);

	FineGrainedSalsaCountSketchHH fg_salsa_hh(salsa_bits, width*(32 / salsa_bits), height, seed, downsamplings_per_merge, HH_num, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		true_sizes[ft_key] += 1;
		baseline.increment(data + i);
		fg_salsa_hh.increment(data + i);
	}

	vector<pair<uint64_t, uint64_t>> true_top_k(HH_num);
	partial_sort_copy(true_sizes.begin(),
		true_sizes.end(),
		true_top_k.begin(),
		true_top_k.end(),
		[](pair<uint64_t, uint64_t> const& l,
			pair<uint64_t, uint64_t> const& r)
	{
		return l.second > r.second;
	});

	uint64_t k_largest = true_top_k[0].second;
	for (int i = 1; i < true_top_k.size(); i++)
	{
		k_largest = (k_largest > true_top_k[i].second) ? true_top_k[i].second : k_largest;
	}

	vector<pair<string, uint64_t>> baseline_top_k = baseline.HH();
	vector<pair<string, uint64_t>> fg_salsa_hh_top_k = fg_salsa_hh.HH();

	int baseline_hits = 0;
	int salsa_fg_hits = 0;

	for (int i = 0; i < HH_num; ++i)
	{
		uint64_t ft_key_baseline = (((uint64_t)ft_to_bobkey_1.run(baseline_top_k[i].first.c_str(), FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(baseline_top_k[i].first.c_str(), FT_SIZE);
		uint64_t ft_key_salsa_fg = (((uint64_t)ft_to_bobkey_1.run(fg_salsa_hh_top_k[i].first.c_str(), FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(fg_salsa_hh_top_k[i].first.c_str(), FT_SIZE);

		baseline_hits += (true_sizes[ft_key_baseline] >= k_largest);
		salsa_fg_hits += (true_sizes[ft_key_salsa_fg] >= k_largest);
	}

	cout << "Baseline hits:" << baseline_hits << " out of " << HH_num << endl;
	cout << "Sasla fg hits:" << salsa_fg_hits << " out of " << HH_num << endl;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_fg_salsa_count_sketch_l2(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, bool maximum_merge)
{

	FineGrainedSalsaCountSketch fg_salsa_hh(salsa_bits, width*(32 / salsa_bits), height, seed, downsamplings_per_merge, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		true_sizes[ft_key] += 1;
		fg_salsa_hh.increment(data + i);
	}

	long double true_l2 = 0;
	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		true_l2 += it->second*it->second;
	}
	true_l2 = sqrtl(true_l2);

	cout << "True l2:" << true_l2 << endl;
	cout << "Sasla fg l2:" << fg_salsa_hh.l2_estimation() << endl;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_fg_salsa_count_sketch_turnstile(int salsa_bits, int N, int width, int height, int seed, const char* data)
{
	// debug!!! 8 and 4 bit not good!

	FineGrainedSalsaCountSketch cms_baseline(salsa_bits, width, height, seed, 0, false);
	FineGrainedSalsaCountSketch cms_baseline2(salsa_bits, width, height, seed, 0, false);

	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> true_sizes2;

	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop / 2; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		ft_key_i_values[ft_key] = i;

		true_sizes[ft_key] += 1;
		cms_baseline.increment(data + i);
	}

	for (int64_t i = stop_loop / 2; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		ft_key_i_values[ft_key] = i;

		true_sizes2[ft_key] += 1;
		cms_baseline2.increment(data + i);
	}

	FineGrainedSalsaCountSketch* diffSketch = cms_baseline - cms_baseline2;

	uint64_t SOS = 0;
	for (auto it = ft_key_i_values.begin(); it != ft_key_i_values.end(); it++)
	{
		uint64_t current_ft_key = it->first;
		int64_t current_index = it->second;

		int64_t true_sizes_val;
		int64_t true_sizes_val2;

		if (true_sizes.find(current_ft_key) == true_sizes.end()) {
			// not found
			true_sizes_val = 0;
		}
		else {
			// found
			true_sizes_val = true_sizes[current_ft_key];
		}

		if (true_sizes2.find(current_ft_key) == true_sizes2.end()) {
			// not found
			true_sizes_val2 = 0;
		}
		else {
			// found
			true_sizes_val2 = true_sizes2[current_ft_key];
		}

		int64_t true_change = true_sizes_val - true_sizes_val2;
		int64_t estimated_change = diffSketch->signedQuery(data + current_index);

		int64_t error = true_change - estimated_change;

		SOS += error * error;

	}

	delete diffSketch;

	long double mse = (long double)SOS / (long double)ft_key_i_values.size();
	long double rmse = sqrtl(mse);

	ofstream results_file;
	string file_name = "test_fg_salsa_count_sketch_turnstile_";
	file_name.append(to_string(salsa_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tRMSE\t" << rmse << endl;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_fg_salsa_count_sketch_all_except_turnstile(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge, int HH_num, bool maximum_merge)
{

	FineGrainedSalsaCountSketchHH cms_baseline(salsa_bits, width, height, seed, downsamplings_per_merge, HH_num, maximum_merge);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		long double point_error = (int64_t)cms_baseline.query(data + i) - (int64_t)true_sizes[ft_key];

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cms_baseline.increment(data + i);
	}

	//////////////////

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file_error_on_arrival;
	string file_name_error_on_arrival = "test_fg_salsa_count_sketch_all_error_on_arrival_tb_";
	file_name_error_on_arrival.append(to_string(salsa_bits));
	file_name_error_on_arrival.append("_seed_");
	file_name_error_on_arrival.append(to_string(seed));
	file_name_error_on_arrival.append("_downsamplings_per_merge_");
	file_name_error_on_arrival.append(to_string(downsamplings_per_merge));
	file_name_error_on_arrival.append(".txt");
	results_file_error_on_arrival.open(file_name_error_on_arrival, ofstream::out | ofstream::app);
	results_file_error_on_arrival << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;

	///////////////////

	vector<pair<uint64_t, uint64_t>> true_top_k(HH_num);
	partial_sort_copy(true_sizes.begin(),
		true_sizes.end(),
		true_top_k.begin(),
		true_top_k.end(),
		[](pair<uint64_t, uint64_t> const& l,
			pair<uint64_t, uint64_t> const& r)
	{
		return l.second > r.second;
	});

	vector<pair<string, uint64_t>> fg_salsa_hh_top_k = cms_baseline.HH();
	vector<pair<string, uint64_t>> sorted_fg_salsa_hh_top_k(HH_num);
	partial_sort_copy(fg_salsa_hh_top_k.begin(),
		fg_salsa_hh_top_k.end(),
		sorted_fg_salsa_hh_top_k.begin(),
		sorted_fg_salsa_hh_top_k.end(),
		[](pair<string, uint64_t> const& l,
			pair<string, uint64_t> const& r)
	{
		return l.second > r.second;
	});

	int* salsa_fg_hits = new int[HH_num]();

	for (int hhn = 0; hhn < HH_num; ++hhn)
	{
		for (int i = 0; i <= hhn; ++i)
		{
			uint64_t ft_key_salsa_fg = (((uint64_t)ft_to_bobkey_1.run(sorted_fg_salsa_hh_top_k[i].first.c_str(), FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(sorted_fg_salsa_hh_top_k[i].first.c_str(), FT_SIZE);
			salsa_fg_hits[hhn] += (true_sizes[ft_key_salsa_fg] >= true_top_k[hhn].second);
		}
	}

	ofstream results_file_hh;
	string file_name_hh = "test_fg_salsa_count_sketch_all_hh_tb_";
	file_name_hh.append(to_string(salsa_bits));
	file_name_hh.append("_seed_");
	file_name_hh.append(to_string(seed));
	file_name_hh.append("_downsamplings_per_merge_");
	file_name_hh.append(to_string(downsamplings_per_merge));
	file_name_hh.append(".txt");
	results_file_hh.open(file_name_hh, ofstream::out | ofstream::app);
	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\t";
	for (int hhn = 0; hhn < HH_num; ++hhn)
	{
		results_file_hh << salsa_fg_hits[hhn] << "\t";
	}
	results_file_hh << endl;

	delete[] salsa_fg_hits;

	///////////////////

	long double true_l2 = 0;
	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		true_l2 += it->second*it->second;
	}
	true_l2 = sqrtl(true_l2);

	ofstream results_file_l2;
	string file_name_l2 = "test_fg_salsa_count_sketch_all_l2_tb_";
	file_name_l2.append(to_string(salsa_bits));
	file_name_l2.append("_seed_");
	file_name_l2.append(to_string(seed));
	file_name_l2.append("_downsamplings_per_merge_");
	file_name_l2.append(to_string(downsamplings_per_merge));
	file_name_l2.append(".txt");
	results_file_l2.open(file_name_l2, ofstream::out | ofstream::app);
	results_file_l2 << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\t";
	results_file_l2 << "TrueL2\t" << true_l2 << "\tSaslaL2\t" << cms_baseline.l2_estimation() << endl;

}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
