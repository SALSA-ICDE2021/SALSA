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

#include "CMS.hpp"
#include "SalsaCMSBaseline.hpp"
#include "TangoCMSBaseline.hpp"
#include "SalsaCMS.hpp"
#include "CountMinAEE.hpp"
#include "TangoCMS.hpp"

#include "CMSTests.hpp"

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data)
{

	CountMinBaselineSanity cms_baseline;
	cms_baseline.initialize(width, height, seed);

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

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_baseline_cms_error_on_arrival_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_baseline_cms_speed(int N, int width, int height, int seed, const char* data)
{

	CountMinBaseline cms_baseline;
	cms_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_baseline_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

/* Weighted */

void test_weighted_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data, const uint16_t* data_weights)
{

	WeightedCountMinBaseline wcms_baseline;
	wcms_baseline.initialize(width, height, seed);

	unordered_map<uint64_t, uint64_t> true_sizes;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	int64_t weight_index = 0;

	int64_t sum_weights = 0;

	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		long double point_error = (int64_t)wcms_baseline.query(data + i) - (int64_t)true_sizes[ft_key];

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += data_weights[weight_index];
		wcms_baseline.add(data + i, data_weights[weight_index]);

		sum_weights += data_weights[weight_index];

		++weight_index;
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_weighted_baseline_cms_error_on_arrival_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tSum\t" << sum_weights << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_weighted_baseline_cms_speed(int N, int width, int height, int seed, const char* data, const uint16_t* data_weights)
{

	WeightedCountMinBaseline wcms_baseline;
	wcms_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;
	int64_t weight_index = -1;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		wcms_baseline.add(data + i, data_weights[++weight_index]);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_weighted_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	results_file.open("test_weighted_baseline_cms_speed.txt", ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_baseline_cus_error_on_arrival(int N, int width, int height, int seed, const char* data)
{

	ConservativeUpdateBaselineSanity cus_baseline;
	cus_baseline.initialize(width, height, seed);

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

		long double point_error = (int64_t)cus_baseline.query(data + i) - (int64_t)true_sizes[ft_key];

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cus_baseline.increment(data + i);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_baseline_cus_error_on_arrival_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_baseline_cus_speed(int N, int width, int height, int seed, const char* data)
{

	ConservativeUpdateBaseline cus_baseline;
	cus_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cus_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_baseline_cus_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_baseline_cus_speed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_salsa_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data)
{

	SalsaCMSBaseline cms_baseline;
	cms_baseline.initialize(width, height, seed);

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

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_salsa_baseline_cms_error_on_arrival_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_salsa_baseline_cms_speed(int N, int width, int height, int seed, const char* data)
{

	SalsaCMSBaseline cms_baseline;
	cms_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_salsa_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_salsa_baseline_cms_speed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_maximum_salsa_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data)
{

	MaximumSalsaCMSBaselineSanity cms_baseline;
	cms_baseline.initialize(width, height, seed);

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

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_maximum_salsa_baseline_cms_error_on_arrival_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_maximum_salsa_baseline_cms_speed(int N, int width, int height, int seed, const char* data)
{

	MaximumSalsaCMSBaseline cms_baseline;
	cms_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_maximum_salsa_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_maximum_salsa_baseline_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_fg_salsa_baseline_cms_error_on_arrival(int salsa_bits, int N, int width, int height, int seed, const char* data)
{

	FineGrainedSalsaCMSBaselineSanity cms_baseline(salsa_bits, width, height, seed);

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

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string file_name = "test_fg_salsa_baseline_cms_error_on_arrival_";
	file_name.append("salsa_bits_");
	file_name.append(to_string(salsa_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_fg_salsa_baseline_cms_speed(int salsa_bits, int N, int width, int height, int seed, const char* data)
{

	FineGrainedSalsaCMSBaseline cms_baseline(salsa_bits, width, height, seed);


	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_fg_salsa_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string file_name = "test_fg_salsa_baseline_cms_speed_";
	file_name.append("salsa_bits_");
	file_name.append(to_string(salsa_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_salsa_cms_error_on_arrival(int N, int width, int height, int seed, const char* data, int downsamplings_per_merge)
{

	SalsaCMSSanity cms_baseline(width, height, seed, downsamplings_per_merge);

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

		//cout << point_error << endl;

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
	string fn = "test_salsa_cms_error_on_arrival_";
	fn.append(to_string(seed));
	fn.append("_downsamplings_per_merge_");
	fn.append(to_string(downsamplings_per_merge));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_salsa_cms_speed(int N, int width, int height, int seed, const char* data, int downsamplings_per_merge)
{

	SalsaCMSSanity cms_baseline(width, height, seed, downsamplings_per_merge);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_salsa_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	results_file.open("test_salsa_cms_speed.txt", ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_fg_salsa_cms_error_on_arrival(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge)
{

	FineGrainedSalsaCMSSanity cms_baseline(salsa_bits, width, height, seed, downsamplings_per_merge, false);

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

		//cout << point_error << " " << ft_key << " " << (int64_t)cms_baseline.query(data + i) << " " << (int64_t)true_sizes[ft_key] << endl;

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cms_baseline.increment(data + i);
		//cms_baseline.query(data + 0);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file; 
	string file_name = "test_fg_salsa_cms_error_on_arrival_tb_";
	file_name.append(to_string(salsa_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append("_downsamplings_per_merge_");
	file_name.append(to_string(downsamplings_per_merge));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_fg_salsa_split_cms_error_on_arrival(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge)
{

	FineGrainedSalsaCMSSanity cms_baseline(salsa_bits, width, height, seed, downsamplings_per_merge, true);

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

		//cout << point_error << " " << ft_key << " " << (int64_t)cms_baseline.query(data + i) << " " << (int64_t)true_sizes[ft_key] << endl;

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cms_baseline.increment(data + i);
		//cms_baseline.query(data + 0);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string file_name = "test_fg_salsa_split_cms_error_on_arrival_tb_";
	file_name.append(to_string(salsa_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append("_downsamplings_per_merge_");
	file_name.append(to_string(downsamplings_per_merge));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_fg_salsa_cms_speed(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge)
{

	FineGrainedSalsaCMSSanity cms_baseline(salsa_bits, width, height, seed, downsamplings_per_merge, false);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_fg_salsa_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	results_file.open("test_fg_salsa_cms_speed.txt", ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_aee_cms_error_on_arrival(int N, int width, int height, int seed, const char* data)
{

	CountMinCCSanity cms_baseline(width, height, seed);

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

		//cout << point_error << endl;

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
	string fn = "test_aee_cms_error_on_arrival_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_aee_cms_speed(int N, int width, int height, int seed, const char* data)
{

	CountMinCC cms_baseline(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_aee_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_aee_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_aee_max_speed_cms_error_on_arrival(int N, int width, int height, int seed, const char* data, double epsilon_cc, double delta_cc)
{

	CountMinCC_MaxSpeed cms_baseline(width, height, seed, epsilon_cc, delta_cc);

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

		//cout << point_error << endl;

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
	string fn = "test_aee_max_speed_cms_error_on_arrival_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_aee_max_speed_cms_speed(int N, int width, int height, int seed, const char* data, double epsilon_cc, double delta_cc)
{

	CountMinCC_MaxSpeed cms_baseline(width, height, seed, epsilon_cc, delta_cc);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_aee_max_speed_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_aee_max_speed_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_maximum_tango_baseline_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data)
{

	TangoBaselineSanity cms_baseline;
	cms_baseline.initialize(tango_bits, width, height, seed);

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

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string file_name = "test_maximum_tango_baseline_cms_error_on_arrival_";
	file_name.append("tango_bits_");
	file_name.append(to_string(tango_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_maximum_tango_baseline_cms_speed(int tango_bits, int N, int width, int height, int seed, const char* data)
{

	TangoBaseline cms_baseline;
	cms_baseline.initialize(tango_bits, width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_maximum_tango_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string file_name = "test_maximum_tango_baseline_cms_speed_";
	file_name.append("tango_bits_");
	file_name.append(to_string(tango_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_maximum_smartango_baseline_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data)
{

	smarTangoBaselineSanity cms_baseline;
	cms_baseline.initialize(tango_bits, width, height, seed);

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

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string file_name = to_string(tango_bits);
	file_name.append("_");
	file_name.append(to_string(seed));
	file_name.append("_");
	file_name.append("test_maximum_smartango_baseline_cms_error_on_arrival.txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_smartango_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge)
{

	smarTango cms_baseline(tango_bits, width, height, seed, downsamplings_per_merge);

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

		//cout << point_error << " " << ft_key << " " << (int64_t)cms_baseline.query(data + i) << " " << (int64_t)true_sizes[ft_key] << endl;

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cms_baseline.increment(data + i);
		//cms_baseline.query(data + 0);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string file_name = "test_smartango_cms_error_on_arrival_tb_";
	file_name.append(to_string(tango_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append("_downsamplings_per_merge_");
	file_name.append(to_string(downsamplings_per_merge));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_tango_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge)
{

	Tango cms_baseline(tango_bits, width, height, seed, downsamplings_per_merge);

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

		//cout << point_error << " " << ft_key << " " << (int64_t)cms_baseline.query(data + i) << " " << (int64_t)true_sizes[ft_key] << endl;

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cms_baseline.increment(data + i);
		//cms_baseline.query(data + 0);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string file_name = "test_tango_cms_error_on_arrival_tb_";
	file_name.append(to_string(tango_bits));
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

void test_maximum_salsa_baseline_cus_error_on_arrival(int N, int width, int height, int seed, const char* data)
{

	MaximumSalsaCUSBaseline cus_baseline;
	cus_baseline.initialize(width, height, seed);

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

		long double point_error = (int64_t)cus_baseline.query(data + i) - (int64_t)true_sizes[ft_key];

		L1e += abs(point_error);
		L2e += (point_error * point_error);
		L_max = (L_max < abs(point_error)) ? abs(point_error) : L_max;

		true_sizes[ft_key] += 1;
		cus_baseline.increment(data + i);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);

	ofstream results_file;
	string fn = "test_maximum_salsa_baseline_cus_error_on_arrival_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_maximum_salsa_baseline_cus_speed(int N, int width, int height, int seed, const char* data)
{

	MaximumSalsaCUSBaseline cus_baseline;
	cus_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cus_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_maximum_salsa_baseline_cus_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_maximum_salsa_baseline_cus_speed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_sanity_cms(int N, int width, int height, int seed, const char* data)
{

	CountMinBaselineSanity cms_baseline_sanity;
	cms_baseline_sanity.initialize(width, height, seed);

	MaximumSalsaCMSBaselineSanity maximum_salsa_cms_baseline_sanity;
	maximum_salsa_cms_baseline_sanity.initialize(width << 2, height, seed);

	TangoBaselineSanity tango_1_bit;
	tango_1_bit.initialize(1, width << 5, height, seed);

	TangoBaselineSanity tango_2_bit;
	tango_2_bit.initialize(2, width << 4, height, seed);

	TangoBaselineSanity tango_4_bit;
	tango_4_bit.initialize(4, width << 3, height, seed);

	TangoBaselineSanity tango_8_bit;
	tango_8_bit.initialize(8, width << 2, height, seed);

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

		long double point_error_1 = (int64_t)cms_baseline_sanity.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_2 = (int64_t)maximum_salsa_cms_baseline_sanity.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_3 = (int64_t)tango_8_bit.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_4 = (int64_t)tango_4_bit.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_5 = (int64_t)tango_2_bit.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_6 = (int64_t)tango_1_bit.query(data + i) - (int64_t)true_sizes[ft_key];

		if ((abs(point_error_2) > abs(point_error_1)) || 
			(abs(point_error_3) > abs(point_error_2)) || 
			(abs(point_error_4) > abs(point_error_2)) ||
			(abs(point_error_5) > abs(point_error_2)) ||
			(abs(point_error_6) > abs(point_error_2)))

		{
			cout << "we are not BOB! BOB is guilty! " << i << endl;
			cout << "we are not BOB! BOB is guilty! " << ft_key << endl;
			(int64_t)tango_4_bit.query(data + i);
			(int64_t)tango_2_bit.query(data + i);
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
		cms_baseline_sanity.increment(data + i);
		maximum_salsa_cms_baseline_sanity.increment(data + i);
		tango_8_bit.increment(data + i);
		tango_4_bit.increment(data + i);
		tango_2_bit.increment(data + i);
		tango_1_bit.increment(data + i);

	}
}

void test_sanity_salsa_cms(int N, int width, int height, int seed, const char* data)
{

	CountMinBaselineSanity cms_baseline_sanity;
	cms_baseline_sanity.initialize(width, height, seed);

	MaximumSalsaCMSBaselineSanity maximum_salsa_cms_baseline_sanity;
	maximum_salsa_cms_baseline_sanity.initialize(width << 2, height, seed);

	FineGrainedSalsaCMSBaselineSanity salsa_1_bit(1, width << 5, height, seed);
	FineGrainedSalsaCMSBaselineSanity salsa_2_bit(2, width << 4, height, seed);
	FineGrainedSalsaCMSBaselineSanity salsa_4_bit(4, width << 3, height, seed);
	FineGrainedSalsaCMSBaselineSanity salsa_8_bit(8, width << 2, height, seed);

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

		long double point_error_1 = (int64_t)cms_baseline_sanity.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_2 = (int64_t)maximum_salsa_cms_baseline_sanity.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_3 = (int64_t)salsa_8_bit.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_4 = (int64_t)salsa_4_bit.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_5 = (int64_t)salsa_2_bit.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_6 = (int64_t)salsa_1_bit.query(data + i) - (int64_t)true_sizes[ft_key];

		if ((abs(point_error_2) > abs(point_error_1))  ||
			(abs(point_error_3) != abs(point_error_2)) ||
			(abs(point_error_4) > abs(point_error_3))  ||
			(abs(point_error_5) > abs(point_error_4))  ||
			(abs(point_error_6) > abs(point_error_5)))

		{
			cout << "we are not BOB! BOB is guilty! " << i << endl;
			//cout << "we are not BOB! BOB is guilty! " << ft_key << endl;
			//(int64_t)salsa_4_bit.query(data + i);
			//(int64_t)salsa_2_bit.query(data + i);
			system("pause");
			//exit(1);
		}
		//int query_i = 126321;
		//cout << i << "\t" << salsa_1_bit.read_row_value(1800, 2) << " " << salsa_1_bit.read_row_value(1804, 2) << endl;
		//cout << i << "\t" << salsa_1_bit.read_row_value_ft(data + query_i, 0) << " "
		//	<< salsa_1_bit.read_row_value_ft(data + query_i, 1) << " "
		//	<< salsa_1_bit.read_row_value_ft(data + query_i, 2) << " "
		//	<< salsa_1_bit.read_row_value_ft(data + query_i, 3) << " "
		//	<< endl;

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
		cms_baseline_sanity.increment(data + i);
		maximum_salsa_cms_baseline_sanity.increment(data + i);
		salsa_8_bit.increment(data + i);
		salsa_4_bit.increment(data + i);
		salsa_2_bit.increment(data + i);
		salsa_1_bit.increment(data + i);

	}
}

void test_sanity_cus(int N, int width, int height, int seed, const char* data)
{

	ConservativeUpdateBaselineSanity cus_baseline_sanity;
	cus_baseline_sanity.initialize(width, height, seed);

	MaximumSalsaCMSBaseline maximum_salsa_cms_baseline;
	maximum_salsa_cms_baseline.initialize(width << 2, height, seed);

	MaximumSalsaCUSBaseline maximum_salsa_cus_baseline;
	maximum_salsa_cus_baseline.initialize(width << 2, height, seed);

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

		long double point_error_1 = (int64_t)cus_baseline_sanity.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_2 = (int64_t)maximum_salsa_cms_baseline.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_3 = (int64_t)maximum_salsa_cus_baseline.query(data + i) - (int64_t)true_sizes[ft_key];

		if ((abs(point_error_3) > abs(point_error_1)) || (abs(point_error_3) > abs(point_error_2)))
		{
			cout << "we are not BOB! BOB is guilty! " << i << endl;
			system("pause");
			exit(1);
		}

		true_sizes[ft_key] += 1;
		cus_baseline_sanity.increment(data + i);
		maximum_salsa_cms_baseline.increment(data + i);
		maximum_salsa_cus_baseline.increment(data + i);
	}
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_fg_salsa_split_counters_sanity_cms(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge)
{

	FineGrainedSalsaCMSSanity w_split(salsa_bits, width, height, seed, downsamplings_per_merge, true);
	FineGrainedSalsaCMSSanity n_split(salsa_bits, width, height, seed, downsamplings_per_merge, false);

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

		long double point_error_1 = (int64_t)n_split.query(data + i) - (int64_t)true_sizes[ft_key];
		long double point_error_2 = (int64_t)w_split.query(data + i) - (int64_t)true_sizes[ft_key];

		if (abs(point_error_2) > abs(point_error_1))

		{
			cout << "we are not BOB! BOB is guilty! " << i << endl;
			cout << "we are not BOB! BOB is guilty! " << ft_key << endl;
			cout << (int64_t)n_split.query(data + i) << endl;
			cout << (int64_t)w_split.query(data + i) << endl;
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

		int p = 21888321;
		
		if (i >= 21880321)
		{
			cout << "no split: " << n_split.read_row_value_ft(data + p, 0) << " ";
			cout << n_split.read_row_value_ft(data + p, 1) << " ";
			cout << n_split.read_row_value_ft(data + p, 2) << " ";
			cout << n_split.read_row_value_ft(data + p, 3) << endl;

			cout << "wi split: " << w_split.read_row_value_ft(data + p, 0) << " ";
			cout << w_split.read_row_value_ft(data + p, 1) << " ";
			cout << w_split.read_row_value_ft(data + p, 2) << " ";
			cout << w_split.read_row_value_ft(data + p, 3) << "\t\t" << i << endl;
		}


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
		n_split.increment(data + i);
		w_split.increment(data + i);

	}
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_fg_salsa_baseline_cms_count_distinct(int salsa_bits, int N, int width, int height, int seed, const char* data)
{
	FineGrainedSalsaCMSBaselineSanity salsa(salsa_bits, width, height, seed);

	unordered_map<uint64_t, uint64_t> true_sizes;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	int mycounter = 0;
	for (int64_t i = 0; i < stop_loop; i += 13)
	{

		uint64_t ft_key = (((uint64_t)ft_to_bobkey_1.run(data + i, FT_SIZE)) << 32) + (uint64_t)ft_to_bobkey_2.run(data + i, FT_SIZE);

		if (true_sizes.find(ft_key) == true_sizes.end())
		{
			true_sizes[ft_key] = 0;
		}

		true_sizes[ft_key] += 1;
		salsa.increment(data + i);
	}

	ofstream results_file;
	string file_name = "test_fg_salsa_baseline_cms_count_distinct_salsa_bits_";
	file_name.append(to_string(salsa_bits));
	file_name.append("_seed_");
	file_name.append(to_string(seed));
	file_name.append(".txt");
	results_file.open(file_name, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height;

	results_file << "\tTrueCountDistinct\t" << true_sizes.size();
	results_file << "\tRowMedian\t" << salsa.provable_query_count_distinct();
	results_file << "\tGlobal\t" << salsa.global_query_count_distinct();
	results_file << "\tBestGuess\t" << salsa.best_guess_query_count_distinct();
	results_file << "\tBoundAverage\t" << salsa.bound_average_query_count_distinct();
	results_file << "\tBest_Guess\t" << salsa.global_best_guess_query_count_distinct();
	results_file << endl;
	
	//cout << "RS: Salsa " << salsa_bits << " Bit Count Distinct:" << salsa.provable_query_count_distinct() << endl;
	//cout << "MM: Salsa " << salsa_bits << " Bit Count Distinct:" << salsa.global_query_count_distinct() << endl;
	//cout << "BG: Salsa " << salsa_bits << " Bit Count Distinct:" << salsa.best_guess_query_count_distinct() << endl;
	//cout << "AV: Salsa " << salsa_bits << " Bit Count Distinct:" << salsa.bound_average_query_count_distinct() << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_analytics_salsa_cms_error_on_arrival(int N, int width, int height, int seed, const char* data, double delta, int downsamplings_before_merge)
{

	SalsaAnalyticalErrorCMS cms_baseline(width, height, seed, delta, downsamplings_before_merge);

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

		//cout << point_error << endl;

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
	string fn = "test_analytics_salsa_cms_error_on_arrival_";
	fn.append(to_string(seed));
	fn.append("_delta_");
	fn.append(to_string(delta));
	fn.append("_dbm_");
	fn.append(to_string(downsamplings_before_merge));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_analytics_salsa_cms_speed(int N, int width, int height, int seed, const char* data, double delta, int downsamplings_before_merge)
{
	SalsaAnalyticalErrorCMS cms_baseline(width, height, seed, delta, downsamplings_before_merge);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_analytics_salsa_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_analytics_salsa_cms_speed_";
	fn.append(to_string(seed));
	fn.append("_delta_");
	fn.append(to_string(delta));
	fn.append("_dbm_");
	fn.append(to_string(downsamplings_before_merge));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_analytics_fg_salsa_cms_error_on_arrival(int N, int counterSize, int width, int height, int seed, const char* data, int downsamplings_before_merge, bool split_counters_flag, double delta)
{

	FineGrainedSalsaAnalyticalErrorCMSSanity cms_baseline(counterSize, width, height, seed, downsamplings_before_merge, split_counters_flag, delta);

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

		//cout << point_error << endl;

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
	string fn = "test_analytics_fg_salsa_cms_error_on_arrival_";
	fn.append(to_string(seed));
	fn.append("_delta_");
	fn.append(to_string(delta));
	fn.append("_dbm_");
	fn.append(to_string(downsamplings_before_merge));
	fn.append("_salsa_bits_");
	fn.append(to_string(counterSize));
	fn.append("_splitting_");
	fn.append(to_string(split_counters_flag));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////

void test_analytics_fg_salsa_cms_final_error(int N, int counterSize, int width, int height, int seed, const char* data, int downsamplings_before_merge, bool split_counters_flag, double delta)
{

	FineGrainedSalsaAnalyticalErrorCMSSanity cms_baseline(counterSize, width, height, seed, downsamplings_before_merge, split_counters_flag, delta);

	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	unordered_map<double, vector<uint64_t>> threshold_to_hh_ft_keys;
	vector<double> thresholds = {0.0001, 0.000178, 0.000316, 0.000562, 0.001, 0.00178, 0.00316, 0.00562, 0.01};
	for (int i = 0; i < thresholds.size(); ++i)
	{
		threshold_to_hh_ft_keys[thresholds[i]] = vector<uint64_t>();
	}

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
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

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = cms_baseline.query(data + ft_key_i_values[it->first]) - it->second;
		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;
		L_max = L_max > error ? L_max : error;

		for (int i = 0; i < thresholds.size(); ++i)
		{
			if (it->second >= thresholds[i] * N)
			{
				threshold_to_hh_ft_keys[thresholds[i]].push_back(it->first);
			}
			else 
			{
				break;
			}
		}
	}

	vector<long double> HHre;
	for (int i = 0; i < thresholds.size(); ++i)
	{
		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			long double current_query = cms_baseline.query(data + ft_key_i_values[current_ftkey]);
			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];
		}
		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);
	RL1e /= N;

	ofstream results_file;
	string fn = "test_analytics_fg_salsa_cms_final_relative_error_";
	fn.append(to_string(seed));
	fn.append("_delta_");
	fn.append(to_string(delta));
	fn.append("_dbm_");
	fn.append(to_string(downsamplings_before_merge));
	fn.append("_salsa_bits_");
	fn.append(to_string(counterSize));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;

	ofstream results_file_hh;
	string fn_hh = "test_analytics_fg_salsa_cms_final_relative_error_hh_";
	fn_hh.append(to_string(seed));
	fn_hh.append("_delta_");
	fn_hh.append(to_string(delta));
	fn_hh.append("_dbm_");
	fn_hh.append(to_string(downsamplings_before_merge));
	fn_hh.append("_salsa_bits_");
	fn_hh.append(to_string(counterSize));
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i];
	}
	results_file_hh << endl;

}

void test_baseline_fg_salsa_cms_final_error(int N, int counterSize, int width, int height, int seed, const char* data)
{

	FineGrainedSalsaCMSBaselineSanity cms_baseline(counterSize, width, height, seed);

	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	unordered_map<double, vector<uint64_t>> threshold_to_hh_ft_keys;
	vector<double> thresholds = { 0.0001, 0.000178, 0.000316, 0.000562, 0.001, 0.00178, 0.00316, 0.00562, 0.01 };
	for (int i = 0; i < thresholds.size(); ++i)
	{
		threshold_to_hh_ft_keys[thresholds[i]] = vector<uint64_t>();
	}

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
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

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = cms_baseline.query(data + ft_key_i_values[it->first]) - it->second;
		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;
		L_max = L_max > error ? L_max : error;

		for (int i = 0; i < thresholds.size(); ++i)
		{
			if (it->second >= thresholds[i] * N)
			{
				threshold_to_hh_ft_keys[thresholds[i]].push_back(it->first);
			}
			else
			{
				break;
			}
		}
	}

	vector<long double> HHre;
	for (int i = 0; i < thresholds.size(); ++i)
	{
		long double current_hh_rl = 0;
		for (int j = 0; j < threshold_to_hh_ft_keys[thresholds[i]].size(); ++j)
		{
			uint64_t current_ftkey = threshold_to_hh_ft_keys[thresholds[i]][j];
			long double current_query = cms_baseline.query(data + ft_key_i_values[current_ftkey]);
			current_hh_rl += (long double)abs(current_query - true_sizes[current_ftkey]) / (long double)true_sizes[current_ftkey];
		}
		current_hh_rl /= threshold_to_hh_ft_keys[thresholds[i]].size();
		HHre.push_back(current_hh_rl);
	}

	L1e /= N;
	L2e /= N;
	L2e = sqrt(L2e);
	RL1e /= N;

	ofstream results_file;
	string fn = "test_baseline_fg_salsa_cms_final_relative_error_";
	fn.append(to_string(seed));
	fn.append("_salsa_bits_");
	fn.append(to_string(counterSize));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;

	ofstream results_file_hh;
	string fn_hh = "test_baseline_fg_salsa_cms_final_relative_error_hh_";
	fn_hh.append(to_string(seed));
	fn_hh.append("_salsa_bits_");
	fn_hh.append(to_string(counterSize));
	fn_hh.append(".txt");
	results_file_hh.open(fn_hh, ofstream::out | ofstream::app);

	results_file_hh << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height;
	for (int i = 0; i < HHre.size(); ++i)
	{
		results_file_hh << "\tThreshold\t" << thresholds[i] << "\tRelError\t" << HHre[i];
	}
	results_file_hh << endl;

}

void test_salsa_cms_final_error_histogram(int N, int width, int height, int seed, const char* data)
{

	MaximumSalsaCMSBaseline cms_baseline;
	cms_baseline.initialize(width, height, seed);

	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
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

	ofstream results_histogram_file;
	string fn_histogram = "test_salsa_cms_final_error_histogram_";
	fn_histogram.append(to_string(seed));
	fn_histogram.append(".txt");
	results_histogram_file.open(fn_histogram, ofstream::out | ofstream::app);

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = (long long)cms_baseline.query(data + ft_key_i_values[it->first]) - (long long)it->second;
		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;
		L_max = L_max > error ? L_max : error;
		results_histogram_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tsize\t" << (long long)it->second << "\tError\t" << error << endl;
	}

	L1e /= true_sizes.size();
	L2e /= true_sizes.size();
	L2e = sqrt(L2e);
	RL1e /= true_sizes.size();

	ofstream results_file;
	string fn = "test_salsa_cms_final_error_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;

}

void test_baseline_cms_final_error_histogram(int N, int width, int height, int seed, const char* data)
{

	CountMinBaseline cms_baseline;
	cms_baseline.initialize(width, height, seed);

	unordered_map<uint64_t, uint64_t> true_sizes;
	unordered_map<uint64_t, uint64_t> ft_key_i_values;

	int64_t stop_loop = N * FT_SIZE;

	long double L1e = 0, L2e = 0, L_max = 0, RL1e = 0;

	BOBHash ft_to_bobkey_1(seed);
	BOBHash ft_to_bobkey_2(seed + 17);

	for (int64_t i = 0; i < stop_loop; i += 13)
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

	ofstream results_histogram_file;
	string fn_histogram = "test_baseline_cms_final_error_histogram_";
	fn_histogram.append(to_string(seed));
	fn_histogram.append(".txt");
	results_histogram_file.open(fn_histogram, ofstream::out | ofstream::app);

	for (auto it = true_sizes.begin(); it != true_sizes.end(); it++)
	{
		long double error = (long long)cms_baseline.query(data + ft_key_i_values[it->first]) - (long long)it->second;
		L1e += abs(error);
		L2e += error * error;
		RL1e += abs(error) / it->second;
		L_max = L_max > error ? L_max : error;
		results_histogram_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tsize\t" << (long long)it->second << "\tError\t" << error << endl;
	}

	L1e /= true_sizes.size();
	L2e /= true_sizes.size();
	L2e = sqrt(L2e);
	RL1e /= true_sizes.size();

	ofstream results_file;
	string fn = "test_baseline_cms_final_error_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << "\tLR1e Error\t" << RL1e << endl;

}




