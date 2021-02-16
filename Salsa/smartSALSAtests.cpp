#include "smartSALSAtests.h"



void test_SmartEncodingSALSA_baseline_cms_sanity(int N, int width, int height, int seed, const char* data)
{
	//SmartEncodingSALSASanity cms_baseline;
	SmartEncodingSALSA smartEncodingSALSA;
	FastEncodingSALSA fastEncodingSALSA;
	MaximumSalsaCMSBaseline cms_baseline;
	StupidEncodingSALSA stupidEncodingSALSA;
	StupiderEncodingSALSA stupiderEncodingSALSA;


	cms_baseline.initialize(width, height, seed);
	smartEncodingSALSA.initialize(width, height, seed);
	fastEncodingSALSA.initialize(width, height, seed);
	stupidEncodingSALSA.initialize(width, height, seed);
	stupiderEncodingSALSA.initialize(width, height, seed);

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
		//if (i == 52845689) {
		//if (fastEncodingSALSA.query(data + 53460238) == 96){
		if (i == 832832) {
			cout << "break here" << endl;
		}

		cms_baseline.increment(data + i);
		smartEncodingSALSA.increment(data + i);
		fastEncodingSALSA.increment(data + i);
		stupidEncodingSALSA.increment(data + i);
		stupiderEncodingSALSA.increment(data + i);
		
		//if (fastEncodingSALSA.indexof(data + i, 0) == fastEncodingSALSA.indexof(data + 51974, 0)) {
		//	cout << "11485, " << i << " " << cms_baseline.query(data + 53460238) << " " << smartEncodingSALSA.query(data + 53460238) << " " << fastEncodingSALSA.query(data + 53460238) << endl;
		//}
		/*
		if (fastEncodingSALSA.indexof(data + i, 0) == fastEncodingSALSA.indexof(data + 51974, 0)) {
			cout << "11486, " << i << " " << cms_baseline.query(data + 53460238) << " " << smartEncodingSALSA.query(data + 53460238) << " " << fastEncodingSALSA.query(data + 53460238) << endl;
		}
		if (fastEncodingSALSA.indexof(data + i, 0) == fastEncodingSALSA.indexof(data + 51974, 0)) {
			cout << "11487, " << i << " " << cms_baseline.query(data + 53460238) << " " << smartEncodingSALSA.query(data + 53460238) << " " << fastEncodingSALSA.query(data + 53460238) << endl;
		}*/

		int cms_baseline_query = cms_baseline.query(data + i);
		int smartEncodingSALSA_query = smartEncodingSALSA.query(data + i);
		int fastEncodingSALSA_query = fastEncodingSALSA.query(data + i);
		int stupidEncodingSALSA_query = stupidEncodingSALSA.query(data + i);
		int stupiderEncodingSALSA_query = stupiderEncodingSALSA.query(data + i);

		if ((cms_baseline_query != smartEncodingSALSA_query) ||
			(cms_baseline_query != fastEncodingSALSA_query)  ||
			(cms_baseline_query != stupiderEncodingSALSA_query) ||
			(cms_baseline_query != stupidEncodingSALSA_query)) {
			cout << i << " " << cms_baseline.query(data + i) << " " << smartEncodingSALSA.query(data + i) 
				<< " " << fastEncodingSALSA.query(data + i) << " " << stupidEncodingSALSA.query(data + i)
				<< " " << stupiderEncodingSALSA.query(data + i) << endl;
			system("pause");
		}
		if ((i % 13000000) == 0) {
			cout << i / 13 << endl;
		}
	}
}

void test_SmartEncodingSALSA_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data)
{

	//SmartEncodingSALSASanity cms_baseline;
	SmartEncodingSALSA cms_baseline;

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
	string fn = "test_SmartEncodingSALSA_baseline_cms_error_on_arrival_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tL1 Error\t" << L1e << "\tL2 Error\t" << L2e << "\tL(inf) Error\t" << L_max << endl;
}

void test_SmartEncodingSALSA_baseline_cms_speed(int N, int width, int height, int seed, const char* data)
{

	SmartEncodingSALSA cms_baseline;
	cms_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_SmartEncodingSALSA_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_SmartEncodingSALSA_baseline_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_FastEncodingSALSA_baseline_cms_speed(int N, int width, int height, int seed, const char* data)
{

	FastEncodingSALSA cms_baseline;
	cms_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_FastEncodingSALSA_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_FastEncodingSALSA_baseline_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_StupidEncodingSALSA_baseline_cms_speed(int N, int width, int height, int seed, const char* data)
{

	StupidEncodingSALSA cms_baseline;
	cms_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_StupidEncodingSALSA_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_StupidEncodingSALSA_baseline_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

void test_StupiderEncodingSALSA_baseline_cms_speed(int N, int width, int height, int seed, const char* data)
{

	StupiderEncodingSALSA cms_baseline;
	cms_baseline.initialize(width, height, seed);

	int64_t stop_loop = N * FT_SIZE;

	auto start = chrono::steady_clock::now();
	for (int64_t i = 0; i < stop_loop; i += FT_SIZE)
	{
		cms_baseline.increment(data + i);
	}
	auto end = chrono::steady_clock::now();

	auto time = chrono::duration_cast<chrono::microseconds>(end - start).count();
	cout << "test_StupiderEncodingSALSA_baseline_cms_speed: Elapsed time in milliseconds : "
		<< time / 1000
		<< " ms" << endl;

	ofstream results_file;
	string fn = "test_StupiderEncodingSALSA_baseline_cms_speed_seed_";
	fn.append(to_string(seed));
	fn.append(".txt");
	results_file.open(fn, ofstream::out | ofstream::app);
	results_file << "N\t" << N << "\tWidth\t" << width << "\tHeight\t" << height << "\tTime\t" << time << endl;
}

//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////