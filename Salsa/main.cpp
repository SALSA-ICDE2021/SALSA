#pragma warning(disable:4996)

#include <iostream>
#include <stdlib.h>
#include <random>
#include <chrono>
#include <math.h>
#include <fstream>
#include <assert.h>
#include <time.h>
#include <sstream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <windows.h>

#include "Defs.hpp"

#include "CMSTests.hpp"
#include "CountSketchTests.hpp"
#include "UnivMonTests.hpp"
#include "smartSALSAtests.h"

#include "genzipf.h"

using namespace std;

int main(int argc, char* argv[])
{

	if (argc < 8) {
		cout << "Usage Error:\n";
		cout << "argv[1]: int N\n";
		cout << "argv[2]: int Seed\n";
		cout << "argv[3]: int Alg\n";
		cout << "argv[4]: bool weighted_experiment\n";
		cout << "argv[5]: String Trace\n";
		cout << "argv[6]: int RowCounterNum\n";
		cout << "argv[7]: int RowsNum\n";
		system("pause");
		return 0;
	}

	// Read arguments
	int N = atoi(argv[1]);
	int Seed = atoi(argv[2]);
	int Alg = atoi(argv[3]);
	bool weighted_experiment = (atoi(argv[4]) == 0) ? false : true;
	string Trace = string(argv[5]);
	int RowCounterNum = atoi(argv[6]);
	int RowsNum = atoi(argv[7]);

	// Input sanity checks
	assert((weighted_experiment == false) && "no weighted simulations");

	// Print input arguments
	cout << "Received input arguments:" << endl;
	cout << "argv[1]: int N = " << N << endl;
	cout << "argv[2]: int Seed = " << Seed << endl;
	cout << "argv[3]: int Alg = " << Alg << endl;
	cout << "argv[4]: bool weighted_experiment = " << weighted_experiment << endl;
	cout << "argv[5]: String Trace = \"" << Trace << "\"" << endl;
	cout << "argv[6]: Counters in row = \"" << RowCounterNum << "\"" << endl;
	cout << "argv[7]: Number of rows = \"" << RowsNum << "\"" << endl;
	cout << endl;
	
	string path = "F:/Traces/";

	if (Trace.rfind("youtube", 0) == 0)
	{
		fstream fileStream;
		path.append("youtube/");
		path.append(Trace);
	}
	else if (Trace.rfind("zipf", 0) == 0)
	{
		fstream fileStream;
		path.append("zipf/");
		path.append(Trace);
		fileStream.open(path);
		if (fileStream.fail()) {
			// file could not be opened, lets generate it
			float alpha = atof(Trace.c_str() + 4);
			genzipf(path.c_str(), 42, alpha, 1 << 24, 99000000);
			return 0;
		}
		else
		{
			fileStream.close();
		}		
	}
	else
	{
		path.append(Trace);
	}

	if (weighted_experiment)
	{
		path.append("_weighted");
	}

	// Asserts
	assert((Seed < 1229) && "Seed too large for BobHash");

	// Seed 
	srand(Seed);

	// Read trace?
	ifstream f(path, ios::binary);

	char* data = NULL;
	uint16_t* data_weights = NULL;

	int64_t read_so_far = 0;
	int64_t remaining = N;

	data = new char[FT_SIZE * N]();

	if (weighted_experiment)
	{
		assert((N <= 98000000) && "N too large for used traces");

		data_weights = new uint16_t[N]();
		char* temp_data_buffer = new char[100000 * FT_SIZE_WEIGHTS]();

		while (remaining > 0)
		{
			int64_t to_read = remaining > 100000 ? 100000 : remaining;
			f.read(temp_data_buffer, to_read * FT_SIZE_WEIGHTS);

			for (int i = 0; i < to_read; ++i)
			{
				memcpy(data + (i + read_so_far)* FT_SIZE, temp_data_buffer + i * FT_SIZE_WEIGHTS, FT_SIZE);

				uint8_t* p1 = (uint8_t*)(temp_data_buffer + FT_SIZE_WEIGHTS * i + FT_SIZE);
				uint8_t* p2 = (uint8_t*)(temp_data_buffer + FT_SIZE_WEIGHTS * i + FT_SIZE + 1);
				data_weights[read_so_far + i] = (*p1) * 256 + *p2;
			}

			remaining -= to_read;
			read_so_far += to_read;
		}

		delete[] temp_data_buffer;
	}

	else
	{
		while (remaining > 0)
		{
			int64_t to_read = remaining > 100000 ? 100000 : remaining;
			f.read(data + read_so_far * FT_SIZE, to_read * FT_SIZE);
			remaining -= to_read;
			read_so_far += to_read;
		}
	}

	assert((read_so_far == N) && "Read trace file error");
	assert((remaining == 0) && "Read trace file error");

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	if (CreateDirectory(Trace.c_str(), NULL) ||
		ERROR_ALREADY_EXISTS == GetLastError())
	{
		// Enter
		if (SetCurrentDirectoryA(Trace.c_str()))
		{
			cout << "Directory change ok" << endl;
		}	
		else
		{
			cout << "Directory change failed" << endl;
			return 0;
		}			
	}
	else
	{
		// Failed to create directory
		cout << "Failed to create trace directory" << endl;
		return 0;
	}

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	if (Alg == 0)
	{
		/* error */

		// CMS baseline (32-bit counters)
		test_baseline_cms_error_on_arrival(N, RowCounterNum, RowsNum, Seed, data);

		// Salsa sum (start with 8-bit counters)
		test_salsa_baseline_cms_error_on_arrival(N, RowCounterNum << 2, RowsNum, Seed, data);

		// Salsa max (start with 8-bit counters)
		test_maximum_salsa_baseline_cms_error_on_arrival(N, RowCounterNum << 2, RowsNum, Seed, data);

		// Salsa Max
		test_fg_salsa_baseline_cms_error_on_arrival(1, N, RowCounterNum << 5, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_error_on_arrival(2, N, RowCounterNum << 4, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_error_on_arrival(4, N, RowCounterNum << 3, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_error_on_arrival(8, N, RowCounterNum << 2, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_error_on_arrival(16, N, RowCounterNum << 1, RowsNum, Seed, data);

		// Tango Max
		test_maximum_tango_baseline_cms_error_on_arrival(1, N, RowCounterNum << 5, RowsNum, Seed, data);
		test_maximum_tango_baseline_cms_error_on_arrival(2, N, RowCounterNum << 4, RowsNum, Seed, data);
		test_maximum_tango_baseline_cms_error_on_arrival(4, N, RowCounterNum << 3, RowsNum, Seed, data);
		test_maximum_tango_baseline_cms_error_on_arrival(8, N, RowCounterNum << 2, RowsNum, Seed, data);
		test_maximum_tango_baseline_cms_error_on_arrival(16, N, RowCounterNum << 1, RowsNum, Seed, data);
	}
	else if (Alg == 1)
	{
		/* speed */

		// CMS baseline (32-bit counters)
		test_baseline_cms_speed(N, RowCounterNum, RowsNum, Seed, data);

		// Salsa sum (start with 8-bit counters)
		test_salsa_baseline_cms_speed(N, RowCounterNum << 2, RowsNum, Seed, data);

		// Salsa max (start with 8-bit counters)
		test_maximum_salsa_baseline_cms_speed(N, RowCounterNum << 2, RowsNum, Seed, data);

		// Salsa Max
		test_fg_salsa_baseline_cms_speed(1, N, RowCounterNum << 5, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_speed(2, N, RowCounterNum << 4, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_speed(4, N, RowCounterNum << 3, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_speed(8, N, RowCounterNum << 2, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_speed(16, N, RowCounterNum << 1, RowsNum, Seed, data);

		// Tango Max
		test_maximum_tango_baseline_cms_speed(1, N, RowCounterNum << 5, RowsNum, Seed, data);
		test_maximum_tango_baseline_cms_speed(2, N, RowCounterNum << 4, RowsNum, Seed, data);
		test_maximum_tango_baseline_cms_speed(4, N, RowCounterNum << 3, RowsNum, Seed, data);
		test_maximum_tango_baseline_cms_speed(8, N, RowCounterNum << 2, RowsNum, Seed, data);
		test_maximum_tango_baseline_cms_speed(16, N, RowCounterNum << 1, RowsNum, Seed, data);
	}
	else if (Alg == 2)
	{
		// estimators

		// our analytical solution (merge-downsample)
		test_analytics_fg_salsa_cms_error_on_arrival(N, 1, RowCounterNum << 5, RowsNum, Seed, data, 0, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 2, RowCounterNum << 4, RowsNum, Seed, data, 0, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 4, RowCounterNum << 3, RowsNum, Seed, data, 0, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 8, RowCounterNum << 2, RowsNum, Seed, data, 0, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 16, RowCounterNum << 1, RowsNum, Seed, data, 0, 0, 0.001); // for sanity

		test_analytics_fg_salsa_cms_error_on_arrival(N, 1, RowCounterNum << 5, RowsNum, Seed, data, 1, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 2, RowCounterNum << 4, RowsNum, Seed, data, 1, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 4, RowCounterNum << 3, RowsNum, Seed, data, 1, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 8, RowCounterNum << 2, RowsNum, Seed, data, 1, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 16, RowCounterNum << 1, RowsNum, Seed, data, 1, 0, 0.001); // for sanity

		test_analytics_fg_salsa_cms_error_on_arrival(N, 1, RowCounterNum << 5, RowsNum, Seed, data, 2, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 2, RowCounterNum << 4, RowsNum, Seed, data, 2, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 4, RowCounterNum << 3, RowsNum, Seed, data, 2, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 8, RowCounterNum << 2, RowsNum, Seed, data, 2, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 16, RowCounterNum << 1, RowsNum, Seed, data, 2, 0, 0.001); // for sanity

		test_analytics_fg_salsa_cms_error_on_arrival(N, 1, RowCounterNum << 5, RowsNum, Seed, data, 4, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 2, RowCounterNum << 4, RowsNum, Seed, data, 4, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 4, RowCounterNum << 3, RowsNum, Seed, data, 4, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 8, RowCounterNum << 2, RowsNum, Seed, data, 4, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 16, RowCounterNum << 1, RowsNum, Seed, data, 4, 0, 0.001); // for sanity

		test_analytics_fg_salsa_cms_error_on_arrival(N, 1, RowCounterNum << 5, RowsNum, Seed, data, 8, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 2, RowCounterNum << 4, RowsNum, Seed, data, 8, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 4, RowCounterNum << 3, RowsNum, Seed, data, 8, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 8, RowCounterNum << 2, RowsNum, Seed, data, 8, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 16, RowCounterNum << 1, RowsNum, Seed, data, 8, 0, 0.001); // for sanity

		test_analytics_fg_salsa_cms_error_on_arrival(N, 1, RowCounterNum << 5, RowsNum, Seed, data, 10, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 2, RowCounterNum << 4, RowsNum, Seed, data, 10, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 4, RowCounterNum << 3, RowsNum, Seed, data, 10, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 8, RowCounterNum << 2, RowsNum, Seed, data, 10, 0, 0.001); // for sanity
		test_analytics_fg_salsa_cms_error_on_arrival(N, 16, RowCounterNum << 1, RowsNum, Seed, data, 10, 0, 0.001); // for sanity

		// AEE
		test_aee_cms_error_on_arrival(N, RowCounterNum << 1, RowsNum, Seed, data);
		test_aee_max_speed_cms_error_on_arrival(N, RowCounterNum << 1, RowsNum, Seed, data, 0.02, 0.001);

		// salsa baseline
		test_fg_salsa_baseline_cms_error_on_arrival(1, N, RowCounterNum << 5, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_error_on_arrival(2, N, RowCounterNum << 4, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_error_on_arrival(4, N, RowCounterNum << 3, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_error_on_arrival(8, N, RowCounterNum << 2, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_error_on_arrival(16, N, RowCounterNum << 1, RowsNum, Seed, data);

	}
	else if (Alg == 3)
	{
		// estimators - speed

		// our analytical solution (merge-downsample)
		test_analytics_salsa_cms_speed(N, RowCounterNum, RowsNum, Seed, data, 0.001, 0);
		test_analytics_salsa_cms_speed(N, RowCounterNum, RowsNum, Seed, data, 0.001, 1);
		test_analytics_salsa_cms_speed(N, RowCounterNum, RowsNum, Seed, data, 0.001, 2);
		test_analytics_salsa_cms_speed(N, RowCounterNum, RowsNum, Seed, data, 0.001, 4);
		test_analytics_salsa_cms_speed(N, RowCounterNum, RowsNum, Seed, data, 0.001, 8);
		test_analytics_salsa_cms_speed(N, RowCounterNum, RowsNum, Seed, data, 0.001, 10);

		// AEE
		test_aee_cms_speed(N, RowCounterNum >> 1, RowsNum, Seed, data);
		test_aee_max_speed_cms_speed(N, RowCounterNum >> 1, RowsNum, Seed, data, 0.02, 0.001);

		// salsa baseline
		test_maximum_salsa_baseline_cms_speed(N, RowCounterNum, RowsNum, Seed, data);
	}
	else if (Alg == 4)
	{
		test_maximum_salsa_baseline_cus_error_on_arrival(N, RowCounterNum << 2, RowsNum, Seed, data);
		test_baseline_cus_error_on_arrival(N, RowCounterNum, RowsNum, Seed, data);
	}
	else if (Alg == 5) 
	{
		test_maximum_salsa_baseline_cus_speed(N, RowCounterNum << 2, RowsNum, Seed, data);
		test_baseline_cus_speed(N, RowCounterNum, RowsNum, Seed, data);
	}
	else if (Alg == 6)
	{
		test_baseline_fg_salsa_cms_final_error(N, 32, RowCounterNum, RowsNum, Seed, data);
		test_baseline_fg_salsa_cms_final_error(N, 16, RowCounterNum << 1, RowsNum, Seed, data);
		test_baseline_fg_salsa_cms_final_error(N, 8, RowCounterNum << 2, RowsNum, Seed, data);
		test_baseline_fg_salsa_cms_final_error(N, 4, RowCounterNum << 3, RowsNum, Seed, data);
		test_baseline_fg_salsa_cms_final_error(N, 2, RowCounterNum << 4, RowsNum, Seed, data);
		test_baseline_fg_salsa_cms_final_error(N, 1, RowCounterNum << 5, RowsNum, Seed, data);

		test_analytics_fg_salsa_cms_final_error(N, 32, RowCounterNum, RowsNum, Seed, data, 0, 0, 0.001);
		test_analytics_fg_salsa_cms_final_error(N, 16, RowCounterNum << 1, RowsNum, Seed, data, 0, 0, 0.001);
		test_analytics_fg_salsa_cms_final_error(N, 8, RowCounterNum << 2, RowsNum, Seed, data, 0, 0, 0.001);
		test_analytics_fg_salsa_cms_final_error(N, 4, RowCounterNum << 3, RowsNum, Seed, data, 0, 0, 0.001);
		test_analytics_fg_salsa_cms_final_error(N, 2, RowCounterNum << 4, RowsNum, Seed, data, 0, 0, 0.001);
		test_analytics_fg_salsa_cms_final_error(N, 1, RowCounterNum << 5, RowsNum, Seed, data, 0, 0, 0.001);
	}
	else if (Alg == 7)
	{
		test_fg_salsa_baseline_cms_count_distinct(32, N, RowCounterNum, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_count_distinct(16, N, RowCounterNum << 1, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_count_distinct(8, N, RowCounterNum << 2, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_count_distinct(4, N, RowCounterNum << 3, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_count_distinct(2, N, RowCounterNum << 4, RowsNum, Seed, data);
		test_fg_salsa_baseline_cms_count_distinct(1, N, RowCounterNum << 5, RowsNum, Seed, data);
	}
	else if (Alg == 8)
	{
		test_fg_salsa_count_sketch_turnstile(32, N, RowCounterNum, RowsNum, Seed, data);
		test_fg_salsa_count_sketch_turnstile(16, N, RowCounterNum << 1, RowsNum, Seed, data);
		test_fg_salsa_count_sketch_turnstile(8, N, RowCounterNum << 2, RowsNum, Seed, data);
		test_fg_salsa_count_sketch_turnstile(4, N, RowCounterNum << 3, RowsNum, Seed, data);
		test_fg_salsa_count_sketch_turnstile(2, N, RowCounterNum << 4, RowsNum, Seed, data);
		test_fg_salsa_count_sketch_turnstile(1, N, RowCounterNum << 5, RowsNum, Seed, data);
	}
	else if (Alg == 9)
	{
		test_fg_salsa_count_sketch_all_except_turnstile(32, N, RowCounterNum, 5, Seed, data, 0, 1024, false);
		test_fg_salsa_count_sketch_all_except_turnstile(16, N, RowCounterNum << 1, 5, Seed, data, 0, 1024, false);
		test_fg_salsa_count_sketch_all_except_turnstile(8, N, RowCounterNum << 2, 5, Seed, data, 0, 1024, false);
		test_fg_salsa_count_sketch_all_except_turnstile(4, N, RowCounterNum << 3, 5, Seed, data, 0, 1024, false);
		test_fg_salsa_count_sketch_all_except_turnstile(2, N, RowCounterNum << 4, 5, Seed, data, 0, 1024, false);
		test_fg_salsa_count_sketch_all_except_turnstile(1, N, RowCounterNum << 5, 5, Seed, data, 0, 1024, false);
	}
	else if (Alg == 10)
	{	
		test_univmon_final_all(32, N, RowCounterNum, RowsNum, Seed, data, 0, (128 * RowCounterNum) / 4096, false);
		test_univmon_final_all(16, N, RowCounterNum << 1, RowsNum, Seed, data, 0, (128 * RowCounterNum) / 4096, false);
		test_univmon_final_all(8, N, RowCounterNum << 2, RowsNum, Seed, data, 0, (128 * RowCounterNum) / 4096, false);
		test_univmon_final_all(4, N, RowCounterNum << 3, RowsNum, Seed, data, 0, (128 * RowCounterNum) / 4096, false);
		test_univmon_final_all(2, N, RowCounterNum << 4, RowsNum, Seed, data, 0, (128 * RowCounterNum) / 4096, false);
	}
	else if (Alg == 11)
	{
		test_analytics_fg_salsa_cms_error_on_arrival(N, 1, RowCounterNum << 5, RowsNum, Seed, data, 0, 1, 0.001); 
		test_analytics_fg_salsa_cms_error_on_arrival(N, 2, RowCounterNum << 4, RowsNum, Seed, data, 0, 1, 0.001); 
		test_analytics_fg_salsa_cms_error_on_arrival(N, 4, RowCounterNum << 3, RowsNum, Seed, data, 0, 1, 0.001); 
		test_analytics_fg_salsa_cms_error_on_arrival(N, 8, RowCounterNum << 2, RowsNum, Seed, data, 0, 1, 0.001); 
		test_analytics_fg_salsa_cms_error_on_arrival(N, 16, RowCounterNum << 1, RowsNum, Seed, data, 0, 1, 0.001);
	}
	else
	{
		system("pause");
	}

	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////	
	///////////////////////////////////////////////////////////////////////
	///////////////////////////////////////////////////////////////////////

	delete[] data;
	delete[] data_weights;

	return 0;

}

