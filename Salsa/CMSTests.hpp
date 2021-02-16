#ifndef CMS_TESTS
#define CMS_TESTS

/* Unweighted */

// cms
void test_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data);
void test_baseline_cms_speed(int N, int width, int height, int seed, const char* data);

// cus
void test_baseline_cus_error_on_arrival(int N, int width, int height, int seed, const char* data);
void test_baseline_cus_speed(int N, int width, int height, int seed, const char* data);

/* Weighted */
void test_weighted_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data, const uint16_t* data_weights);
void test_weighted_baseline_cms_speed(int N, int width, int height, int seed, const char* data, const uint16_t* data_weights);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

// cms
void test_salsa_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data);
void test_salsa_baseline_cms_speed(int N, int width, int height, int seed, const char* data);

void test_maximum_salsa_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data);
void test_maximum_salsa_baseline_cms_speed(int N, int width, int height, int seed, const char* data);

void test_fg_salsa_baseline_cms_error_on_arrival(int salsa_bits, int N, int width, int height, int seed, const char* data);
void test_fg_salsa_baseline_cms_speed(int salsa_bits, int N, int width, int height, int seed, const char* data);

void test_fg_salsa_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge = 3);
void test_fg_salsa_cms_speed(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge = 3);

void test_fg_salsa_split_cms_error_on_arrival(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge = 3);

void test_salsa_cms_error_on_arrival(int N, int width, int height, int seed, const char* data, int downsamplings_per_merge = 3);
void test_salsa_cms_speed(int N, int width, int height, int seed, const char* data, int downsamplings_per_merge = 3);

void test_analytics_salsa_cms_error_on_arrival(int N, int width, int height, int seed, const char* data, double delta, int downsamplings_before_merge);
void test_analytics_salsa_cms_speed(int N, int width, int height, int seed, const char* data, double delta, int downsamplings_before_merge);

void test_analytics_fg_salsa_cms_error_on_arrival(int N, int counterSize, int width, int height, int seed, const char* data, int downsamplings_before_merge, bool split_counters_flag, double delta);

void test_aee_cms_error_on_arrival(int N, int width, int height, int seed, const char* data);
void test_aee_cms_speed(int N, int width, int height, int seed, const char* data);

void test_aee_max_speed_cms_error_on_arrival(int N, int width, int height, int seed, const char* data, double epsilon_cc, double delta_cc);
void test_aee_max_speed_cms_speed(int N, int width, int height, int seed, const char* data, double epsilon_cc, double delta_cc);

// cus
void test_maximum_salsa_baseline_cus_error_on_arrival(int N, int width, int height, int seed, const char* data);
void test_maximum_salsa_baseline_cus_speed(int N, int width, int height, int seed, const char* data);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

// cms
void test_maximum_tango_baseline_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data);
void test_maximum_tango_baseline_cms_speed(int tango_bits, int N, int width, int height, int seed, const char* data);

void test_maximum_smartango_baseline_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data);
void test_tango_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge = 3);
void test_smartango_cms_error_on_arrival(int tango_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge = 3);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

// sanity
void test_sanity_cms(int N, int width, int height, int seed, const char* data);
void test_sanity_salsa_cms(int N, int width, int height, int seed, const char* data);
void test_sanity_cus(int N, int width, int height, int seed, const char* data);
void test_fg_salsa_split_counters_sanity_cms(int salsa_bits, int N, int width, int height, int seed, const char* data, int downsamplings_per_merge);

// count distinct with fg-salsa
void test_fg_salsa_baseline_cms_count_distinct(int salsa_bits, int N, int width, int height, int seed, const char* data);

// relative (final) error
void test_analytics_fg_salsa_cms_final_error(int N, int counterSize, int width, int height, int seed, const char* data, int downsamplings_before_merge, bool split_counters_flag, double delta);
void test_baseline_fg_salsa_cms_final_error(int N, int counterSize, int width, int height, int seed, const char* data);

///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////

void test_salsa_cms_final_error_histogram(int N, int width, int height, int seed, const char* data);

void test_baseline_cms_final_error_histogram(int N, int width, int height, int seed, const char* data);

#endif // !CMS_TESTS