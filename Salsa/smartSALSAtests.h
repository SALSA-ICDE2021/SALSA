#pragma once

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
#include "smartSALSA.hpp"
#include "SalsaCMSBaseline.hpp"

#ifndef SMART_SALSA_TESTS_H
#define SMART_SALSA_TESTS_H

void test_SmartEncodingSALSA_baseline_cms_sanity(int N, int width, int height, int seed, const char* data);
void test_SmartEncodingSALSA_baseline_cms_error_on_arrival(int N, int width, int height, int seed, const char* data);
void test_SmartEncodingSALSA_baseline_cms_speed(int N, int width, int height, int seed, const char* data);
void test_FastEncodingSALSA_baseline_cms_speed(int N, int width, int height, int seed, const char* data);

void test_StupidEncodingSALSA_baseline_cms_speed(int N, int width, int height, int seed, const char* data);
void test_StupiderEncodingSALSA_baseline_cms_speed(int N, int width, int height, int seed, const char* data);

#endif // SMART_SALSA_TESTS_H