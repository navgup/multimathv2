#include "mathomp.h"
#include "/Users/navgup/repos/benchmark/include/benchmark/benchmark.h"
#include "/Users/navgup/repos/benchmark/include/benchmark/export.h"
#include <cmath>

const static int size = 50;
double TestingData[size] = {};

static void BM_StdSin(benchmark::State& state) {
    srand48(8597243);
    for(int i=0; i<size; i++) {
        TestingData[i] = drand48();
    }

    for (auto _ : state) {
        //TestingData[state.range(0)] = sin(TestingData[state.range(0)]);
        benchmark::DoNotOptimize(sin(TestingData[state.range(0)]));
    }
}

static void BM_InternetSinus(benchmark::State& state) {
    srand48(91723412);
    for(int i=0; i<size; i++) {
        TestingData[i] = drand48();
    }

    for (auto _ : state) {
        //TestingData[state.range(0)] = sinusInternet(TestingData[state.range(0)]);
        benchmark::DoNotOptimize(sinusInternet(TestingData[state.range(0)]));

    }
}

static void BM_CustomSinus(benchmark::State& state) {
    srand48(1389434);
    for(int i=0; i<size; i++) {
        TestingData[i] = drand48();
    }

    for (auto _ : state) {
        //TestingData[state.range(0)] = sinus(TestingData[state.range(0)]);
        benchmark::DoNotOptimize(sinus(TestingData[state.range(0)]));
    }
}

BENCHMARK(BM_StdSin)->Arg(15)->Arg(30)->Arg(50);

BENCHMARK(BM_InternetSinus)->Arg(15)->Arg(30)->Arg(50);

BENCHMARK(BM_CustomSinus)->Arg(15)->Arg(30)->Arg(50);

BENCHMARK_MAIN();
