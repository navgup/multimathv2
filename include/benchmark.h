
#ifndef MATHOMP_BENCHMARK_H
#define MATHOMP_BENCHMARK_H

#include <vector>
#include <chrono>
#include <atomic>
#include <cstring>

namespace HuBench {
    class HuTime {
        std::chrono::high_resolution_clock _clock;

    public:
        static inline uint64_t millisec() {
            return std::chrono::duration_cast<std::chrono::milliseconds>
                    (std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        }

        static inline uint64_t microsec() {
            return std::chrono::duration_cast<std::chrono::microseconds>
                    (std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        }

        static inline uint64_t nanosec() {
            return std::chrono::duration_cast<std::chrono::nanoseconds>
                    (std::chrono::high_resolution_clock::now().time_since_epoch()).count();
        }
    };


    class Benchmark {
#define BENCH_VECTOR_CAPACITY 250
#define T_ENTRY_LABEL_SIZE 50
        typedef struct TimeEntry {
            uint64_t time = 0U;
            char label[T_ENTRY_LABEL_SIZE]{};

            TimeEntry() = default;
        } TimeEntry;
        using darray = std::vector<TimeEntry>;
    public:
        Benchmark() : _ActiveTime(0) {
            _Timings.reserve(BENCH_VECTOR_CAPACITY);
        }

        inline void startTiming(const char *msg_fmt) {
            if (!_IsRunning) {
                _IsRunning.store(true);
                strncpy(_Timings[_CurIndex].label, msg_fmt, T_ENTRY_LABEL_SIZE);

                _ActiveTime = HuTime::nanosec();
            }
        }

        inline void stopTiming() {
            uint64_t time = HuTime::nanosec();
            if (_IsRunning) {
                _Timings[_CurIndex].time = time - _ActiveTime;

                _IsRunning.store(_CurIndex >= BENCH_VECTOR_CAPACITY - 1);
                _CurIndex += _CurIndex < BENCH_VECTOR_CAPACITY - 1;
            }
        }

        void printAll() {
            for (int i = 0; i <= _CurIndex; i++) {
                printf(_Timings[i].label, _Timings[i].time);
            }
            reset();
        }

        void reset() {
            size = _Timings.size();
            _CurIndex = 0;
            _Timings.clear();
        }

        void printAt(int index) const {
            if (index <= _CurIndex) {
                printf(_Timings[index].label, _Timings[index].time);
            } else {
                printf("Access out of bounds! No entry for this index available.\n");
            }
        }

    private:
        int size = BENCH_VECTOR_CAPACITY;
        int _CurIndex = 0;
        darray _Timings;

        std::atomic_bool _IsRunning{false};
        uint64_t _ActiveTime;
    };

}; // End of namespace HuBench

#endif //MATHOMP_BENCHMARK_H
