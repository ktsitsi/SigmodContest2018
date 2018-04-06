#pragma once
#include "joiner.h"
#include <fstream>
#include <future>
#include <queue>
#include <string>
#include <tuple>

// concurrency and synchronisation
#include <condition_variable>
#include <mutex>
//#include <shared_mutex>
#include <thread>

namespace Granma {

using JobInfo = std::tuple<Joiner, QueryInfo, std::promise<std::string>>;

class Concurrency {
  public:
    Concurrency() = default;
    ~Concurrency() = default;

    static void ProduceQueries(std::ifstream &work_file, Joiner &joiner);
    static void ConsumeQueries();

    static const unsigned NUM_OF_THREADS;

  private:
    static std::queue<JobInfo> query_queue_;
    static std::mutex mtx_;
    static std::condition_variable cond_;
    static bool finished_;

    static std::vector<std::thread> threads_;
};

} // namespace Granma
