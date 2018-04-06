#include "include/concurrency.h"

using namespace std;

namespace Granma {

const unsigned Concurrency::NUM_OF_THREADS = 16u;

queue<JobInfo> Concurrency::query_queue_;
mutex Concurrency::mtx_;
condition_variable Concurrency::cond_;
bool Concurrency::finished_;

vector<thread> Concurrency::threads_;

void Concurrency::ProduceQueries(ifstream &work_file, Joiner &joiner) {
    QueryInfo query_info, query_info_rewrite, query_info_optimize;

    // results will be populated here
    vector<future<string>> results;

    // create and start threads
    for (unsigned i = 0; i != NUM_OF_THREADS; ++i) {
        thread thr(Concurrency::ConsumeQueries);
        thr.detach();
        threads_.emplace_back(std::move(thr));
    }

    string line;
    while (getline(work_file, line)) {
        if (line == "F") { // End of a batch, print results
            for (auto &futurama : results) {
                futurama.wait();
                cout << futurama.get() << endl; // print result
            }
            results.clear(); // clear results to prepare for new batch
            continue;
        }
        query_info.parseQuery(line);    // parse query
        query_info_rewrite.rewriteQuery(query_info);
        query_info_optimize = joiner.optimizing(query_info_rewrite);
        promise<string> result_promise; // create new promise
        // append pending result to 'results'
        results.emplace_back(result_promise.get_future());
        {
            // take queue lock
            lock_guard<mutex> lock(mtx_);
            // add new job to queue
            query_queue_.emplace(joiner, query_info_optimize, std::move(result_promise));
            // mtx_ is unlocked once lock goes out of scope
        }
        // notify threads after unlocking mtx_, to minimize lock contention
        cond_.notify_all();
        // make QueryInfo object ready for next query
        query_info.clear();
    }
    {
        lock_guard<mutex> lock(mtx_);
        finished_ = true;
        // mtx_ is unlocked once lock goes out of scope
    }
    // notify threads
    cond_.notify_all();
}

void Concurrency::ConsumeQueries() {
    unique_lock<mutex> lock(mtx_, std::defer_lock);
    while (true) {
        // obtain the mutex lock via the unique_lock object
        lock.lock();
        // wait until queue is not empty
        cond_.wait(lock, [&] { return !query_queue_.empty() || finished_; });
        if (finished_ && query_queue_.empty())
            break;
        // get next pending job
        JobInfo job = std::move(query_queue_.front());
        query_queue_.pop();
        // unlock mutex, to allow other threads to obtain queries
        lock.unlock();
        // extract job info
        Joiner &joiner = std::get<0>(job);
        QueryInfo &qinfo = std::get<1>(job);
        promise<string> &prom = std::get<2>(job);
        // run join
        //QueryInfo j;
        //j.rewriteQuery(qinfo);
        //joiner.join(std::move(joiner.optimizing(j)));
        joiner.join(std::move(qinfo));
        prom.set_value(joiner.result());
    }
}

} // namespace Granma
