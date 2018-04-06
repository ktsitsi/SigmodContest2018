#include "include/scheduler.h"

JobScheduler::JobScheduler(unsigned int tnum) : queue_(POOL_SIZE), threads_(tnum) {
    finished_ = 0;
    imported_ = false;
    done_ = false;

    // Initialise mutexes and condition variables
    pthread_mutex_init(&mutex_, NULL);
    pthread_mutex_init(&burst_, NULL);
    pthread_cond_init(&cond_nonempty_, NULL);
    pthread_cond_init(&cond_empty_, NULL);

    tids_ = (pthread_t*)malloc(threads_*sizeof(pthread_t));
    pkg_ = (Package*)malloc(threads_*sizeof(Package));

    // Create packages (arguments) and threads
    for (int i = 0; i < threads_; i++) {

        pkg_[i].job_scheduler = this; 
        pkg_[i].id = i;

        pthread_create (&tids_[i], NULL, worker, &pkg_[i]);
    }
}

void JobScheduler::submit_job (Job *job) 
{
    // Safe enqueue
    pthread_mutex_lock(&mutex_);

    queue_.Enqueue(job);

    pthread_mutex_unlock(&mutex_);
}


void JobScheduler::execute_all_jobs () 
{
    pthread_mutex_lock(&mutex_);

    // Finished loading
    imported_ = true;
    threshold_ = queue_.get_length();

    pthread_mutex_unlock(&mutex_);

    // Wake up an available thread
    pthread_cond_signal(&cond_nonempty_);

}

void JobScheduler::wait_all_tasks_finish() 
{
    pthread_mutex_lock(&burst_);

    // Blocked until every task finished
    while (finished_ < threshold_) {
        pthread_cond_wait(&cond_empty_, &burst_);
    }

    pthread_mutex_unlock(&burst_);

    finished_ = 0;
    imported_=false;

    queue_.clear();
}

void JobScheduler::register_job () 
{
    pthread_mutex_lock(&burst_);

    int32_t cur_finished = ++finished_;

    pthread_mutex_unlock(&burst_);

    // Wake up scheduler if everyone finished
    if(cur_finished == threshold_) pthread_cond_signal(&cond_empty_);
}

// Attempt to get a group of jobs
// If workload finished or jobs are available returns 0 or the number respectively
// else it is blocked waiting for a job
unsigned int JobScheduler::get_n_jobs (Job_ptr** elements) {
    unsigned int group_num;

    pthread_mutex_lock(&mutex_);

    // Block until workload done or jobs available
    while ((queue_.is_empty() || !imported_) && !done_) {
        pthread_cond_wait(&cond_nonempty_, &mutex_);
    }

    // Unblocked because workload finished
    if (done_) {
        group_num = 0;
        pthread_mutex_unlock(&mutex_);
    } else {
        group_num = queue_.get_front_n(GROUPJ_NUM,elements); 
        queue_.PopN(group_num);

        bool more_allowed = !queue_.is_empty();

        pthread_mutex_unlock(&mutex_);
        // Wake up the next available thread
        if (more_allowed) {
            pthread_cond_signal(&cond_nonempty_);
        }
    }
    return group_num;
}

JobScheduler::~JobScheduler() {
    pthread_mutex_lock(&mutex_);

    done_ = true;   // Workload done

    pthread_mutex_unlock(&mutex_);

    // Wake up threads to finish
    pthread_cond_broadcast(&cond_nonempty_);

    void* ptr;

    for (int i = 0; i < threads_; i++)
        pthread_join (tids_[i], &ptr);

    pthread_mutex_destroy(&mutex_);
    pthread_mutex_destroy(&burst_);
    pthread_cond_destroy(&cond_nonempty_);
    pthread_cond_destroy(&cond_empty_);

    free(tids_);
    free(pkg_);
}

