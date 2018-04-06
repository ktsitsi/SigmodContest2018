#ifndef SCHEDULER_H_
#define SCHEDULER_H_

#include "queue.h"
#include "job.h"
#include "worker.h"
#include <pthread.h>
#include <cstdlib>
#include <stdint.h>

#define RESULTS_NUM 4000
#define POOL_SIZE 1000 

// Maximum number of jobs in group
#define GROUPJ_NUM 1

typedef Job *Job_ptr;

class JobScheduler;

struct Package                                                                                                        
{
    uint32_t id;
    JobScheduler *job_scheduler;
};

class JobScheduler {
    private:
        pthread_t *tids_;    // Ids of threads
        Package *pkg_;       // Packages of threads

        int32_t finished_;      // Number of jobs finished yet
        bool imported_;         // Wheterh batch loading finished
        bool done_;             // Whether workload finished

        int32_t threshold_;     // Batch length
        Queue<Job_ptr> queue_;  // Queue of jobs (pointers)

        pthread_mutex_t mutex_;  // Controls queue
        pthread_mutex_t burst_;  // Controls finished counter
        pthread_cond_t cond_nonempty_;   // Signaled to wake up a thread
        pthread_cond_t cond_empty_;      // Signaled to wake up scheduler checking

        unsigned int threads_;   // Number of threads
    public:
        JobScheduler(unsigned int tnum);
        void submit_job (Job *job);
        void execute_all_jobs ();   // Execution of jobs begins
        void wait_all_tasks_finish();   // Wait all submitted jobs to finish
        void register_job ();       // Worker reports finished job
        unsigned int get_n_jobs (Job_ptr** elements);   // Worker gets jobs if available
        ~JobScheduler();
};

#endif
