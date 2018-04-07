#include "include/worker.h"
#include "include/scheduler.h"
#include <pthread.h>
#include <sched.h>

// Thread work 
void* worker(void* arg) 
{
    Package* pkg = (Package*) arg;
    JobScheduler* job_scheduler = pkg->job_scheduler;

    Job_ptr *jobs;

    unsigned int order = 0;

    //printf ("%d\n", pkg->id);
    cpu_set_t *cpusetp;
    size_t num_cpus = 40;

    cpusetp = CPU_ALLOC(num_cpus);
    if (cpusetp == NULL)
        printf("Houston, we got a problem\n");

    size_t sizep = CPU_ALLOC_SIZE(num_cpus);

    CPU_ZERO_S(sizep, cpusetp);
    CPU_SET_S(((pkg->id) * 2 + 1) % num_cpus, sizep, cpusetp);

    //printf ("%d\n", CPU_COUNT_S(sizep, cpusetp));

    pthread_setaffinity_np(pthread_self(), sizep, cpusetp);

    while (1) {
        // Attempt to get jobs available
        unsigned int group_num = job_scheduler->get_n_jobs(&jobs);
        // Failed to get any jobs and returned : workload finished
        if (!group_num) break;

        for (unsigned int i=0;i<group_num;++i) {
            // Run job task : virtual function so the task of the respective subclass is called
            jobs[i]->RunTask(pkg->id,++order);
            delete jobs[i];
            // Inform scheduler that task finished
            job_scheduler->register_job();
        }
    }

    CPU_FREE(cpusetp);
}
