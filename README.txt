SIGMOD PROGRAMMING COMPETITION

TEAM GRANMA
Ecole Polytechnique Federale de Lausanne

Ergys Dona
Panagiotis Sioulas
Konstantinos Tsitsimpikos

In this project, we are given batches of queries (which consist of a series of joins and selections) and want to compute the results required as fast as possible. The next batch is only given after the results of the previous one are returned.

First we are given 1 second to prepare for query execution. During this time period, we scan the data to compute some basic statistics (min and max value) and sample the data in order to estimate distinct values. This is run in a multithreaded manner, otherwise it wouldn't be computed in a timely manner for large datasets.

After that, there is the loop of processing queries. Here each query is parsed, rewritten to another equivalent so that the requested checksums will be computed once (some queries request checksums of the same column multiple times or on columns that have been joined). This is done with the QueryInfo::rewriteQuery function. After that, the rewritten query is passed to the optimizer.

The optimizer uses the statistics gathered in the beginning along with Selinger's algorithm to compute the least expensive join order.

Using the optimized query, we construct a query plan which is essentially a tree of operators. We also use passes in order to modify the plan (for example pushing down projections or materializer utilities). Once we have the final query, we can execute the plan by calling next for the root.

The executor works in a vector at a time mode. Each operator calls Next on its children to fetch the next vector of 1024 elements. However, we can fetch vectors of multiple columns. Projection and Scan are zero-copy in that they just manipulate pointers to underlying objects. Selections fill a vector with qualifying data and forward it once the vector is full. A Self-Join operator is used to close a cycle of joins by esssentially doing an equality comparison between two columns (works like Selections). For Joins, we use Radix Hash Join with two passes. Again, the operator does build and probe for each copartition pair and generates results until it fills the buffer. Selections, Joins and Self-Joins need to retain state from their previous next calls so that they are re-entrant and might have to append result that couldn't be attached to the previous vector. Finally, checksum just accumulates the columns involved.

The query plans can be passed to a multithreading pool as jobs to introduce parallelism. But our implementation stays single threaded because we didnt have the time to implement a good concurrency pan for our implementatio. The main idea was to run the operators in different threads.

We didn't have time for custom memory allocaion, numa-awareness, cache-awareness, finer-grained concurrency and work sharing.

In the implementation we have also ready a late materialization module that was hooked in the scheme but didnt help so much because of the its memory copies. We could refine this even better if we had the time , because our idea was to import this module only in cases that it would be helpful in the query plan without so big impact in memory accesses. 

