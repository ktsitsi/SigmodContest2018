#ifndef QUEUE_H
#define QUEUE_H

#include <cstdlib>

#define MIN(a,b) (((a)<(b)) ? (a) : (b))

// Array-based waiting queue
// Initialised with a max size and resizable when necessary
// Example use:
// Queue<int> q(10);
// q.Enqueue(1);
// q.Pop();
template <typename E>
class Queue
{
    public:
        Queue(size_t s)
        {
            front_ = 0;
            rear_ = 0;
            max_size_ = s;
            len_ = 0;
            buffer_ = (E*)malloc(s*sizeof(E));
        };

        ~Queue() 
        {
            free(buffer_);
            buffer_ = NULL;

        };

        bool is_empty() const { return !len_; }

        size_t get_length() const { return len_; }

        const E& get_front() const { return buffer_[front_]; }

        unsigned int get_front_n(unsigned int num, E **dest) const
        {
            unsigned int available = MIN(num,len_);
            *dest = buffer_ + front_;
            return available;
        }

        void clear() { rear_ = front_ = len_ = 0; }

        void Enqueue(const E& element)
        {
            if (rear_ == max_size_){
                // queue at max capacity, allocate double space
                buffer_ = (E*)realloc(buffer_,(max_size_ <<= 1)*sizeof(E));
            }

            buffer_[rear_++]=element;
            ++len_;
        }


        void Pop() 
        {
            if (len_) {
                ++front_;
                --len_;
            };
        }

        void PopN(unsigned int num)
        {
            if (len_>= num) {
                front_ += num;
                len_ -= num;
            }
        }

    private:
        E *buffer_;
        size_t front_;  // index to pop from
        size_t rear_;   // index to enqueue
        size_t len_;
        size_t max_size_;
};

#endif
