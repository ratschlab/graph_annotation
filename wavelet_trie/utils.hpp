#ifndef __UTILS_HPP__
#define __UTILS_HPP__

#include <cstdint>
#include <string>
#include <vector>
#include <deque>
#include <thread>
#include <mutex>
#include <vector>
#include <condition_variable>
#include <future>
#include <functional>
#include <stdexcept>
#include <queue>
#include <atomic>
#include <cassert>

namespace utils {

    /** This returns the currently used memory by the process.
     *
     * The code was copied and has been modified from:
     * https://github.com/progschj/ThreadPool/blob/master/ThreadPool.h
     */
    class ThreadPool {
      public:
        ThreadPool(size_t num_workers);

        template <class F, typename... Args>
        auto enqueue(F&& f, Args&&... args) -> std::future<decltype(f(args...))> {
            using return_type = typename std::result_of<F(Args...)>::type;
            auto task = std::make_shared<std::packaged_task<return_type()>>(
                std::bind(std::forward<F>(f), std::forward<Args>(args)...)
            );

            if (!workers.size()) {
                (*task)();
                return task->get_future();
            } else {
                std::unique_lock<std::mutex> lock(queue_mutex);
                tasks.emplace([task](){ (*task)(); });
            }
            condition.notify_one();

            return task->get_future();
        }

        void join();

        ~ThreadPool();

      private:
        void initialize(size_t num_threads);

        std::vector<std::thread> workers;
        std::queue<std::function<void()>> tasks;

        std::mutex queue_mutex;
        std::condition_variable condition;

        bool joining_;
        bool stop_;
    };

} // namespace utils

#endif // __UTILS_HPP__
