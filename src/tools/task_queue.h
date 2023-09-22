

#ifndef TASK_QUEUE_H_
#define TASK_QUEUE_H_

template<typename _t, typename _callback>
struct Task_queue
{

	Task_queue(unsigned limit, _callback &callback):
		queue_ (limit),
		state_ (limit),
		head_ (0),
		tail_ (0),
		limit_ (limit),
		head_idx_ (0),
		queued_ (0),
		at_end_ (false),
		callback_ (callback)
	{ }

	volatile bool waiting() const
	{ return tail_ - head_ >= limit_; }

	template<typename _init>
	volatile bool get(size_t &n, _t*& res, _init &init)
	{
		{
			tthread::lock_guard<tthread::mutex> lock (mtx_);
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " waiting=" << waiting() << " head=" << head_ << " tail=" << tail_ << endl;
#endif
			while(waiting() && !at_end_ && !exception_state())
				cond_.wait(mtx_);
			if(exception_state() || at_end_) {
#ifdef ENABLE_LOGGING
				log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " quit" << endl;
#endif
				return false;
			}
			n = tail_++;
			res = &slot(n);
			if(!init())
				at_end_ = true;
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue get() thread=" << omp_get_thread_num() << " n=" << n << endl;
#endif
		}
		if(at_end_)
			cond_.notify_all();
		return true;
	}

	void wake_all()
	{ cond_.notify_all(); }

	volatile void push(unsigned n)
	{
		mtx_.lock();
		if(n == head_) {
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue flush() thread=" << omp_get_thread_num() << " n=" << n << endl;
#endif
			mtx_.unlock();
			flush();
		} else {
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue push() thread=" << omp_get_thread_num() << "n=" << n << " head=" << head_ << endl;
#endif
			state_[idx(n)] = true;
			++queued_;
			mtx_.unlock();
		}
	}

	volatile unsigned flush()
	{
		callback_(slot(head_));
		bool next=false;
		{
			tthread::lock_guard<tthread::mutex> lock(mtx_);
#ifdef ENABLE_LOGGING
			log_stream << "Task_queue flush() thread=" << omp_get_thread_num() << " head=" << head_ << " waiting=" << notify << endl;
#endif
			state_[idx(head_)] = false;
			++head_;
			head_idx_ = (head_idx_+1)%limit_;
			if(state_[idx(head_)]) {
				next = true;
#ifdef ENABLE_LOGGING
				log_stream << "Task_queue n=" << queued_ << endl;
#endif
				--queued_;
			}
		}
		cond_.notify_one();
		if(next)
			return flush()+1;
		return 1;
	}

private:

	unsigned idx(unsigned n) const
	{ return (head_idx_ + n - head_)%limit_; }

	_t& slot(unsigned n)
	{ return queue_[idx(n)]; }

	vector<_t> queue_;
	vector<bool> state_;
	tthread::mutex mtx_;
	tthread::condition_variable cond_;
	volatile unsigned head_, tail_, limit_, head_idx_, queued_;
	bool at_end_;
	_callback &callback_;

};

#endif
