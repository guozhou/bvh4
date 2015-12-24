#pragma once
class Util
{
public:
	Util();
	~Util();
	// http://stackoverflow.com/a/21995693/695128
	template<typename TimeT = std::chrono::milliseconds>
	struct measure
	{
		template<typename F, typename ...Args>
		static typename TimeT::rep execution(F func, Args&&... args)
		{
			auto start = std::chrono::system_clock::now();
			func(std::forward<Args>(args)...);
			auto duration = std::chrono::duration_cast<TimeT>
				(std::chrono::system_clock::now() - start);
			return duration.count();
		}
	};
};

