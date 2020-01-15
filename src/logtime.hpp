#ifndef __LOGTIME_H__
#define __LOGTIME_H__

#include <ctime>
#include <ostream>


struct logtime { };

inline std::ostream& operator<<(std::ostream& os, const logtime& l) {
	time_t CurrentTime;
	struct tm LocalTime;

	time(&CurrentTime);
	if(localtime_r(&CurrentTime, &LocalTime)) {
	  char buf[128];
	  if(strftime(buf, sizeof(buf), "[%c]", &LocalTime))
	    return os << buf;
	}

	// If get there, something failed
	return os << "[???]";
}
#endif /* __LOGTIME_H__ */
