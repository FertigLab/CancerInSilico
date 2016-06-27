#ifndef EXCEPTION_HANDLING_HPP
#define EXCEPTION_HANDLING_HPP

#ifndef _WIN32

#include <execinfo.h>

inline void RCPP_STOP_TRACE(std::string message) {
    int nptrs;
    void *buffer[100];
    nptrs = backtrace(buffer, 100);
    backtrace_symbols_fd(buffer, nptrs, STDOUT_FILENO);
    Rcpp::stop(message);
}

#else

inline void RCPP_STOP_TRACE(std::string message) {
    Rcpp::stop(message);
}

#endif

#endif
