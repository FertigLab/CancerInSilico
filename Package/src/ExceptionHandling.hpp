#include <execinfo.h> 

#ifndef EXCEPTION_HANDLING_HPP
#define EXCEPTION_HANDLING_HPP

inline void RCPP_STOP_TRACE(std::string message) {
  
  int nptrs;
  void *buffer[100];
  nptrs = backtrace(buffer, 100);
  backtrace_symbols_fd(buffer, nptrs, STDOUT_FILENO);

  Rcpp::stop(message);

}

#endif
