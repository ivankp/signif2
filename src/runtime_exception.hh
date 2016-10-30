#ifndef IVANP_RUNTIME_EXCEPTION
#define IVANP_RUNTIME_EXCEPTION

#include <string>
#include <sstream>
#include <exception>

class rte : public std::exception {
  std::string str;
public:
  rte() noexcept { }

  template <typename... Args>
  rte(Args&& ...args) noexcept {
    std::stringstream ss;
    (ss << ... << args);
    str = ss.str();
  }

  rte(const char* cstr) noexcept : str(cstr) { }
  rte& operator= (const char* cstr) noexcept {
    str = cstr;
    return *this;
  }

  rte(const std::string& s) noexcept : str(s) { }
  rte(std::string&& s) noexcept : str(std::move(s)) { }
  rte& operator= (const std::string& s) noexcept {
    str = s;
    return *this;
  }
  rte& operator= (std::string&& s) noexcept {
    str = std::move(s);
    return *this;
  }

  rte(const rte& e) noexcept : str(e.str) { }
  rte(rte&& e) noexcept : str(std::move(e.str)) { }
  rte& operator= (const rte& e) noexcept {
    str = e.str;
    return *this;
  }
  rte& operator= (rte&& e) noexcept {
    str = std::move(e.str);
    return *this;
  }

  rte(const std::exception& e) noexcept : str(e.what()) { }
  rte& operator= (const std::exception& e) noexcept {
    str = e.what();
    return *this;
  }

  virtual const char* what() const noexcept {
    return str.c_str();
  }
};

#endif
