#include <iostream>
#include <algorithm>
#include "universal_error.hpp"

using std::pair;
using std::for_each;

UniversalError::UniversalError(string const& err_msg):
  err_msg_(err_msg),
  data_() {}

void UniversalError::appendToMessage(string const& em_add)
{
  err_msg_ += em_add;
}

void UniversalError::operator()(string const& field, double val)
{
  data_.push_back(std::pair<string,double>(field,val));
}

string const& UniversalError::getErrorMessage(void) const
{
  return err_msg_;
}

const vector<std::pair<string,double> >& UniversalError::getData(void) const
{
  return data_;
}

void report_error(const UniversalError& eo,
		   std::ostream& s)
{
  s << eo.getErrorMessage() << "\n";
  for_each(eo.getData().begin(),
	   eo.getData().end(),
	   [&s](const pair<string, double>& x)
	   {s << x.first << ": " << x.second << "\n";});
}
