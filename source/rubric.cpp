#include "rubric.hpp"
#include "utilities.hpp"

Rubric::Rubric(const string& prefix,
	       const string& postfix):
  prefix_(prefix), postfix_(postfix) {}

string Rubric::operator()(int index)
{
  return prefix_+int2str(index)+postfix_;
}
