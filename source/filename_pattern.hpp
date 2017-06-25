/*! \file filename_pattern.hpp
  \brief Abstract class for file name generator
  \author Almog Yalinewich
*/

#ifndef FILENAME_PATTERN_HPP
#define FILENAME_PATTERN_HPP 1

#include <string>

using std::string;

//! \brief Abstract class for file name pattern
class FileNamePattern
{
public:

  /*! \brief Converts an index to a file name
    \param index Index
    \return File name
   */
  virtual string operator()(int index) = 0;
  
  virtual ~FileNamePattern(void);
};

#endif // FILENAME_PATTERN_HPP
