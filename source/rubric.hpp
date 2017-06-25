/*! 
  \brief A rubric based 
  \author Almog Yalinewich
*/

#ifndef RUBRIC_HPP
#define RUBRIC_HPP 1

#include "filename_pattern.hpp"

class Rubric: public FileNamePattern
{
public:

  /*! \brief Class constructor
    \param prefix Prefix
    \param postfix Postfix
   */
  Rubric(const string& prefix,
	 const string& postfix);

  string operator()(int index);

private:

  //! \brief Prefix
  const string prefix_;

  //! \brief Postfix
  const string postfix_;
};

#endif // RUBRIC_HPP
