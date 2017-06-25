#ifndef UNIVERSAL_ERROR
#define UNIVERSAL_ERROR 1

#include <vector>
#include <string>

using std::vector;
using std::string;

/*! \brief A universal method to convey exception data
  \details The class contains an error message and entries that each layer of the code
  can add to, but not to delete.
 */
class UniversalError
{
public:

  /*! \brief Class constructor
    \param err_msg Error message
   */
  UniversalError(string const& err_msg);
  
  /*! \brief Appends string to error message
    \param em_add String to attach
   */
  void appendToMessage(string const& em_add);
  
  /*! \brief appends an entry
    \param field Parameter name
    \param val Parameter value
   */
  void operator()(string const& field, double val);
  
  /*! \brief Returns the error message
   */
  string const& getErrorMessage(void) const;

  /*! \brief Access to error data
    \return Reference to error data
   */
  const vector<std::pair<string, double> >& getData(void) const;

private:
  
  //! \brief Description of the exception
  string err_msg_;

  //! \brief Relevant data
  vector<std::pair<string, double> > data_;
};

void report_error(const UniversalError& eo, std::ostream& s);

#endif // UNIVERSAL_ERROR
