#ifndef VECTOR_INITIALIZER_HPP
#define VECTOR_INITIALIZER_HPP 1

//! \brief Class for initialising vectors
template<class T> class VectorInitializer
{
public:

  /*! \brief Class constructor
    \param t First member
   */
  VectorInitializer(T const& t):
    list_(1,t) {}

  /*! \brief Appends member to initializer
    \param t Appendage
    \return Self reference
   */
  VectorInitializer& operator()(T const& t)
  {
    list_.push_back(t);
    return *this;
  }

  /*! \brief Returns complete vector
    \return stl vector
   */
  vector<T> const& operator()(void)
  {
    return list_;
  }

private:

  //! \brief stl vector
  vector<T> list_;
};

#endif // VECTOR_INITIALIZER_HPP
