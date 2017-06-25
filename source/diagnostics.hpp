/*! \file diagnostics.hpp
  \brief Diagnostic functions for srhd_sim
  \author Almog Yalinewich
*/

#ifndef DIAGNOSTICS_HPP
#define DIAGNOSTICS_HPP 1

#include <fstream>
#include "hydrodynamic_variables.hpp"
#include "srhd_sim.hpp"

/*! \brief Approximate comparison
  \details Verifies that two values are very close to one another
  \param v1 First value
  \param v2 Second value
  \param thres Threshold
  \return True if the difference between v1 and v2 is below thres, false otherwise
 */
bool ApproxCompare(double v1, double v2, double thres);

/*! \brief Verifies that the primitive variables are consistent with the conserved variable, up to a specified accuracy
  \param p Primitive variables
  \param c Conserved variables
  \param eos Equation of state
  \param thres Threshold
  \return True if both parameters agree, false otherwise
 */
bool ConservedPrimitiveConsistency(Primitive const& p,
				   NewConserved const& c,
				   const EquationOfState& eos,
				   double thres);

/*! \brief Verifies that the primitive variables in all the cells are consistent with the conserved variables
  \param sim Hydrodynamic simulation
  \param thres Threshold
  \return True if all cells are consistent, false otherwise
 */
bool ConservedPrimitiveConsistency(SRHDSimulation const& sim, 
				   double thres);

/*! \brief Returns the total energy
  \param sim Hydrodynamic simulation
  \return Total energy
 */
double TotalEnergy(SRHDSimulation const& sim);

/*! \brief Returns the total momentum
  \param sim Hydrodynamic simulation
  \return Total momentum
 */
double TotalMomentum(SRHDSimulation const& sim);

/*! \brief Checks that the vertices are in increasing order
  \param sim Reference to simulation
  \return True if they are in increasing order, false otherwise
 */
bool VerticesIncreasingOrder(SRHDSimulation const& sim);

/*! \brief Retrieves a property of a hydrodynamic cell
  \param sim Hydrodynamic simulation
  \param i Cell index
  \param property Name of propery
  \return Value of property
 */
double get_cell_property(SRHDSimulation const& sim,
			 size_t i, string const& property);

/*! \brief Writes list of data to file as ascii
  \param v Data
  \param fname Name of file
  \param precision Number of digits
 */
template <class T> void write_to_file(vector<T> const& v, string const& fname,
				      int precision = 6)
{
  std::ofstream f(fname.c_str());
  f.precision(precision);
  for(size_t i=0;i<v.size();++i)
    f << v[i] << "\n";
  f.close();
}

/*! \brief Writes a certain cell property to a file
  \param sim Hydrodynamic simulation
  \param property Name of propery
  \param fname Name of output file
  \param precision Number of digits
 */
void write_cells_property(SRHDSimulation const& sim,
			  string const& property,
			  string const& fname,
			  int precision = 6);

/*! \brief Writes grid and cells to a file
  \param sim Hydrodynamic simulation
  \param fname Name of file
  \param precision Number of digits
 */
void write_snapshot(SRHDSimulation const& sim,
		    string const& fname,
		    int precision=6);

/*! \brief Writes a single number to a file
  \param num Number
  \param fname Name of file
  \param precision Number of digits to write
 */
void write_number(double num,
		  string const& fname,
		  int precision=14);

//! \brief Retrieves properties from a list of primitive variables
class PrimitivePropertyGetter: public Index2Member<double>
{
public:

  /*! \brief Class constructor
    \param sim Simulation
    \param ppm Pointer to member
  */
  PrimitivePropertyGetter(const SRHDSimulation& sim,
			  double Primitive::* ppm);
  
  size_t getLength(void) const;

  double operator()(size_t i) const;

private:
  //! \brief List of primitive variables
  const vector<Primitive>& cells_;

  //! \brief Pointer to member of primitives
  double Primitive::* ppm_;
};

/*! \brief Writes all vector member to a file
  \param list Data to be written
  \param fname Name of file
  \param precision Precision
 */
void write_vector(const vector<double>& list,
		  const string& fname,
		  int precision=14);

#endif // DIAGNOSTICS_HPP
