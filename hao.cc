// LIC// ====================================================================
// LIC// This file forms part of oomph-lib, the object-oriented,
// LIC// multi-physics finite-element library, available
// LIC// at http://www.oomph-lib.org.
// LIC//
// LIC// Copyright (C) 2006-2023 Matthias Heil and Andrew Hazel
// LIC//
// LIC// This library is free software; you can redistribute it and/or
// LIC// modify it under the terms of the GNU Lesser General Public
// LIC// License as published by the Free Software Foundation; either
// LIC// version 2.1 of the License, or (at your option) any later version.
// LIC//
// LIC// This library is distributed in the hope that it will be useful,
// LIC// but WITHOUT ANY WARRANTY; without even the implied warranty of
// LIC// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// LIC// Lesser General Public License for more details.
// LIC//
// LIC// You should have received a copy of the GNU Lesser General Public
// LIC// License along with this library; if not, write to the Free Software
// LIC// Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
// LIC// 02110-1301  USA.
// LIC//
// LIC// The authors may be contacted at oomph-lib@maths.man.ac.uk.
// LIC//
// LIC//====================================================================
// Driver function for a simple beam problem

// OOMPH-LIB includes
#include "generic.h"
#include "beam.h"
#include "meshes/one_d_lagrangian_mesh.h"

using namespace std;
using namespace oomph;

//===================================================================
/// Function-type-object to perform comparison of complex data types
/// Needed to sort the complex eigenvalues into order based on the
/// size of the real part.
//==================================================================
template<class T>
class ComplexLess
{
public:
  /// Comparison. Are the values identical or not?
  bool operator()(const complex<T>& x, const complex<T>& y) const
  {
    return x.real() < y.real();
  }
};

//========start_of_namespace========================
/// Namespace for physical parameters
//==================================================
namespace Global_Physical_Variables
{
  /// Non-dimensional thickness
  double H = 0.0;

  /// Non-dimensional coefficient (FSI)
  double I = 0.0;

  /// Angle between the two arms of the beam
  double Alpha = 0.0;


  // These are parameters that can be set from the command line:

  /// Aspect ratio: Don't change this on the fly; should only be assigned
  /// once before the mesh is generated
  /// first arm length = |q+0.5|, second arm length = |q-0.5|
  // double Q = 0.4;

  // Define the length of the beam in the GeomObejct
  double Stretch_ratio = 0.4;

  /// Initial value for theta_eq in the Newton solve
  double Initial_value_for_theta_eq = 1.57;

  /// Default value for desired ds
  double Ds_default = 1.0e-4;

  // To prevent large solution jumps in critical intervals for I, try to reduce
  // ds specifically in those areas for smoother and more precise results

  /// End point for the first I interval [Interval1_start,Interval1_end]
  double Interval1_start = 0.0;
  double Interval1_end = 0.01;

  /// Start and end points for the second I interval
  /// [Interval2_start,Interval2_end]
  double Interval2_start = 0.07;
  double Interval2_end = 0.08;

  // If the interval_ds is smaller than the default, it will automatically
  // revert to using the default interval

  /// Value of ds for first interval
  double Ds_interval1 = 10.0;

  /// Value of ds for second interval
  double Ds_interval2 = 10.0;

  /// Timescale ratio
  double Lambda_sq = 0.0;

  // Fulid traction from the slender body theory
  // Input: R is the final deformed position after rotation and translation,
  // dR_dt is the first time derivative of R, N is the normal to the deformed
  // configuration, X is the movement of x direction, Y is the movement of y
  // direction, Theta is the orientaion, dX_dt is the first time derivative of
  // X, dY_dt is the first time derivative of Y, dTheta_dt is the first time
  // derivative of Theta,
  // Output: fluid traction from slender body theory
  void fluid_traction(const Vector<double>& R,
                      const Vector<double>& dR_dt,
                      const Vector<double>& N,
                      Vector<double>& traction)
  {
    // Particle velocity from kinematics
    Vector<double> U(2);
    U[0] = dR_dt[0];
    U[1] = dR_dt[1];

    // Compute the traction into the element
    traction[0] =
      (1.0 - 0.5 * N[1] * N[1]) * (R[1] - U[0]) - 0.5 * (N[1] * N[0] * U[1]);

    traction[1] =
      0.5 * N[1] * N[0] * (R[1] - U[0]) - (1.0 - 0.5 * N[0] * N[0]) * U[1];
  }

  // Compute the first derivative of fluid traction with respect to unknowns
  // (element level)
  // traction_f is the fulid traction computed from the slender body theory
  // Input: R_0 is the deformed position before rotation and translation, N is
  // the normal to the deformed configuration, Theta_add_theta_initial is the
  // orientaion plus theta_initial (first arm: theta_initial=0.0, second arm:
  // theta_initial=alpha)
  // Output:
  // dtraction_f_ddotX,dtraction_f_ddotY,dtraction_f_ddotTheta,dtraction_f_ddotR_0_1,dtraction_f_ddotR_0_2
  void compute_dtraction_f_ddot_unknowns(
    const Vector<double>& R_0, // hierher check
    const Vector<double>& N, // hierher check
    const double& Theta_plus_theta_initial,
    Vector<double>& dtraction_f_ddotX,
    Vector<double>& dtraction_f_ddotY,
    Vector<double>& dtraction_f_ddotTheta,
    Vector<double>& dtraction_f_ddotR_0_1,
    Vector<double>& dtraction_f_ddotR_0_2)
  {
    // Assign values by using the expression from Maple
    dtraction_f_ddotX[0] = N[1] * N[1] / 0.2e1 - 0.1e1;
    dtraction_f_ddotX[1] = -N[1] * N[0] / 0.2e1;

    dtraction_f_ddotY[0] = -N[1] * N[0] / 0.2e1;
    dtraction_f_ddotY[1] = N[0] * N[0] / 0.2e1 - 0.1e1;

    dtraction_f_ddotTheta[0] =
      (-R_0[1] * N[1] * N[1] - N[0] * N[1] * R_0[0] + 0.2e1 * R_0[1]) *
        cos(Theta_plus_theta_initial) / 0.2e1 +
      (R_0[1] * N[0] * N[1] - N[1] * N[1] * R_0[0] + 0.2e1 * R_0[0]) *
        sin(Theta_plus_theta_initial) / 0.2e1;

    dtraction_f_ddotTheta[1] =
      (R_0[1] * N[0] * N[1] + N[0] * N[0] * R_0[0] - 0.2e1 * R_0[0]) *
        cos(Theta_plus_theta_initial) / 0.2e1 +
      (-R_0[1] * N[0] * N[0] + N[0] * N[1] * R_0[0] + 0.2e1 * R_0[1]) *
        sin(Theta_plus_theta_initial) / 0.2e1;


    dtraction_f_ddotR_0_1[0] =
      (N[1] * N[1] - 0.2e1) * cos(Theta_plus_theta_initial) / 0.2e1 -
      N[1] * N[0] * sin(Theta_plus_theta_initial) / 0.2e1;

    dtraction_f_ddotR_0_1[1] =
      -N[1] * N[0] * cos(Theta_plus_theta_initial) / 0.2e1 +
      (N[0] * N[0] - 0.2e1) * sin(Theta_plus_theta_initial) / 0.2e1;


    dtraction_f_ddotR_0_2[0] =
      -N[1] * N[0] * cos(Theta_plus_theta_initial) / 0.2e1 +
      (-N[1] * N[1] + 0.2e1) * sin(Theta_plus_theta_initial) / 0.2e1;

    dtraction_f_ddotR_0_2[1] =
      (N[0] * N[0] - 0.2e1) * cos(Theta_plus_theta_initial) / 0.2e1 +
      N[1] * N[0] * sin(Theta_plus_theta_initial) / 0.2e1;
  }


} // namespace Global_Physical_Variables


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=========================================================================
/// RigidBodyElement
//=========================================================================
class RigidBodyElement : public GeneralisedElement
{
public:
  /// Constructor: Pass initial values for rigid body parameters (pinned
  /// by default)
  RigidBodyElement(const double& V,
                   const double& U0,
                   const double& Theta_eq,
                   const double& X0,
                   const double& Y0)
  {
    // Create internal data which contains the "rigid body" parameters
    for (unsigned i = 0; i < 5; i++)
    {
      // Create data: One value, no timedependence, free by default
      add_internal_data(new Data(1));
    }

    // Give them a value:
    internal_data_pt(0)->set_value(0, V);
    internal_data_pt(1)->set_value(0, U0);
    internal_data_pt(2)->set_value(0, Theta_eq);

    // These are just initial values so pin
    internal_data_pt(3)->set_value(0, X0);
    internal_data_pt(3)->pin(0);
    internal_data_pt(4)->set_value(0, Y0);
    internal_data_pt(4)->pin(0);
  }


  /// Function that returns the Vector of pointers to the "rigid body"
  /// parameters
  Vector<Data*> rigid_body_parameters()
  {
    Vector<Data*> tmp_pt(5);
    for (unsigned i = 0; i < 5; i++)
    {
      tmp_pt[i] = internal_data_pt(i);
    }
    return tmp_pt;
  }


  /// Helper function to compute the meaningful parameter values
  /// from enumerated data
  void get_parameters(
    double& V, double& U0, double& Theta_eq, double& X0, double& Y0)
  {
    V = internal_data_pt(0)->value(0);
    U0 = internal_data_pt(1)->value(0);
    Theta_eq = internal_data_pt(2)->value(0);
    X0 = internal_data_pt(3)->value(0);
    Y0 = internal_data_pt(4)->value(0);
  }


  /// Pass pointer to the Mesh of HaoHermiteBeamElements
  /// and add their unknowns to be external data for this element
  void set_pointer_to_beam_meshes(const Vector<SolidMesh*>& beam_mesh_pt)
  {
    // Store the pointer for future reference
    Beam_mesh_pt = beam_mesh_pt;

    // Loop over the nodes in the all mesh and add them as external Data
    // because they affect the traction and therefore the total drag
    // and torque on the object
    unsigned npointer = beam_mesh_pt.size();
    for (unsigned i = 0; i < npointer; i++)
    {
      unsigned nnode = beam_mesh_pt[i]->nnode();
      for (unsigned j = 0; j < nnode; j++)
      {
        add_external_data(beam_mesh_pt[i]->node_pt(j)->variable_position_pt());
      }
    }
  }


  /// Compute the beam's centre of mass
  void compute_centre_of_mass(Vector<double>& sum_r_centre);


  /// Compute the drag and torque on the entire beam structure according
  /// to slender body theory
  void compute_drag_and_torque(Vector<double>& sum_total_drag,
                               double& sum_total_torque);


  /// Output the Theta_eq, Theta_eq_orientation (make comparision with paper's
  /// results), drag and torque on the entire beam structure
  void output(std::ostream& outfile)
  {
    Vector<double> sum_total_drag(2);
    double sum_total_torque = 0.0;

    // Compute the drag and torque on the entire beam structure
    compute_drag_and_torque(sum_total_drag, sum_total_torque);

    // Output Theta_eq
    double Theta_eq = internal_data_pt(2)->value(0);
    outfile << fmod(Theta_eq, 2 * acos(-1.0)) << "  ";

    // Make a transformation from Theta_eq to Theta_eq_orientation
    // Note that here Theta_eq_orientation is controlled in the range of
    // [-2*PI,2*PI]
    double Theta_eq_orientation = fmod(
      fmod(Theta_eq, 2.0 * acos(-1.0)) + acos(-1.0) / 2.0, 2.0 * acos(-1.0));

    // To escape the jump of the solutions
    if (fabs(Theta_eq_orientation) > 1.5 * acos(-1.0))
    {
      if (Theta_eq_orientation > 0)
      {
        outfile << Theta_eq_orientation - 2 * acos(-1.0) << "  ";
      }
      else
      {
        outfile << Theta_eq_orientation + 2 * acos(-1.0) << "  ";
      }
    }
    else
    {
      outfile << Theta_eq_orientation << "  ";
    }

    // Output drag and torque on the entire beam structure
    outfile << sum_total_drag[0] << "  ";
    outfile << sum_total_drag[1] << "  ";
    outfile << sum_total_torque << "  ";
  }


protected:
  // Fill in contribution to residuals
  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    // oomph_info << "ndof in element: " << residuals.size() << std::endl;

    // Get current total drag and torque
    Vector<double> sum_total_drag(2);
    double sum_total_torque = 0.0;
    compute_drag_and_torque(sum_total_drag, sum_total_torque);

    unsigned n_internal = ninternal_data();
    for (unsigned i = 0; i < n_internal; i++)
    {
      // Get the local equation number of the zeroth dof
      // associated with this internal Data object
      unsigned j = 0;
      int eqn_number = internal_local_eqn(i, j);

      // Is it an actual dof
      if (eqn_number >= 0)
      {
        if (i == 0)
        {
          // Eqn for V:
          residuals[eqn_number] = sum_total_drag[0];
          // internal_data_pt(i)->value(j)-
          // Global_Physical_Variables::V;
        }
        else if (i == 1)
        {
          // Eqn for U0:
          residuals[eqn_number] = sum_total_drag[1];
          // internal_data_pt(i)->value(j)-
          // Global_Physical_Variables::U0;
        }
        else if (i == 2)
        {
          // Eqn for Theta_eq:
          residuals[eqn_number] = sum_total_torque;
          // internal_data_pt(i)->value(j)-
          // Global_Physical_Variables::Theta_eq;
        }
        else
        {
          oomph_info << "Never get here\n";
          abort();
        }

        // std::cout << "internal data " << i << " is not pinned\n";
      }
      else
      {
        // std::cout << "internal data " << i << " is pinned\n";
      }
    }
  }

private:
  /// Pointer to the Mesh of HaoHermiteBeamElements
  Vector<SolidMesh*> Beam_mesh_pt;
};


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


//=====================================================================
/// Upgraded Hermite Beam Element to incorporate slender body traction
//=====================================================================
class HaoHermiteBeamElement : public virtual HermiteBeamElement
{
public:
  /// Constructor: Initialise private member data
  HaoHermiteBeamElement()
    : Rigid_body_element_pt(0), I_pt(0), Theta_initial_pt(0)
  {
  }


  /// Pass pointer to RigidBodyElement that contains the rigid body parameters
  void set_pointer_to_rigid_body_element(
    RigidBodyElement* rigid_body_element_pt)
  {
    // Store the pointer for future reference
    Rigid_body_element_pt = rigid_body_element_pt;

    // Get the rigid body parameters
    Vector<Data*> rigid_body_data_pt =
      Rigid_body_element_pt->rigid_body_parameters();

#ifdef PARANOID
    if (rigid_body_data_pt.size() != 5)
    {
      std::ostringstream error_message;
      error_message << "rigid_body_data_pt should have size 5, not "
                    << rigid_body_data_pt.size() << std::endl;

      // loop over all entries
      for (unsigned i = 0; i < 5; i++)
      {
        if (rigid_body_data_pt[i]->nvalue() != 1)
        {
          error_message << "rigid_body_data_pt[" << i
                        << "] should have 1 value, not "
                        << rigid_body_data_pt[i]->nvalue() << std::endl;
        }
      }

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Add the rigid body parameters as the external data for this element
    for (unsigned i = 0; i < 5; i++)
    {
      add_external_data(rigid_body_data_pt[i]);
    }
  }


  /// Pointer to non-dimensional coefficient (FSI)
  double*& i_pt()
  {
    return I_pt;
  }


  /// Pointer to initial angle
  void theta_initial_pt(const double* theta_initial_pt)
  {
    Theta_initial_pt = theta_initial_pt;
  }

  /// Initial angle
  double theta_initial() const
  {
    if (Theta_initial_pt == 0)
    {
      return 0.0;
    }
    else
    {
      return *Theta_initial_pt;
    }
  }


  /// Compute the element's contribution to the (\int r ds) and length of beam
  void compute_contribution_to_int_r_and_length(Vector<double>& int_r,
                                                double& length)
  {
#ifdef PARANOID
    if (int_r.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "int_r should have size 2, not " << int_r.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise
    int_r[0] = 0.0;
    int_r[1] = 0.0;
    length = 0.0;

    // Local coordinate (1D!)
    Vector<double> s(1);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Return local coordinate s[j]  of i-th integration point.
      unsigned j = 0;
      s[j] = integral_pt()->knot(ipt, j);

      // Get position vector to and non-unit tangent vector on wall:
      // dr/ds. NOTE: This is before we apply the rigid body motion!
      // so in terms of the write-up the position vector is R_0
      Vector<double> R_0(2);
      Vector<double> drds(2);
      get_non_unit_tangent(s, R_0, drds);

      // Jacobian of mapping between local and global coordinates
      double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Translate rigid body parameters into meaningful variables
      // (Type=0: first arm, Type=1: second arm.)
      double V = 0.0;
      double U0 = 0.0;
      double Theta_eq = 0.0;
      double X0 = 0.0;
      double Y0 = 0.0;
      Rigid_body_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

      // Note that we're looking for an pseudo "equilibrium position"
      // where the angle (and the traction!) remain constant while
      // the beam still moves as a rigid body!
      double t = 0.0;


      // hierher use Theta_initial everywhere whenever you're processing
      // Theta_eq

      // Apply rigid body translation and rotation to get the actual
      // shape of the deformed body in the fluid
      Vector<double> R(2);
      R[0] = cos(Theta_eq + theta_initial()) * R_0[0] -
             sin(Theta_eq + theta_initial()) * R_0[1] + 0.5 * V * t * t +
             U0 * t + X0;
      R[1] = sin(Theta_eq + theta_initial()) * R_0[0] +
             cos(Theta_eq + theta_initial()) * R_0[1] + V * t + Y0;

      // Add 'em.
      length += W;
      int_r[0] += R[0] * W;
      int_r[1] += R[1] * W;
    }
  }


  /// Compute the slender body traction acting on the actual beam onto the
  /// element at local coordinate s
  void compute_slender_body_traction_on_actual_beam(const Vector<double>& s,
                                                    Vector<double>& traction)
  {
#ifdef PARANOID
    if (traction.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "traction should have size 2, not " << traction.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get the Eulerian position and the unit normal.
    // NOTE: This is before we apply the rigid body motion!
    // so in terms of the write-up the position vector is R_0 and N_0
    Vector<double> R_0(2);
    Vector<double> N_0(2);
    get_normal(s, R_0, N_0);

    // Translate rigid body parameters into meaningful variables
    // (Type=0: first arm, Type=1: second arm.)
    double V = 0.0;
    double U0 = 0.0;
    double Theta_eq = 0.0;
    double X0 = 0.0;
    double Y0 = 0.0;
    Rigid_body_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

    // Note that we're looking for an pseudo "equilibrium position"
    // where the angle (and the traction!) remain constant while
    // the beam still moves as a rigid body!
    double t = 0.0;

    // Compute R which is after translation and rotation
    Vector<double> R(2);
    R[0] = cos(Theta_eq + theta_initial()) * R_0[0] -
           sin(Theta_eq + theta_initial()) * R_0[1] + 0.5 * V * t * t + U0 * t +
           X0;
    R[1] = sin(Theta_eq + theta_initial()) * R_0[0] +
           cos(Theta_eq + theta_initial()) * R_0[1] + V * t + Y0;

    // Compute normal N which is after translation and rotation
    Vector<double> N(2);
    N[0] = cos(Theta_eq + theta_initial()) * N_0[0] -
           sin(Theta_eq + theta_initial()) * N_0[1];
    N[1] = sin(Theta_eq + theta_initial()) * N_0[0] +
           cos(Theta_eq + theta_initial()) * N_0[1];

    // Compute the traction onto the element at local coordinate s
    traction[0] = 0.5 * (V * t - R[1] + U0) * N[1] * N[1] -
                  0.5 * N[1] * N[0] * V - V * t - U0 + R[1];

    traction[1] =
      0.5 * V * N[0] * N[0] - 0.5 * N[1] * (V * t - R[1] + U0) * N[0] - V;
  }


  /// Compute the slender body traction acting on the beam in the reference
  /// configuration (i.e. without rigid body motion!) at local coordinate s
  void compute_slender_body_traction_on_beam_in_reference_configuration(
    const Vector<double>& s, Vector<double>& traction_0)
  {
#ifdef PARANOID
    if (traction_0.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "traction_0 should have size 2, not "
                    << traction_0.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Translate rigid body parameters into meaningful variables
    double V = 0.0;
    double U0 = 0.0;
    double Theta_eq = 0.0;
    double X0 = 0.0;
    double Y0 = 0.0;
    Rigid_body_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

    // Compute the slender body traction acting on the actual beam onto the
    // element at local coordinate s
    Vector<double> traction(2);
    compute_slender_body_traction_on_actual_beam(s, traction);

    // Rotate the traction from the actual beam back to the reference
    // configuration.
    traction_0[0] = traction[0] * cos(Theta_eq + theta_initial()) +
                    traction[1] * sin(Theta_eq + theta_initial());
    traction_0[1] = -traction[0] * sin(Theta_eq + theta_initial()) +
                    traction[1] * cos(Theta_eq + theta_initial());
  }


  // overloaded load_vector to apply the computed traction_0 (i.e. the
  // traction acting on the beam before its rigid body motion is applied)
  // including the non-dimensional coefficient I (FSI)
  void load_vector(const unsigned& intpt,
                   const Vector<double>& xi,
                   const Vector<double>& x,
                   const Vector<double>& N,
                   Vector<double>& load)
  {
    /// Return local coordinate s[j] at the specified integration point.
    Vector<double> s(1);
    unsigned j = 0;
    s[j] = integral_pt()->knot(intpt, j);

    compute_slender_body_traction_on_beam_in_reference_configuration(s, load);
    load[0] = *(i_pt()) * load[0];
    load[1] = *(i_pt()) * load[1];
  }


  // Compute the element's contribution to the total drag and torque on
  // the entire beam structure according to slender body theory
  void compute_contribution_to_drag_and_torque(Vector<double>& drag,
                                               double& torque)
  {
#ifdef PARANOID
    if (drag.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "drag should have size 2, not " << drag.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise
    drag[0] = 0.0;
    drag[1] = 0.0;
    torque = 0.0;

    // Compute the beam's positon of centre of mass
    Vector<double> sum_r_centre(2);
    Rigid_body_element_pt->compute_centre_of_mass(sum_r_centre);

    // Local coordinate (1D!)
    Vector<double> s(1);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Translate rigid body parameters into meaningful variables
    // (Type=0: first arm, Type=1: second arm.)
    double V = 0.0;
    double U0 = 0.0;
    double Theta_eq = 0.0;
    double X0 = 0.0;
    double Y0 = 0.0;
    Rigid_body_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

    // Note that we're looking for an pseudo "equilibrium position"
    // where the angle (and the traction!) remain constant while
    // the beam still moves as a rigid body!
    double t = 0.0;

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      /// Return local coordinate s[j] of i-th integration point.
      unsigned j = 0;
      s[j] = integral_pt()->knot(ipt, j);


      // Get position vector to and non-unit tangent vector on wall:
      // dr/ds
      // NOTE: This is before we apply the rigid body motion!
      // so in terms of the write-up the position vector is R_0
      Vector<double> R_0(2);
      Vector<double> drds(2);
      get_non_unit_tangent(s, R_0, drds);

      // Jacobian. Since Jacobian is the same for R, still use it here.
      double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Compute the slender body traction on actual beam; note this is
      // inefficient since we've already computed certain quantities that will
      // be needed in this function
      Vector<double> traction(2);
      compute_slender_body_traction_on_actual_beam(s, traction);

      // Compute R (after translation and rotation)
      Vector<double> R(2);
      R[0] = cos(Theta_eq + theta_initial()) * R_0[0] -
             sin(Theta_eq + theta_initial()) * R_0[1] + 0.5 * V * t * t +
             U0 * t + X0;
      R[1] = sin(Theta_eq + theta_initial()) * R_0[0] +
             cos(Theta_eq + theta_initial()) * R_0[1] + V * t + Y0;

      // calculate the contribution to torque
      double local_torque =
        (R[0] - (0.5 * V * t * t + U0 * t + X0)) * traction[1] -
        (R[1] - (V * t + Y0)) * traction[0];

      // Add 'em
      drag[0] += traction[0] * W;
      drag[1] += traction[1] * W;
      torque += local_torque * W;
    }
  }


  /// Overloaded output function
  void output(std::ostream& outfile, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_dofs = nnodal_position_type();

    Vector<double> R_0(n_dim);

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_dofs);

    // Loop over element plot points
    for (unsigned l1 = 0; l1 < n_plot; l1++)
    {
      s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

      // Get shape functions
      shape(s, psi);

      Vector<double> interpolated_xi(n_lagrangian);
      interpolated_xi[0] = 0.0;

      // Loop over coordinate directions/components of Vector
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Initialise
        R_0[i] = 0.0;
      }

      // Calculate positions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_dofs; k++)
        {
          // Loop over Lagrangian coordinate directions [xi_gen[] are the
          // the *gen*eralised Lagrangian coordinates: node, type, direction]
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over components of the deformed position Vector
          for (unsigned i = 0; i < n_dim; i++)
          {
            R_0[i] += raw_dnodal_position_gen_dt(0, l, k, i) * psi(l, k);
          }
        }
      }

      // Get the normal vector N0 at each plotted point
      Vector<double> N_0(n_dim);
      get_normal(s, N_0);

      // Compute slender body traction acting on the actual beam
      Vector<double> traction(n_dim);
      compute_slender_body_traction_on_actual_beam(s, traction);

      // Compute slender body traction acting on the beam in the reference
      // configuration
      Vector<double> traction_0(n_dim);
      compute_slender_body_traction_on_beam_in_reference_configuration(
        s, traction_0);

      // Translate rigid body parameters into meaningful variables
      double V = 0.0;
      double U0 = 0.0;
      double Theta_eq = 0.0;
      double X0 = 0.0;
      double Y0 = 0.0;
      Rigid_body_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

      // Note that we're looking for an pseudo "equilibrium position"
      // where the angle (and the traction!) remain constant while
      // the beam still moves as a rigid body!
      double t = 0.0;

      // Compute R after translation and rotation
      Vector<double> R(n_dim);
      R[0] = cos(Theta_eq + theta_initial()) * R_0[0] -
             sin(Theta_eq + theta_initial()) * R_0[1] + 0.5 * V * t * t +
             U0 * t + X0;
      R[1] = sin(Theta_eq + theta_initial()) * R_0[0] +
             cos(Theta_eq + theta_initial()) * R_0[1] + V * t + Y0;

      // Compute normal N after translation and rotation
      Vector<double> N(n_dim);
      N[0] = cos(Theta_eq + theta_initial()) * N_0[0] -
             sin(Theta_eq + theta_initial()) * N_0[1];
      N[1] = sin(Theta_eq + theta_initial()) * N_0[0] +
             cos(Theta_eq + theta_initial()) * N_0[1];

      // Output R0 which is clamped at the origin
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << R_0[i] << " ";
      }

      // Output R which is after translation and rotation
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << R[i] << " ";
      }

      // Output unit normal N0
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << N_0[i] << " ";
      }

      // Output unit normal N which is after translation and rotation
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << N[i] << " ";
      }

      // Output traction acting on the beam in the reference configuration
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << traction_0[i] << " ";
      }

      // Output traction acting on the actual beam
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << traction[i] << " ";
      }

      // Output the velocity of the background
      outfile << R[1] << "  " << 0;
      outfile << std::endl;
    }
  }

private:
  /// Pointer to element that controls the rigid body motion
  RigidBodyElement* Rigid_body_element_pt;

  /// Pointer to non-dimensional coefficient (FSI)
  double* I_pt;

  /// Pointer to initial rotation of the element when it's in its (otherwise)
  /// undeformed configuration
  const double* Theta_initial_pt;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////


//=============================================================================
/// Compute the beam's centre of mass (defined outside class to avoid
/// forward references)
//=============================================================================
void RigidBodyElement::compute_centre_of_mass(Vector<double>& sum_r_centre)
{
#ifdef PARANOID
  if (sum_r_centre.size() != 2)
  {
    std::ostringstream error_message;
    error_message << "sum_r_centre should have size 2, not "
                  << sum_r_centre.size() << std::endl;

    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Initialise
  sum_r_centre[0] = 0.0;
  sum_r_centre[1] = 0.0;
  Vector<double> int_r(2);
  double length = 0.0;

  // Find number of beam meshes
  unsigned npointer = Beam_mesh_pt.size();

  // Loop over the beam meshes to compute the centre of mass of the entire beam
  for (unsigned i = 0; i < npointer; i++)
  {
    // Initialise
    Vector<double> total_int_r(2);
    double total_length = 0.0;

    // Find number of elements in the mesh
    unsigned n_element = Beam_mesh_pt[i]->nelement();

    // Loop over the elements to compute the sum of elements' contribution to
    // the (\int r ds) and the length of beam
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast to the specific element type
      HaoHermiteBeamElement* elem_pt =
        dynamic_cast<HaoHermiteBeamElement*>(Beam_mesh_pt[i]->element_pt(e));

      // Compute contribution to the the (\int r ds) and length of beam within
      // the e-th element
      elem_pt->compute_contribution_to_int_r_and_length(int_r, length);

      // Sum the elements' contribution to the (\int r ds) and length of beam
      total_int_r[0] += int_r[0];
      total_int_r[1] += int_r[1];
      total_length += length;
    } // end of loop over elements

    // Assemble the (\int r ds) and beam length to get the centre of mass for
    // one arm
    Vector<double> r_centre(2);
    r_centre[0] = (1.0 / total_length) * total_int_r[0];
    r_centre[1] = (1.0 / total_length) * total_int_r[1];

    // Compute the centre of mass of the entire beam
    sum_r_centre[0] = sum_r_centre[0] + r_centre[0];
    sum_r_centre[1] = sum_r_centre[1] + r_centre[1];
  }

  // Get centre of mass
  sum_r_centre[0] = 0.5 * sum_r_centre[0];
  sum_r_centre[1] = 0.5 * sum_r_centre[1];
}


//=============================================================================
/// Compute the drag and torque on the entire beam structure according to
/// slender body theory (Type=0: first arm, Type=1: second arm.)
//=============================================================================
void RigidBodyElement::compute_drag_and_torque(Vector<double>& sum_total_drag,
                                               double& sum_total_torque)
{
#ifdef PARANOID
  if (sum_total_drag.size() != 2)
  {
    std::ostringstream error_message;
    error_message << "sum_total_drag should have size 2, not "
                  << sum_total_drag.size() << std::endl;

    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Initialise
  sum_total_drag[0] = 0.0;
  sum_total_drag[1] = 0.0;
  sum_total_torque = 0.0;

  Vector<double> drag(2);
  double torque = 0.0;

  // Find number of beam meshes
  unsigned npointer = Beam_mesh_pt.size();

  // Loop over the beam meshes to compute the drag and torque of the entire beam
  for (unsigned i = 0; i < npointer; i++)
  {
    // Initialise
    Vector<double> total_drag(2);
    double total_torque = 0.0;

    // Find number of elements in the mesh
    unsigned n_element = Beam_mesh_pt[i]->nelement();

    // Loop over the elements to compute the sum of elements' contribution to
    // the drag and torque on the entire beam structure
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast to the specific element type
      HaoHermiteBeamElement* elem_pt =
        dynamic_cast<HaoHermiteBeamElement*>(Beam_mesh_pt[i]->element_pt(e));

      // Compute contribution to the drag and torque within the e-th element
      elem_pt->compute_contribution_to_drag_and_torque(drag, torque);

      // Sum the elements' contribution to the drag and torque
      total_drag[0] += drag[0];
      total_drag[1] += drag[1];
      total_torque += torque;
    } // end of loop over elements

    // Compute the drag and torque of the entire beam
    sum_total_drag[0] = sum_total_drag[0] + total_drag[0];
    sum_total_drag[1] = sum_total_drag[1] + total_drag[1];
    sum_total_torque = sum_total_torque + total_torque;
  }
}


//======start_of_problem_class==========================================
/// Beam problem object
//======================================================================
class ElasticBeamProblem : public Problem
{
public:
  /// Constructor: The arguments are the number of elements and the parameter to
  /// determine the length of the beam
  ElasticBeamProblem(const unsigned& n_elem1, const unsigned& n_elem2);

  /// Conduct a parameter study
  // void parameter_study();

  /// No actions need to be performed after a solve
  void actions_after_newton_solve() {}

  /// No actions need to be performed before a solve
  void actions_before_newton_solve() {}

  /// Dump problem data to allow for later restart
  void dump_it(ofstream& dump_file)
  {
    // Current value of the FSI parameter
    dump_file << Global_Physical_Variables::I << " # FSI parameter"
              << std::endl;

    // hierher maybe add q and alpha but then issue warning
    // if it differs from the one specified on the command line

    // Dump the refinement pattern and the generic problem data
    Problem::dump(dump_file);
  }

  /// Read problem data for restart
  void restart(ifstream& restart_file)
  {
    // Read line up to termination sign
    string input_string;
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read in FSI parameter
    Global_Physical_Variables::I = double(atof(input_string.c_str()));

    // Refine the mesh and read in the generic problem data
    Problem::read(restart_file);
  }

  // Get arms pointers
  void get_steady_problem_beam_meshes(
    SolidMesh*& steady_beam_mesh_first_arm_pt,
    SolidMesh*& steady_beam_mesh_second_arm_pt)
  {
    steady_beam_mesh_first_arm_pt = Beam_mesh_first_arm_pt;
    steady_beam_mesh_second_arm_pt = Beam_mesh_second_arm_pt;
  }

  // Get Rigid_body_element pointers
  void get_steady_problem_rigid_body_element_pt(
    RigidBodyElement*& steady_rigid_body_element_pt)
  {
    steady_rigid_body_element_pt = Rigid_body_element_pt;
  }

private:
  /// Pointer to geometric object that represents the beam's undeformed shape
  GeomObject* Undef_beam_pt1;

  GeomObject* Undef_beam_pt2;

  /// Pointer to RigidBodyElement that actually contains the rigid body data
  RigidBodyElement* Rigid_body_element_pt;

  /// Pointer to beam mesh (first arm)
  OneDLagrangianMesh<HaoHermiteBeamElement>* Beam_mesh_first_arm_pt;

  /// Pointer to beam mesh (second arm)
  OneDLagrangianMesh<HaoHermiteBeamElement>* Beam_mesh_second_arm_pt;

  /// Pointer to mesh containing the rigid body element
  Mesh* Rigid_body_element_mesh_pt;

}; // end of problem class


//=========================================================================
/// Steady, straight 1D line in 2D space with stretch_ratio
///  \f[ x = 0.0 \f]
///  \f[ y = \zeta*stretch_ratio  \f]
//=========================================================================
class StraightLineVertical_new : public GeomObject
{
public:
  /// Constructor derives from GeomObject(1, 2)
  /// Constructor: Pass stretch_ratio
  StraightLineVertical_new(const double& stretch_ratio) : GeomObject(1, 2)
  {
    Stretch_ratio = stretch_ratio;
  }

  /// Broken copy constructor
  StraightLineVertical_new(const StraightLineVertical_new& dummy) = delete;

  /// Broken assignment operator
  void operator=(const StraightLineVertical_new&) = delete;

  /// Position Vector at Lagrangian coordinate zeta
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    r[0] = 0.0;
    r[1] = zeta[0] * Stretch_ratio;
  }


  /// Derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// Evaluated at current time.
  virtual void dposition(const Vector<double>& zeta,
                         DenseMatrix<double>& drdzeta) const
  {
    // Tangent vector
    drdzeta(0, 0) = 0.0;
    drdzeta(0, 1) = Stretch_ratio;
  }


  /// 2nd derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i). Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }


  /// Posn Vector and its  1st & 2nd derivatives
  /// w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i).
  /// Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          Vector<double>& r,
                          DenseMatrix<double>& drdzeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Position Vector
    r[0] = 0.0;
    r[1] = zeta[0] * Stretch_ratio;

    // Tangent vector
    drdzeta(0, 0) = 0.0;
    drdzeta(0, 1) = Stretch_ratio;

    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }

private:
  /// Define the length of the beam
  double Stretch_ratio;
};


//=============start_of_constructor=====================================
/// Constructor for elastic beam problem
//======================================================================
ElasticBeamProblem::ElasticBeamProblem(const unsigned& n_elem1,
                                       const unsigned& n_elem2)
{
  // Drift speed and acceleration of horizontal motion
  double v = 0.0;

  // Speed of horizontal motion
  double u0 = 0.0;

  // Beam's inclination
  double theta_eq = Global_Physical_Variables::Initial_value_for_theta_eq;

  // x position of clamped point
  double x0 = 0.0;

  // y position of clamped point
  double y0 = 0.0;

  // Make the RigidBodyElement that stores the parameters for the rigid body
  // motion
  Rigid_body_element_pt = new RigidBodyElement(v, u0, theta_eq, x0, y0);

  // Add the rigid body element to its own mesh
  Rigid_body_element_mesh_pt = new Mesh;
  Rigid_body_element_mesh_pt->add_element_pt(Rigid_body_element_pt);

  // Set the undeformed beam shape for two arms (in the reference orientation
  // before applying the rigid body motion)
  // Undef_beam_pt = new StraightLineVertical();

  // Aspect ratio to determine the length of the beam
  // first arm length = |q+0.5|, second arm length = |q-0.5|
  // double* q_pt = &Global_Physical_Variables::Q;

  // Still use the same expression of q as before to represent the length of the
  // two arms
  double* stretch_ratio_pt = &Global_Physical_Variables::Stretch_ratio;
  Undef_beam_pt1 = new StraightLineVertical_new(fabs(*stretch_ratio_pt + 0.5));
  Undef_beam_pt2 = new StraightLineVertical_new(fabs(*stretch_ratio_pt - 0.5));

  // Create the (Lagrangian!) mesh, using the StraightLineVertical object
  // Undef_beam_pt to specify the initial (Eulerian) position of the
  // nodes. (first arm)
  // double length_1 = fabs(*q_pt + 0.5);
  double length_1 = 1.0;
  Beam_mesh_first_arm_pt = new OneDLagrangianMesh<HaoHermiteBeamElement>(
    n_elem1, length_1, Undef_beam_pt1);

  // Create the (Lagrangian!) mesh, using the StraightLineVertical object
  // Undef_beam_pt to specify the initial (Eulerian) position of the
  // nodes. (second arm)
  // double length_2 = fabs(*q_pt - 0.5);
  double length_2 = 1.0;
  Beam_mesh_second_arm_pt = new OneDLagrangianMesh<HaoHermiteBeamElement>(
    n_elem2, length_2, Undef_beam_pt2);

  // Pass the pointer of the mesh to the RigidBodyElement class
  // so it can work out the drag and torque on the entire structure
  Vector<SolidMesh*> Beam_mesh_pt(2);
  Beam_mesh_pt[0] = Beam_mesh_first_arm_pt;
  Beam_mesh_pt[1] = Beam_mesh_second_arm_pt;
  Rigid_body_element_pt->set_pointer_to_beam_meshes(Beam_mesh_pt);

  // Build the problem's global mesh
  add_sub_mesh(Beam_mesh_first_arm_pt);
  add_sub_mesh(Beam_mesh_second_arm_pt);
  add_sub_mesh(Rigid_body_element_mesh_pt);
  build_global_mesh();

  // Set the boundary conditions: One end of the beam is clamped in space
  // Pin displacements in both x and y directions, and pin the derivative of
  // position Vector w.r.t. to coordinates in x direction. (first arm)
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(0);
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(1);
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(1, 0);

  // Find number of elements in the mesh (first arm)
  unsigned n_element = Beam_mesh_first_arm_pt->nelement();

  // Loop over the elements to set physical parameters etc. (first arm)
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast to the specific element type
    HaoHermiteBeamElement* elem_pt = dynamic_cast<HaoHermiteBeamElement*>(
      Beam_mesh_first_arm_pt->element_pt(e));

    // Pass the pointer of RigidBodyElement to the each element
    // so we can work out the rigid body motion
    elem_pt->set_pointer_to_rigid_body_element(Rigid_body_element_pt);

    // Set physical parameters for each element:
    elem_pt->h_pt() = &Global_Physical_Variables::H;
    elem_pt->i_pt() = &Global_Physical_Variables::I;
    elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

    // Note: no rotation!

    // Set the undeformed shape for each element
    elem_pt->undeformed_beam_pt() = Undef_beam_pt1;

  } // end of loop over elements


  // Set the boundary conditions: One end of the beam is clamped in space
  // Pin displacements in both x and y directions, and pin the derivative of
  // position Vector w.r.t. to coordinates in x direction. (second arm)
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(0);
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(1);
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(1, 0);

  // Find number of elements in the mesh (second arm)
  n_element = Beam_mesh_second_arm_pt->nelement();

  // Loop over the elements to set physical parameters etc. (second arm)
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast to the specific element type
    HaoHermiteBeamElement* elem_pt = dynamic_cast<HaoHermiteBeamElement*>(
      Beam_mesh_second_arm_pt->element_pt(e));

    // Pass the pointer of RigidBodyElement to the each element
    // so we can work out the rigid body motion
    elem_pt->set_pointer_to_rigid_body_element(Rigid_body_element_pt);

    // Set physical parameters for each element:
    elem_pt->h_pt() = &Global_Physical_Variables::H;
    elem_pt->i_pt() = &Global_Physical_Variables::I;
    elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

    // Rotate by opening angle
    elem_pt->theta_initial_pt(&Global_Physical_Variables::Alpha);

    // Set the undeformed shape for each element
    elem_pt->undeformed_beam_pt() = Undef_beam_pt2;

  } // end of loop over elements

  // Assign the global and local equation numbers
  cout << "# of dofs " << assign_eqn_numbers() << std::endl;

  // Problem::Max_residuals = 1.0e10;
  //  Problem::Max_newton_iterations = 20;
  Problem::Always_take_one_newton_step = true;
  Problem::Scale_arc_length = false;
  Problem::Theta_squared = 0.3;

} // end of constructor


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////


//=========================================================================
/// UnsteadyRigidBodyElement
//=========================================================================
class UnsteadyRigidBodyElement : public GeneralisedElement
{
public:
  /// Constructor: Pass initial values for rigid body parameters (pinned
  /// by default) and time stepper pointer
  UnsteadyRigidBodyElement(const double& X,
                           const double& Y,
                           const double& Theta,
                           TimeStepper*& time_stepper_pt)
    : Beam_mesh_pt(0),
      Time_stepper_pt(0),
      Steady_beam_mesh_pt(0),
      Steady_rigid_body_element_pt(0)
  {
    // Pass time stepper pointer
    Time_stepper_pt = time_stepper_pt;

    // Create internal data which contains the "rigid body" parameters
    for (unsigned i = 0; i < 3; i++)
    {
      // Create data: One value, timedependence, free by default
      add_internal_data(new Data(Time_stepper_pt, 1));
    }

    // Assign these values as internal data with elements
    internal_data_pt(0)->set_value(0, X);
    internal_data_pt(1)->set_value(0, Y);
    internal_data_pt(2)->set_value(0, Theta);
  }


  /// Function that returns the Vector of pointers to the "rigid body"
  /// parameters
  Vector<Data*> rigid_body_parameters()
  {
    Vector<Data*> tmp_pt(3);
    for (unsigned i = 0; i < 3; i++)
    {
      tmp_pt[i] = internal_data_pt(i);
    }
    return tmp_pt;
  }


  /// Helper function to compute the meaningful parameter values
  /// from enumerated data
  void get_parameters(double& X, double& Y, double& Theta)
  {
    X = internal_data_pt(0)->value(0);
    Y = internal_data_pt(1)->value(0);
    Theta = internal_data_pt(2)->value(0);
  }


  /// Helper function to compute the first time derivative of the meaningful
  /// parameter values
  void get_first_time_derivative_of_parameters(double& dX_dt,
                                               double& dY_dt,
                                               double& dTheta_dt)
  {
    dX_dt = Time_stepper_pt->time_derivative(1, internal_data_pt(0), 0);
    dY_dt = Time_stepper_pt->time_derivative(1, internal_data_pt(1), 0);
    dTheta_dt = Time_stepper_pt->time_derivative(1, internal_data_pt(2), 0);
  }


  /// Pass pointer to the Mesh of UnsteadyHaoHermiteBeamElements
  /// and add their unknowns to be external data for this element
  void set_pointer_to_beam_meshes(const Vector<SolidMesh*>& beam_mesh_pt);


  /// Pass steady beam meshes pointer to the Mesh of
  /// UnsteadyHaoHermiteBeamElements
  void set_pointer_to_steady_beam_meshes_and_steady_rigid_body_element(
    const Vector<SolidMesh*>& steady_beam_mesh_pt,
    RigidBodyElement*& steady_rigid_body_element_pt)
  {
    // Store the pointer for future reference
    Steady_beam_mesh_pt = steady_beam_mesh_pt;
    Steady_rigid_body_element_pt = steady_rigid_body_element_pt;
  }


  /// Compute the beam's centre of mass
  void compute_centre_of_mass(Vector<double>& sum_r_centre);


  /// Compute the drag and torque on the entire beam structure according
  /// to slender body theory
  void compute_drag_and_torque(Vector<double>& sum_total_drag,
                               double& sum_total_torque);


  /// Output the Theta_eq, Theta_eq_orientation (make comparision with paper's
  /// results), drag and torque on the entire beam structure
  void output(std::ostream& outfile)
  {
    Vector<double> sum_total_drag(2);
    double sum_total_torque = 0.0;

    // Compute the drag and torque on the entire beam structure
    compute_drag_and_torque(sum_total_drag, sum_total_torque);

    // Output Theta_eq
    double Theta = internal_data_pt(2)->value(0);
    outfile << fmod(Theta, 2.0 * acos(-1.0)) << "  ";

    // Make a transformation from Theta_eq to Theta_eq_orientation
    // Note that here Theta_eq_orientation is controlled in the range of
    // [-2*PI,2*PI]
    double theta_orientation =
      fmod(fmod(Theta, 2.0 * acos(-1.0)) + acos(-1.0) / 2.0, 2.0 * acos(-1.0));

    // To escape the jump of the solutions
    if (fabs(theta_orientation) > 1.5 * acos(-1.0))
    {
      if (theta_orientation > 0)
      {
        outfile << theta_orientation - 2.0 * acos(-1.0) << "  ";
      }
      else
      {
        outfile << theta_orientation + 2.0 * acos(-1.0) << "  ";
      }
    }
    else
    {
      outfile << theta_orientation << "  ";
    }

    // Output X,Y
    double X = internal_data_pt(0)->value(0);
    double Y = internal_data_pt(1)->value(0);
    outfile << X << "  ";
    outfile << Y << "  ";

    // Output drag and torque on the entire beam structure
    outfile << sum_total_drag[0] << "  ";
    outfile << sum_total_drag[1] << "  ";
    outfile << sum_total_torque << "  ";
  }

  // Output the eigenvectors
  void output_eigenvectors(std::ostream& outfile)
  {
    // X
    outfile << internal_data_pt(0)->value(0) << std::endl;
    // Y
    outfile << internal_data_pt(1)->value(0) << std::endl;
    // Theta
    outfile << internal_data_pt(2)->value(0) << std::endl;
  }

  // Find the external data index l and its value index v in terms of mesh g,
  // element e, local node n, position direction i, and position type k
  void external_data_index_and_value_index(const unsigned& g,
                                           const unsigned& e,
                                           const unsigned& n,
                                           const unsigned& i,
                                           const unsigned& k,
                                           unsigned& l,
                                           unsigned& v)
  {
    // Read the external data label from the Mapping
    for (unsigned c = 0; c < Mapping.ncol(); c++)
    {
      if (Mapping(0, c) == g && Mapping(1, c) == e && Mapping(2, c) == n)
      {
        l = Mapping(3, c);
      }
    }

    // Determine the value label v in terms of position direction i, and
    // position type k
    // The order of the values stored at exterrnal data is
    // v=0:x; v=1:dx; v=2:y; v=3:dy
    if (i == 0 && k == 0)
    {
      v = 0;
    }
    else if (i == 0 && k == 1)
    {
      v = 1;
    }
    else if (i == 1 && k == 0)
    {
      v = 2;
    }
    else
    {
      v = 3;
    }
  }

protected:
  // Fill in contribution to residuals
  void fill_in_contribution_to_residuals(Vector<double>& residuals)
  {
    // Get current total drag and torque
    Vector<double> sum_total_drag(2);
    double sum_total_torque = 0.0;
    compute_drag_and_torque(sum_total_drag, sum_total_torque);

    unsigned n_internal = ninternal_data();
    for (unsigned i = 0; i < n_internal; i++)
    {
      // Get the local equation number of the zeroth dof
      // associated with this internal Data object
      unsigned j = 0;
      int eqn_number = internal_local_eqn(i, j);

      // Is it an actual dof
      if (eqn_number >= 0)
      {
        if (i == 0)
        {
          // Eqn for drag_x:
          residuals[eqn_number] = sum_total_drag[0];
        }
        else if (i == 1)
        {
          // Eqn for drag_y:
          residuals[eqn_number] = sum_total_drag[1];
        }
        else if (i == 2)
        {
          // Eqn for torque:
          residuals[eqn_number] = sum_total_torque;
        }
        else
        {
          oomph_info << "Never get here\n";
          abort();
        }

        // std::cout << "internal data " << i << " is not pinned\n";
      }
      else
      {
        // std::cout << "internal data " << i << " is pinned\n";
      }
    }
  }

  /// Assemble the contributions to the jacobian and mass matrices
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix);

private:
  /// Pointer to the Mesh of UnsteadyHaoHermiteBeamElements
  Vector<SolidMesh*> Beam_mesh_pt;

  /// Time stepper pointer
  TimeStepper* Time_stepper_pt;

  /// Pointer to the Mesh of UnsteadyHaoHermiteBeamElements
  Vector<SolidMesh*> Steady_beam_mesh_pt;

  /// Pointer to steady Rigid_body_element_pt
  RigidBodyElement* Steady_rigid_body_element_pt;

  /// Mapping the external data label l with the mesh, element, local node
  DenseMatrix<unsigned> Mapping;
};


/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////


//=====================================================================
/// Upgraded Hermite Beam Element to incorporate slender body traction
//=====================================================================
class UnsteadyHaoHermiteBeamElement : public virtual HermiteBeamElement
{
public:
  /// Constructor: Initialise private member data
  UnsteadyHaoHermiteBeamElement()
    : Rigid_body_element_pt(0),
      I_pt(0),
      Theta_initial_pt(0),
      Steady_rigid_body_element_pt(0)
  {
  }


  /// Pass pointer to UnsteadyRigidBodyElement that contains the rigid body
  /// parameters
  void set_pointer_to_rigid_body_element(
    UnsteadyRigidBodyElement* rigid_body_element_pt)
  {
    // Store the pointer for future reference
    Rigid_body_element_pt = rigid_body_element_pt;

    // Get the rigid body parameters
    Vector<Data*> rigid_body_data_pt =
      Rigid_body_element_pt->rigid_body_parameters();

#ifdef PARANOID
    if (rigid_body_data_pt.size() != 3)
    {
      std::ostringstream error_message;
      error_message << "rigid_body_data_pt should have size 3, not "
                    << rigid_body_data_pt.size() << std::endl;

      // loop over all entries
      for (unsigned i = 0; i < 3; i++)
      {
        if (rigid_body_data_pt[i]->nvalue() != 1)
        {
          error_message << "rigid_body_data_pt[" << i
                        << "] should have 1 value, not "
                        << rigid_body_data_pt[i]->nvalue() << std::endl;
        }
      }

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Add the rigid body parameters as the external data for this element
    for (unsigned i = 0; i < 3; i++)
    {
      add_external_data(rigid_body_data_pt[i]);
    }
  }


  /// Pass pointer to steady RigidBodyElement
  void set_pointer_to_steady_rigid_body_element(
    RigidBodyElement*& steady_rigid_body_element_pt)
  {
    // Store the pointer for future reference
    Steady_rigid_body_element_pt = steady_rigid_body_element_pt;
  }


  /// Pointer to non-dimensional coefficient (FSI)
  double*& i_pt()
  {
    return I_pt;
  }


  /// Pointer to initial angle
  void theta_initial_pt(const double* theta_initial_pt)
  {
    Theta_initial_pt = theta_initial_pt;
  }

  /// Initial angle
  double theta_initial() const
  {
    if (Theta_initial_pt == 0)
    {
      return 0.0;
    }
    else
    {
      return *Theta_initial_pt;
    }
  }


  /// Compute the element's contribution to the (\int r ds) and length of beam
  void compute_contribution_to_int_r_and_length(Vector<double>& int_r,
                                                double& length)
  {
#ifdef PARANOID
    if (int_r.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "int_r should have size 2, not " << int_r.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise
    int_r[0] = 0.0;
    int_r[1] = 0.0;
    length = 0.0;

    // Local coordinate (1D!)
    Vector<double> s(1);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Return local coordinate s[j]  of i-th integration point.
      unsigned j = 0;
      s[j] = integral_pt()->knot(ipt, j);

      // Get position vector to and non-unit tangent vector on wall:
      // dr/ds. NOTE: This is before we apply the rigid body motion!
      // so in terms of the write-up the position vector is R_0
      Vector<double> R_0(2);
      Vector<double> drds(2);
      get_non_unit_tangent(s, R_0, drds);

      // Jacobian of mapping between local and global coordinates
      double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Translate rigid body parameters into meaningful variables
      double X = 0.0;
      double Y = 0.0;
      double Theta = 0.0;
      Rigid_body_element_pt->get_parameters(X, Y, Theta);


      // hierher use Theta_initial everywhere whenever you're processing
      // Theta_eq

      // Apply rigid body translation and rotation to get the actual
      // shape of the deformed body in the fluid
      Vector<double> R(2);
      R[0] = cos(Theta + theta_initial()) * R_0[0] -
             sin(Theta + theta_initial()) * R_0[1] + X;
      R[1] = sin(Theta + theta_initial()) * R_0[0] +
             cos(Theta + theta_initial()) * R_0[1] + Y;

      // Add 'em.
      length += W;
      int_r[0] += R[0] * W;
      int_r[1] += R[1] * W;
    }
  }


  /// Compute the slender body traction acting on the actual beam onto the
  /// element at local coordinate s
  void compute_slender_body_traction_on_actual_beam(const Vector<double>& s,
                                                    Vector<double>& traction)
  {
#ifdef PARANOID
    if (traction.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "traction should have size 2, not " << traction.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Get the Eulerian position and the unit normal.
    // NOTE: This is before we apply the rigid body motion!
    // so in terms of the write-up the position vector is R_0 and N_0
    Vector<double> R_0(2);
    Vector<double> N_0(2);
    get_normal(s, R_0, N_0);

    // Translate rigid body parameters into meaningful variables
    double X = 0.0;
    double Y = 0.0;
    double Theta = 0.0;
    Rigid_body_element_pt->get_parameters(X, Y, Theta);


    // Compute R which is after translation and rotation
    Vector<double> R(2);
    R[0] = cos(Theta + theta_initial()) * R_0[0] -
           sin(Theta + theta_initial()) * R_0[1] + X;
    R[1] = sin(Theta + theta_initial()) * R_0[0] +
           cos(Theta + theta_initial()) * R_0[1] + Y;

    // Compute normal N which is after translation and rotation
    Vector<double> N(2);
    N[0] = cos(Theta + theta_initial()) * N_0[0] -
           sin(Theta + theta_initial()) * N_0[1];
    N[1] = sin(Theta + theta_initial()) * N_0[0] +
           cos(Theta + theta_initial()) * N_0[1];


    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_dofs = nnodal_position_type();

    Vector<double> dR0_dt(n_dim);

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_dofs);

    // Get shape functions
    shape(s, psi);

    Vector<double> interpolated_xi(n_lagrangian);
    interpolated_xi[0] = 0.0;

    // Loop over coordinate directions/components of Vector
    for (unsigned i = 0; i < n_dim; i++)
    {
      // Initialise time derivative of R_0
      dR0_dt[i] = 0.0;
    }


    // Calculate spatial derivatives
    for (unsigned l = 0; l < n_node; l++)
    {
      // Loop over positional dofs
      for (unsigned k = 0; k < n_position_dofs; k++)
      {
        // Loop over Lagrangian coordinate directions [xi_gen[] are the
        // the *gen*eralised Lagrangian coordinates: node, type, direction]
        for (unsigned i = 0; i < n_lagrangian; i++)
        {
          interpolated_xi[i] +=
            raw_lagrangian_position_gen(l, k, i) * psi(l, k);
        }

        // Loop over components of the deformed position Vector
        for (unsigned i = 0; i < n_dim; i++)
        {
          dR0_dt[i] += raw_dnodal_position_gen_dt(1, l, k, i) * psi(l, k);
        }
      }
    }

    // Translate rigid body parameters into meaningful variables (first time
    // derivative)
    double dX_dt = 0.0;
    double dY_dt = 0.0;
    double dTheta_dt = 0.0;
    Rigid_body_element_pt->get_first_time_derivative_of_parameters(
      dX_dt, dY_dt, dTheta_dt);

    // Translate rigid body parameters into meaningful variables
    double V = 0.0;
    double U0 = 0.0;
    double Theta_eq = 0.0;
    double X0 = 0.0;
    double Y0 = 0.0;
    Steady_rigid_body_element_pt->get_parameters(V, U0, Theta_eq, X0, Y0);

    // Note that here is to sove the eigenproblem, thus we set dX_dt=Vt+U0 and
    // dY_dt=V for computing the Jacobian matirx
    double t = 0.0;
    dX_dt = V * t + U0;
    dY_dt = V;

    // The first time derivative of the R (R is the final deformed position
    // after rotation and translation)
    Vector<double> dR_dt(2);
    dR_dt[0] = dX_dt - dTheta_dt * sin(Theta + theta_initial()) * R_0[0] +
               cos(Theta + theta_initial()) * dR0_dt[0] -
               dTheta_dt * cos(Theta + theta_initial()) * R_0[1] -
               sin(Theta + theta_initial()) * dR0_dt[1];
    dR_dt[1] = dY_dt + dTheta_dt * cos(Theta + theta_initial()) * R_0[0] +
               sin(Theta + theta_initial()) * dR0_dt[0] -
               dTheta_dt * sin(Theta + theta_initial()) * R_0[1] +
               cos(Theta + theta_initial()) * dR0_dt[1];

    // Compute the traction onto the element at local coordinate s from slender
    // body theory
    Global_Physical_Variables::fluid_traction(R, dR_dt, N, traction);
  }


  /// Compute the slender body traction acting on the beam in the reference
  /// configuration (i.e. without rigid body motion!) at local coordinate s
  void compute_slender_body_traction_on_beam_in_reference_configuration(
    const Vector<double>& s, Vector<double>& traction_0)
  {
#ifdef PARANOID
    if (traction_0.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "traction_0 should have size 2, not "
                    << traction_0.size() << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Translate rigid body parameters into meaningful variables
    double X = 0.0;
    double Y = 0.0;
    double Theta = 0.0;
    Rigid_body_element_pt->get_parameters(X, Y, Theta);

    // Compute the slender body traction acting on the actual beam onto the
    // element at local coordinate s
    Vector<double> traction(2);
    compute_slender_body_traction_on_actual_beam(s, traction);

    // Rotate the traction from the actual beam back to the reference
    // configuration.
    traction_0[0] = traction[0] * cos(Theta + theta_initial()) +
                    traction[1] * sin(Theta + theta_initial());
    traction_0[1] = -traction[0] * sin(Theta + theta_initial()) +
                    traction[1] * cos(Theta + theta_initial());
  }


  // overloaded load_vector to apply the computed traction_0 (i.e. the
  // traction acting on the beam before its rigid body motion is applied)
  // including the non-dimensional coefficient I (FSI)
  void load_vector(const unsigned& intpt,
                   const Vector<double>& xi,
                   const Vector<double>& x,
                   const Vector<double>& N,
                   Vector<double>& load)
  {
    /// Return local coordinate s[j] at the specified integration point.
    Vector<double> s(1);
    unsigned j = 0;
    s[j] = integral_pt()->knot(intpt, j);

    compute_slender_body_traction_on_beam_in_reference_configuration(s, load);
    load[0] = *(i_pt()) * load[0];
    load[1] = *(i_pt()) * load[1];
  }


  // Compute the element's contribution to the total drag and torque on
  // the entire beam structure according to slender body theory
  void compute_contribution_to_drag_and_torque(Vector<double>& drag,
                                               double& torque)
  {
#ifdef PARANOID
    if (drag.size() != 2)
    {
      std::ostringstream error_message;
      error_message << "drag should have size 2, not " << drag.size()
                    << std::endl;

      throw OomphLibError(
        error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
#endif

    // Initialise
    drag[0] = 0.0;
    drag[1] = 0.0;
    torque = 0.0;

    /*     // Compute the beam's positon of centre of mass
        Vector<double> sum_r_centre(2);
        Rigid_body_element_pt->compute_centre_of_mass(sum_r_centre); */

    // Local coordinate (1D!)
    Vector<double> s(1);

    // Set # of integration points
    const unsigned n_intpt = integral_pt()->nweight();

    // Translate rigid body parameters into meaningful variables
    double X = 0.0;
    double Y = 0.0;
    double Theta = 0.0;
    Rigid_body_element_pt->get_parameters(X, Y, Theta);

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      /// Return local coordinate s[j] of i-th integration point.
      unsigned j = 0;
      s[j] = integral_pt()->knot(ipt, j);


      // Get position vector to and non-unit tangent vector on wall:
      // dr/ds
      // NOTE: This is before we apply the rigid body motion!
      // so in terms of the write-up the position vector is R_0
      Vector<double> R_0(2);
      Vector<double> drds(2);
      get_non_unit_tangent(s, R_0, drds);

      // Jacobian. Since Jacobian is the same for R, still use it here.
      double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Compute the slender body traction on actual beam; note this is
      // inefficient since we've already computed certain quantities that will
      // be needed in this function
      Vector<double> traction(2);
      compute_slender_body_traction_on_actual_beam(s, traction);

      // Compute R (after translation and rotation)
      Vector<double> R(2);
      R[0] = cos(Theta + theta_initial()) * R_0[0] -
             sin(Theta + theta_initial()) * R_0[1] + X;
      R[1] = sin(Theta + theta_initial()) * R_0[0] +
             cos(Theta + theta_initial()) * R_0[1] + Y;

      // calculate the contribution to torque
      double local_torque = (R[0] - X) * traction[1] - (R[1] - Y) * traction[0];

      // Add 'em
      drag[0] += traction[0] * W;
      drag[1] += traction[1] * W;
      torque += local_torque * W;
    }
  }


  /// Overloaded output function
  void output(std::ostream& outfile, const unsigned& n_plot)
  {
    // Local variables
    Vector<double> s(1);

    // Tecplot header info
    outfile << "ZONE I=" << n_plot << std::endl;

    // Set the number of lagrangian coordinates
    unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // Set the dimension of the global coordinates
    unsigned n_dim = Undeformed_beam_pt->ndim();

    // Find out how many nodes there are
    unsigned n_node = nnode();

    // Find out how many positional dofs there are
    unsigned n_position_dofs = nnodal_position_type();

    Vector<double> R_0(n_dim);

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_dofs);

    // Loop over element plot points
    for (unsigned l1 = 0; l1 < n_plot; l1++)
    {
      s[0] = -1.0 + l1 * 2.0 / (n_plot - 1);

      // Get shape functions
      shape(s, psi);

      Vector<double> interpolated_xi(n_lagrangian);
      interpolated_xi[0] = 0.0;

      // Loop over coordinate directions/components of Vector
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Initialise
        R_0[i] = 0.0;
      }

      // Calculate positions
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_dofs; k++)
        {
          // Loop over Lagrangian coordinate directions [xi_gen[] are the
          // the *gen*eralised Lagrangian coordinates: node, type, direction]
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over components of the deformed position Vector
          for (unsigned i = 0; i < n_dim; i++)
          {
            R_0[i] += raw_dnodal_position_gen_dt(0, l, k, i) * psi(l, k);
          }
        }
      }

      // Get the normal vector N0 at each plotted point
      Vector<double> N_0(n_dim);
      get_normal(s, N_0);

      // Compute slender body traction acting on the actual beam
      Vector<double> traction(n_dim);
      compute_slender_body_traction_on_actual_beam(s, traction);

      // Compute slender body traction acting on the beam in the reference
      // configuration
      Vector<double> traction_0(n_dim);
      compute_slender_body_traction_on_beam_in_reference_configuration(
        s, traction_0);

      // Translate rigid body parameters into meaningful variables
      double X = 0.0;
      double Y = 0.0;
      double Theta = 0.0;
      Rigid_body_element_pt->get_parameters(X, Y, Theta);

      // Compute R after translation and rotation
      Vector<double> R(n_dim);
      R[0] = cos(Theta + theta_initial()) * R_0[0] -
             sin(Theta + theta_initial()) * R_0[1] + X;
      R[1] = sin(Theta + theta_initial()) * R_0[0] +
             cos(Theta + theta_initial()) * R_0[1] + Y;

      // Compute normal N after translation and rotation
      Vector<double> N(n_dim);
      N[0] = cos(Theta + theta_initial()) * N_0[0] -
             sin(Theta + theta_initial()) * N_0[1];
      N[1] = sin(Theta + theta_initial()) * N_0[0] +
             cos(Theta + theta_initial()) * N_0[1];

      // Output R0 which is clamped at the origin
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << R_0[i] << " ";
      }

      // Output R which is after translation and rotation
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << R[i] << " ";
      }

      // Output unit normal N0
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << N_0[i] << " ";
      }

      // Output unit normal N which is after translation and rotation
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << N[i] << " ";
      }

      // Output traction acting on the beam in the reference configuration
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << traction_0[i] << " ";
      }

      // Output traction acting on the actual beam
      for (unsigned i = 0; i < n_dim; i++)
      {
        outfile << traction[i] << " ";
      }

      // Output the velocity of the background
      outfile << R[1] << "  " << 0;
      outfile << std::endl;
    }
  }

protected:
  /// Assemble the contributions to the jacobian and mass
  /// matrices
  void fill_in_contribution_to_jacobian_and_mass_matrix(
    Vector<double>& residuals,
    DenseMatrix<double>& jacobian,
    DenseMatrix<double>& mass_matrix)
  {
    // Get the jacobian
    fill_in_contribution_to_jacobian(residuals, jacobian);

    // Mass matirx should be evaluted at the base state

    // Find out how many nodes within this element
    const unsigned n_node = nnode();

    // Find out how many positional dofs within this element
    const unsigned n_position_type = nnodal_position_type();

    // Set the dimension of the global coordinates
    const unsigned n_dim = Undeformed_beam_pt->ndim();

    // Set the number of lagrangian coordinates
    const unsigned n_lagrangian = Undeformed_beam_pt->nlagrangian();

    // # of nodes, # of positional dofs
    Shape psi(n_node, n_position_type);

    // # of nodes, # of positional dofs, # of lagrangian coords (for deriv)
    DShape dpsidxi(n_node, n_position_type, n_lagrangian);

    // # of nodes, # of positional dofs, # of derivs)
    DShape d2psidxi(n_node, n_position_type, n_lagrangian);

    // Set the number of integration points
    unsigned n_intpt = integral_pt()->nweight();

    // Integers to store the local equation and unknown numbers
    int local_eqn = 0;
    int local_unknown = 0;

    // Local coordinate (1D!)
    Vector<double> s(1);

    // Thickness h_0/R
    const double HoR_0 = h(); // i.e. refers to reference thickness 'h_0'

    // The number of the external data
    unsigned n_external = nexternal_data();

    // Precomputed steady orientation
    double theta_eq =
      Steady_rigid_body_element_pt->internal_data_pt(2)->value(0);

    std::cout << "theta_eq(solid)=" << theta_eq << std::endl;

    // The derivative of solid traction with respect to fluid fraction
    Vector<double> dtraction_s_dtraction_f1(2);
    Vector<double> dtraction_s_dtraction_f2(2);

    dtraction_s_dtraction_f1[0] = (*(i_pt())) * cos(theta_eq + theta_initial());
    dtraction_s_dtraction_f1[1] =
      -(*(i_pt())) * sin(theta_eq + theta_initial());

    dtraction_s_dtraction_f2[0] = (*(i_pt())) * sin(theta_eq + theta_initial());
    dtraction_s_dtraction_f2[1] = (*(i_pt())) * cos(theta_eq + theta_initial());

    // Loop over the integration points
    for (unsigned ipt = 0; ipt < n_intpt; ipt++)
    {
      // Get the integral weight
      double w = integral_pt()->weight(ipt);

      // Return local coordinate s[j]  of i-th integration point.
      unsigned j = 0;
      s[j] = integral_pt()->knot(ipt, j);

      // Call the derivatives of the shape functions w.r.t. Lagrangian coords
      double J = d2shape_lagrangian_at_knot(ipt, psi, dpsidxi, d2psidxi);

      // Premultiply the weights and the Jacobian
      double W = w * J;

      // Get the Eulerian position and the unit normal.
      // NOTE: This is before we apply the rigid body motion!
      // so in terms of the write-up the position vector is R_0 and N_0
      Vector<double> N_0(2);
      Vector<double> R_0(2);
      get_normal(s, R_0, N_0);

      // Compute normal N which is after translation and rotation
      Vector<double> N(2);
      N[0] = cos(theta_eq + theta_initial()) * N_0[0] -
             sin(theta_eq + theta_initial()) * N_0[1];
      N[1] = sin(theta_eq + theta_initial()) * N_0[0] +
             cos(theta_eq + theta_initial()) * N_0[1];

      // Orientaion plus theta_initial (first arm: theta_initial=0.0, second
      // arm: theta_initial=alpha)
      double Theta_plus_theta_initial = theta_eq + theta_initial();

      // The first derivative of fluid traction with respect to
      // unknowns (element level)
      Vector<double> dtraction_f_ddotX(2);
      Vector<double> dtraction_f_ddotY(2);
      Vector<double> dtraction_f_ddotTheta(2);
      Vector<double> dtraction_f_ddotR_0_1(2);
      Vector<double> dtraction_f_ddotR_0_2(2);

      // Compute the first derivative of fluid traction with respect to
      // unknowns (element level)
      // traction_f is the fulid traction computed from the slender body
      // theory
      // Input: R_0 is the deformed position before rotation and translation,
      // N is the normal to the deformed configuration,
      // Theta_add_theta_initial is the orientaion plus theta_initial (first
      // arm: theta_initial=0.0, second arm: theta_initial=alpha)
      // Output:
      // dtraction_f_ddotX,dtraction_f_ddotY,dtraction_f_ddotTheta,dtraction_f_ddotR_0_1,dtraction_f_ddotR_0_2
      Global_Physical_Variables::compute_dtraction_f_ddot_unknowns(
        R_0,
        N,
        Theta_plus_theta_initial,
        dtraction_f_ddotX,
        dtraction_f_ddotY,
        dtraction_f_ddotTheta,
        dtraction_f_ddotR_0_1,
        dtraction_f_ddotR_0_2);

      // Calculate local values of lagrangian position and
      // the derivative of global position wrt lagrangian coordinates
      Vector<double> interpolated_xi(n_lagrangian, 0.0),
        interpolated_x(n_dim, 0.0);
      DenseMatrix<double> interpolated_A(n_lagrangian, n_dim);

      // Initialise to zero
      // Loop over coordinate directions/components of Vector
      for (unsigned i = 0; i < n_dim; i++)
      {
        // Loop over derivatives/base Vectors (just one here)
        for (unsigned j = 0; j < n_lagrangian; j++)
        {
          interpolated_A(j, i) = 0.0;
        }
      }

      // Calculate displacements, accelerations and spatial derivatives
      for (unsigned l = 0; l < n_node; l++)
      {
        // Loop over positional dofs
        for (unsigned k = 0; k < n_position_type; k++)
        {
          // Loop over Lagrangian coordinate directions [xi_gen[] are the
          // the *gen*eralised Lagrangian coordinates: node, type, direction]
          for (unsigned i = 0; i < n_lagrangian; i++)
          {
            interpolated_xi[i] +=
              raw_lagrangian_position_gen(l, k, i) * psi(l, k);
          }

          // Loop over components of the deformed position Vector
          for (unsigned i = 0; i < n_dim; i++)
          {
            interpolated_x[i] += raw_nodal_position_gen(l, k, i) * psi(l, k);


            // Loop over derivative directions (just one here)
            for (unsigned j = 0; j < n_lagrangian; j++)
            {
              interpolated_A(j, i) +=
                raw_nodal_position_gen(l, k, i) * dpsidxi(l, k, j);
            }
          }
        }
      }

      // Declare and calculate the undeformed and deformed metric tensor
      // and the strain tensor (these are just scalars)
      double Amet = 0.0;

      // Now calculate the dot product
      for (unsigned k = 0; k < n_dim; k++)
      {
        Amet += interpolated_A(0, k) * interpolated_A(0, k);
      }

      // Calculate the contravariant metric tensors
      double Adet = Amet; // double Aup = 1.0/Amet;

      // Define the wall thickness ratio profile
      double h_ratio = 0.0;

      // Get wall thickness ratio profile
      wall_profile(interpolated_xi, interpolated_x, h_ratio);

      // Thickness h/R
      double HoR = HoR_0 * h_ratio;

      // Assemble the contributions to the mass matrix
      // Loop over the number of nodes within this element
      for (unsigned n2 = 0; n2 < n_node; n2++)
      {
        // Loop over the type of position
        for (unsigned k2 = 0; k2 < n_position_type; k2++)
        {
          // Loop over the coordinate directions
          for (unsigned i2 = 0; i2 < n_dim; i2++)
          {
            // Find the equation number
            local_eqn = position_local_eqn(n2, k2, i2);

            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // The derivative of "External forcing" with respect to solid
              // traction
              Vector<double> dE_dtraction_s1(2);
              Vector<double> dE_dtraction_s2(2);

              dE_dtraction_s1[0] =
                -0.10e1 * h_ratio / HoR * psi(n2, k2) * sqrt(Adet);
              dE_dtraction_s1[1] = 0.0;

              dE_dtraction_s2[0] = 0.0;
              dE_dtraction_s2[1] =
                -0.10e1 * h_ratio / HoR * psi(n2, k2) * sqrt(Adet);

              // Loop over the local nodes
              for (unsigned n = 0; n < n_node; n++)
              {
                // Loop over the type of position
                for (unsigned k = 0; k < n_position_type; k++)
                {
                  // Loop over the coordinate directions
                  for (unsigned i = 0; i < n_dim; i++)
                  {
                    local_unknown = position_local_eqn(n, k, i);

                    // If at a non-zero degree of freedom add in the entry
                    if (local_unknown >= 0)
                    {
                      // Note all entries are added a minus sign compared with
                      // Maple's expressions

                      // First component of External_forcing ("traction")
                      if (i2 == 0)
                      {
                        // Coeff of dR_1nk/dt (x direction)
                        if (i == 0)
                        {
                          mass_matrix(local_eqn, local_unknown) += -(
                            (dE_dtraction_s1[0] * dtraction_s_dtraction_f1[0] *
                               dtraction_f_ddotR_0_1[0] * psi(n, k) +
                             dE_dtraction_s1[0] * dtraction_s_dtraction_f2[0] *
                               dtraction_f_ddotR_0_1[1] * psi(n, k) +
                             dE_dtraction_s2[0] * dtraction_s_dtraction_f1[1] *
                               dtraction_f_ddotR_0_1[0] * psi(n, k) +
                             dE_dtraction_s2[0] * dtraction_s_dtraction_f2[1] *
                               dtraction_f_ddotR_0_1[1] * psi(n, k)) *
                            W);
                        }
                        // Coeff of dR_2nk/dt (y direction)
                        else
                        {
                          mass_matrix(local_eqn, local_unknown) += -(
                            (dE_dtraction_s1[0] * dtraction_s_dtraction_f1[0] *
                               dtraction_f_ddotR_0_2[0] * psi(n, k) +
                             dE_dtraction_s1[0] * dtraction_s_dtraction_f2[0] *
                               dtraction_f_ddotR_0_2[1] * psi(n, k) +
                             dE_dtraction_s2[0] * dtraction_s_dtraction_f1[1] *
                               dtraction_f_ddotR_0_2[0] * psi(n, k) +
                             dE_dtraction_s2[0] * dtraction_s_dtraction_f2[1] *
                               dtraction_f_ddotR_0_2[1] * psi(n, k)) *
                            W);
                        }
                      }
                      // Second component of External_forcing ("traction")
                      else
                      {
                        // Coeff of dR_1nk/dt (x direction)
                        if (i == 0)
                        {
                          mass_matrix(local_eqn, local_unknown) += -(
                            (dE_dtraction_s1[1] * dtraction_s_dtraction_f1[0] *
                               dtraction_f_ddotR_0_1[0] * psi(n, k) +
                             dE_dtraction_s1[1] * dtraction_s_dtraction_f2[0] *
                               dtraction_f_ddotR_0_1[1] * psi(n, k) +
                             dE_dtraction_s2[1] * dtraction_s_dtraction_f1[1] *
                               dtraction_f_ddotR_0_1[0] * psi(n, k) +
                             dE_dtraction_s2[1] * dtraction_s_dtraction_f2[1] *
                               dtraction_f_ddotR_0_1[1] * psi(n, k)) *
                            W);
                        }
                        // Coeff of dR_2nk/dt (y direction)
                        else
                        {
                          mass_matrix(local_eqn, local_unknown) += -(
                            (dE_dtraction_s1[1] * dtraction_s_dtraction_f1[0] *
                               dtraction_f_ddotR_0_2[0] * psi(n, k) +
                             dE_dtraction_s1[1] * dtraction_s_dtraction_f2[0] *
                               dtraction_f_ddotR_0_2[1] * psi(n, k) +
                             dE_dtraction_s2[1] * dtraction_s_dtraction_f1[1] *
                               dtraction_f_ddotR_0_2[0] * psi(n, k) +
                             dE_dtraction_s2[1] * dtraction_s_dtraction_f2[1] *
                               dtraction_f_ddotR_0_2[1] * psi(n, k)) *
                            W);
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }


      // The derivatives of beam equation respect to the external data
      // (dS/dX,dS/dY,dS/dtheta)

      // Loop over the local nodes
      for (unsigned n2 = 0; n2 < n_node; n2++)
      {
        // Loop over the type of position
        for (unsigned k2 = 0; k2 < n_position_type; k2++)
        {
          // Loop over the coordinate directions
          for (unsigned i2 = 0; i2 < n_dim; i2++)
          {
            // Find the equation number
            local_eqn = position_local_eqn(n2, k2, i2);


            /*IF it's not a boundary condition*/
            if (local_eqn >= 0)
            {
              // The derivative of "External forcing" with respect to solid
              // traction
              Vector<double> dE_dtraction_s1(2);
              Vector<double> dE_dtraction_s2(2);

              dE_dtraction_s1[0] =
                -0.10e1 * h_ratio / HoR * psi(n2, k2) * sqrt(Adet);
              dE_dtraction_s1[1] = 0.0;

              dE_dtraction_s2[0] = 0.0;
              dE_dtraction_s2[1] =
                -0.10e1 * h_ratio / HoR * psi(n2, k2) * sqrt(Adet);

              // Loop over the external data
              for (unsigned m = 0; m < n_external; m++)
              {
                local_unknown = external_local_eqn(m, 0);

                // If at a non-zero degree of freedom add in the entry
                if (local_unknown >= 0)
                {
                  // First component of External_forcing ("traction")
                  if (i2 == 0)
                  {
                    // Coeff of dX0/dt
                    if (m == 0)
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        -((dE_dtraction_s1[0] * dtraction_s_dtraction_f1[0] *
                             dtraction_f_ddotX[0] +
                           dE_dtraction_s1[0] * dtraction_s_dtraction_f2[0] *
                             dtraction_f_ddotX[1] +
                           dE_dtraction_s2[0] * dtraction_s_dtraction_f1[1] *
                             dtraction_f_ddotX[0] +
                           dE_dtraction_s2[0] * dtraction_s_dtraction_f2[1] *
                             dtraction_f_ddotX[1]) *
                          W);
                    }
                    // Coeff of dY0/dt
                    else if (m == 1)
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        -((dE_dtraction_s1[0] * dtraction_s_dtraction_f1[0] *
                             dtraction_f_ddotY[0] +
                           dE_dtraction_s1[0] * dtraction_s_dtraction_f2[0] *
                             dtraction_f_ddotY[1] +
                           dE_dtraction_s2[0] * dtraction_s_dtraction_f1[1] *
                             dtraction_f_ddotY[0] +
                           dE_dtraction_s2[0] * dtraction_s_dtraction_f2[1] *
                             dtraction_f_ddotY[1]) *
                          W);
                    }
                    // Coeff of dtheta/dt
                    else
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        -((dE_dtraction_s1[0] * dtraction_s_dtraction_f1[0] *
                             dtraction_f_ddotTheta[0] +
                           dE_dtraction_s1[0] * dtraction_s_dtraction_f2[0] *
                             dtraction_f_ddotTheta[1] +
                           dE_dtraction_s2[0] * dtraction_s_dtraction_f1[1] *
                             dtraction_f_ddotTheta[0] +
                           dE_dtraction_s2[0] * dtraction_s_dtraction_f2[1] *
                             dtraction_f_ddotTheta[1]) *
                          W);
                    }
                  }
                  // Second component of External_forcing ("traction")
                  else
                  {
                    // Coeff of dX0/dt
                    if (m == 0)
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        -((dE_dtraction_s1[1] * dtraction_s_dtraction_f1[0] *
                             dtraction_f_ddotX[0] +
                           dE_dtraction_s1[1] * dtraction_s_dtraction_f2[0] *
                             dtraction_f_ddotX[1] +
                           dE_dtraction_s2[1] * dtraction_s_dtraction_f1[1] *
                             dtraction_f_ddotX[0] +
                           dE_dtraction_s2[1] * dtraction_s_dtraction_f2[1] *
                             dtraction_f_ddotX[1]) *
                          W);
                    }
                    // Coeff of dY0/dt
                    else if (m == 1)
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        -((dE_dtraction_s1[1] * dtraction_s_dtraction_f1[0] *
                             dtraction_f_ddotY[0] +
                           dE_dtraction_s1[1] * dtraction_s_dtraction_f2[0] *
                             dtraction_f_ddotY[1] +
                           dE_dtraction_s2[1] * dtraction_s_dtraction_f1[1] *
                             dtraction_f_ddotY[0] +
                           dE_dtraction_s2[1] * dtraction_s_dtraction_f2[1] *
                             dtraction_f_ddotY[1]) *
                          W);
                    }
                    // Coeff of dtheta/dt
                    else
                    {
                      mass_matrix(local_eqn, local_unknown) +=
                        -((dE_dtraction_s1[1] * dtraction_s_dtraction_f1[0] *
                             dtraction_f_ddotTheta[0] +
                           dE_dtraction_s1[1] * dtraction_s_dtraction_f2[0] *
                             dtraction_f_ddotTheta[1] +
                           dE_dtraction_s2[1] * dtraction_s_dtraction_f1[1] *
                             dtraction_f_ddotTheta[0] +
                           dE_dtraction_s2[1] * dtraction_s_dtraction_f2[1] *
                             dtraction_f_ddotTheta[1]) *
                          W);
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    // Output the information of jacobian
    std::cout << "number of rows of the solid element jacobian: "
              << jacobian.nrow() << std::endl;
    std::cout << "number of columns of the solid element jacobian: "
              << jacobian.ncol() << std::endl;

    // Output jacobian
    std::cout << "solid_element_jacobian: " << std::endl;
    for (unsigned i = 0; i < jacobian.nrow(); i++)
    {
      // Output the i-th row entries
      std::cout << "[";
      for (unsigned j = 0; j < jacobian.ncol(); j++)
      {
        std::cout << jacobian(i, j) << ", ";
      }
      std::cout << "]" << std::endl;
    }

    // Output mass matrix
    std::cout << "solid_element_mass_matrix: " << std::endl;
    for (unsigned i = 0; i < mass_matrix.nrow(); i++)
    {
      // Output the i-th row entries
      std::cout << "[";
      for (unsigned j = 0; j < mass_matrix.ncol(); j++)
      {
        std::cout << mass_matrix(i, j) << ", ";
      }
      std::cout << "]" << std::endl;
    }
  }


private:
  /// Pointer to element that controls the rigid body motion
  UnsteadyRigidBodyElement* Rigid_body_element_pt;

  /// Pointer to non-dimensional coefficient (FSI)
  double* I_pt;

  /// Pointer to initial rotation of the element when it's in its
  /// (otherwise) undeformed configuration
  const double* Theta_initial_pt;

  /// Pointer to steady Rigid_body_element_pt
  RigidBodyElement* Steady_rigid_body_element_pt;
};


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////

//=============================================================================
/// Pass pointer to the Mesh of UnsteadyHaoHermiteBeamElements
/// and add their unknowns to be external data for this element
//=============================================================================
void UnsteadyRigidBodyElement::set_pointer_to_beam_meshes(
  const Vector<SolidMesh*>& beam_mesh_pt)
{
  // Store the pointer for future reference
  Beam_mesh_pt = beam_mesh_pt;

  // Find number of elements in the first and second meshes
  unsigned n_element1 = Beam_mesh_pt[0]->nelement();
  unsigned n_element2 = Beam_mesh_pt[1]->nelement();

  // Resize and initilize the Mapping which is the relationships between the
  // external data label l and the mesh g, element e, local node n
  // The order of the column is g,e,n,l (can check the document for specific
  // information)
  Mapping.resize(4, 2 * (n_element1 + n_element2), 0);

  unsigned column_label = 0;

  // Loop over every local node in all elements and arms and add them as
  // external data because they affect the traction and therefore the total drag
  // and torque on the object

  // Loop over the all arms
  unsigned n_mesh = beam_mesh_pt.size();
  for (unsigned g = 0; g < n_mesh; g++)
  {
    // Find number of elements in the mesh
    unsigned n_element = Beam_mesh_pt[g]->nelement();

    // Loop over the elements in this mesh
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast to the specific element type
      UnsteadyHaoHermiteBeamElement* elem_pt =
        dynamic_cast<UnsteadyHaoHermiteBeamElement*>(
          Beam_mesh_pt[g]->element_pt(e));

      // Find out how many nodes within this element
      const unsigned n_node = elem_pt->nnode();

      // Loop over the local nodes
      for (unsigned n = 0; n < n_node; n++)
      {
        // Since variable_position_pt() is only the member function of
        // SolidNode, thus transform the type of node into SolidNode
        SolidNode* solid_node_pt =
          dynamic_cast<SolidNode*>(elem_pt->node_pt(n));

        // Assemble the Mapping matrix
        Mapping(0, column_label) = g;
        Mapping(1, column_label) = e;
        Mapping(2, column_label) = n;

        // An error may be reported if a node has already been added, but it can
        // still return the index of the external data and skip this node
        Mapping(3, column_label) =
          add_external_data(solid_node_pt->variable_position_pt());

        // Move to the next column of the Mapping
        column_label++;
      }
    }
  }

  // Output the information of Mapping
  std::cout << "number of rows of the Mapping: " << Mapping.nrow() << std::endl;
  std::cout << "number of columns of the Mapping: " << Mapping.ncol()
            << std::endl;

  // Output Mapping
  std::cout << "Mapping: " << std::endl;
  for (unsigned i = 0; i < Mapping.nrow(); i++)
  {
    // Output the i-th row entries
    std::cout << "[";
    for (unsigned j = 0; j < Mapping.ncol(); j++)
    {
      std::cout << Mapping(i, j) << ", ";
    }
    std::cout << "]" << std::endl;
  }
}


//=============================================================================
/// Compute the beam's centre of mass (defined outside class to avoid
/// forward references)
//=============================================================================
void UnsteadyRigidBodyElement::compute_centre_of_mass(
  Vector<double>& sum_r_centre)
{
#ifdef PARANOID
  if (sum_r_centre.size() != 2)
  {
    std::ostringstream error_message;
    error_message << "sum_r_centre should have size 2, not "
                  << sum_r_centre.size() << std::endl;

    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Initialise
  sum_r_centre[0] = 0.0;
  sum_r_centre[1] = 0.0;
  Vector<double> int_r(2);
  double length = 0.0;

  // Find number of beam meshes
  unsigned npointer = Beam_mesh_pt.size();

  // Loop over the beam meshes to compute the centre of mass of the
  // entire beam
  for (unsigned i = 0; i < npointer; i++)
  {
    // Initialise
    Vector<double> total_int_r(2);
    double total_length = 0.0;

    // Find number of elements in the mesh
    unsigned n_element = Beam_mesh_pt[i]->nelement();

    // Loop over the elements to compute the sum of elements'
    // contribution to the (\int r ds) and the length of beam
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast to the specific element type
      UnsteadyHaoHermiteBeamElement* elem_pt =
        dynamic_cast<UnsteadyHaoHermiteBeamElement*>(
          Beam_mesh_pt[i]->element_pt(e));

      // Compute contribution to the the (\int r ds) and length of
      // beam within the e-th element
      elem_pt->compute_contribution_to_int_r_and_length(int_r, length);

      // Sum the elements' contribution to the (\int r ds) and length
      // of beam
      total_int_r[0] += int_r[0];
      total_int_r[1] += int_r[1];
      total_length += length;
    } // end of loop over elements

    // Assemble the (\int r ds) and beam length to get the centre of
    // mass for one arm
    Vector<double> r_centre(2);
    r_centre[0] = (1.0 / total_length) * total_int_r[0];
    r_centre[1] = (1.0 / total_length) * total_int_r[1];

    // Compute the centre of mass of the entire beam
    sum_r_centre[0] = sum_r_centre[0] + r_centre[0];
    sum_r_centre[1] = sum_r_centre[1] + r_centre[1];
  }

  // Get the centre of mass
  sum_r_centre[0] = 0.5 * sum_r_centre[0];
  sum_r_centre[1] = 0.5 * sum_r_centre[1];
}


//=============================================================================
/// Compute the drag and torque on the entire beam structure according
/// to slender body theory (Type=0: first arm, Type=1: second arm.)
//=============================================================================
void UnsteadyRigidBodyElement::compute_drag_and_torque(
  Vector<double>& sum_total_drag, double& sum_total_torque)
{
#ifdef PARANOID
  if (sum_total_drag.size() != 2)
  {
    std::ostringstream error_message;
    error_message << "sum_total_drag should have size 2, not "
                  << sum_total_drag.size() << std::endl;

    throw OomphLibError(
      error_message.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
  }
#endif

  // Initialise
  sum_total_drag[0] = 0.0;
  sum_total_drag[1] = 0.0;
  sum_total_torque = 0.0;

  Vector<double> drag(2);
  double torque = 0.0;

  // Find number of beam meshes
  unsigned npointer = Beam_mesh_pt.size();

  // Loop over the beam meshes to compute the drag and torque of the
  // entire beam
  for (unsigned i = 0; i < npointer; i++)
  {
    // Initialise
    Vector<double> total_drag(2);
    double total_torque = 0.0;

    // Find number of elements in the mesh
    unsigned n_element = Beam_mesh_pt[i]->nelement();

    // Loop over the elements to compute the sum of elements'
    // contribution to the drag and torque on the entire beam
    // structure
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast to the specific element type
      UnsteadyHaoHermiteBeamElement* elem_pt =
        dynamic_cast<UnsteadyHaoHermiteBeamElement*>(
          Beam_mesh_pt[i]->element_pt(e));

      // Compute contribution to the drag and torque within the e-th
      // element
      elem_pt->compute_contribution_to_drag_and_torque(drag, torque);

      // Sum the elements' contribution to the drag and torque
      total_drag[0] += drag[0];
      total_drag[1] += drag[1];
      total_torque += torque;
    } // end of loop over elements

    // Compute the drag and torque of the entire beam
    sum_total_drag[0] = sum_total_drag[0] + total_drag[0];
    sum_total_drag[1] = sum_total_drag[1] + total_drag[1];
    sum_total_torque = sum_total_torque + total_torque;
  }
}


//=============================================================================
/// Assemble the contributions to the jacobian and mass matrices
//=============================================================================
void UnsteadyRigidBodyElement::fill_in_contribution_to_jacobian_and_mass_matrix(
  Vector<double>& residuals,
  DenseMatrix<double>& jacobian,
  DenseMatrix<double>& mass_matrix)
{
  // Get the jacobian
  fill_in_contribution_to_jacobian(residuals, jacobian);

  // Mass matrix is evaluated at base state

  // Get the steady variables
  // Translate rigid body parameters into meaningful variables
  double V = 0.0;
  double U0 = 0.0;
  double theta_eq = 0.0;
  double X0 = 0.0;
  double Y0 = 0.0;
  Steady_rigid_body_element_pt->get_parameters(V, U0, theta_eq, X0, Y0);

  std::cout << "theta_eq(fluid)=" << theta_eq << std::endl;

  // Find number of beam meshes
  unsigned n_mesh = Beam_mesh_pt.size();

  // Integers to store the local equation and unknown numbers
  int local_eqn = 0;
  int local_unknown = 0;

  // Local coordinate (1D!)
  Vector<double> s(1);

  // Find number of internal_data in the mesh
  unsigned n_internal = ninternal_data();

  // Assemble the entries of d(drag,torque)/d(three motion unknowns)

  // Loop over the beam meshes
  for (unsigned g = 0; g < n_mesh; g++)
  {
    // Find number of elements in the mesh
    unsigned n_element = Beam_mesh_pt[g]->nelement();

    // Loop over the elements
    for (unsigned e = 0; e < n_element; e++)
    {
      // Upcast to the specific element type
      UnsteadyHaoHermiteBeamElement* elem_pt =
        dynamic_cast<UnsteadyHaoHermiteBeamElement*>(
          Beam_mesh_pt[g]->element_pt(e));

      // Set the number of integration points
      unsigned n_intpt = elem_pt->integral_pt()->nweight();

      // Find out how many nodes within this element
      const unsigned n_node = elem_pt->nnode();

      // Find out how many positional dofs within this element
      const unsigned n_position_type = elem_pt->nnodal_position_type();

      // # of nodes, # of positional dofs
      Shape psi(n_node, n_position_type);

      // Loop over the integration points
      for (unsigned ipt = 0; ipt < n_intpt; ipt++)
      {
        // Get the integral weight
        double w = elem_pt->integral_pt()->weight(ipt);

        // Return local coordinate s[j]  of i-th integration point.
        unsigned j = 0;
        s[j] = elem_pt->integral_pt()->knot(ipt, j);

        // Get shape functions
        elem_pt->shape(s, psi);

        // Get position vector to and non-unit tangent vector on wall:
        // dr/ds. NOTE: This is before we apply the rigid body motion!
        // so in terms of the write-up the position vector is R_0
        Vector<double> R_0(2);
        Vector<double> drds(2);
        elem_pt->get_non_unit_tangent(s, R_0, drds);

        // Jacobian of mapping between local and global coordinates
        double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

        // Premultiply the weights and the Jacobian
        double W = w * J;

        // Get the Eulerian position and the unit normal.
        // NOTE: This is before we apply the rigid body motion!
        // so in terms of the write-up the position vector is R_0 and N_0
        Vector<double> N_0(2);
        elem_pt->get_normal(s, R_0, N_0);

        // Compute normal N which is after translation and rotation
        Vector<double> N(2);
        N[0] = cos(theta_eq + elem_pt->theta_initial()) * N_0[0] -
               sin(theta_eq + elem_pt->theta_initial()) * N_0[1];
        N[1] = sin(theta_eq + elem_pt->theta_initial()) * N_0[0] +
               cos(theta_eq + elem_pt->theta_initial()) * N_0[1];

        // hierher!
        // Evalutate at t=0 since this is in steady state
        double t = 0.0;

        // Steady state for X and Y
        double steady_X = 0.5 * V * t * t + U0 * t + X0;
        double steady_Y = V * t + Y0;

        // hierher!
        // Compute R after translation and rotation
        Vector<double> R(2);
        R[0] = cos(theta_eq + elem_pt->theta_initial()) * R_0[0] -
               sin(theta_eq + elem_pt->theta_initial()) * R_0[1] + steady_X;
        R[1] = sin(theta_eq + elem_pt->theta_initial()) * R_0[0] +
               cos(theta_eq + elem_pt->theta_initial()) * R_0[1] + steady_Y;

        // Orientaion plus theta_initial (first arm: theta_initial=0.0, second
        // arm: theta_initial=alpha)
        double Theta_plus_theta_initial = theta_eq + elem_pt->theta_initial();

        // The first derivative of fluid traction with respect to
        // unknowns (element level)
        Vector<double> dtraction_f_ddotX(2);
        Vector<double> dtraction_f_ddotY(2);
        Vector<double> dtraction_f_ddotTheta(2);
        Vector<double> dtraction_f_ddotR_0_1(2);
        Vector<double> dtraction_f_ddotR_0_2(2);

        // Compute the first derivative of fluid traction with respect to
        // unknowns (element level)
        // traction_f is the fulid traction computed from the slender body
        // theory
        // Input: R_0 is the deformed position before rotation and translation,
        // N is the normal to the deformed configuration,
        // Theta_add_theta_initial is the orientaion plus theta_initial (first
        // arm: theta_initial=0.0, second arm: theta_initial=alpha)
        // Output:
        // dtraction_f_ddotX,dtraction_f_ddotY,dtraction_f_ddotTheta,dtraction_f_ddotR_0_1,dtraction_f_ddotR_0_2
        Global_Physical_Variables::compute_dtraction_f_ddot_unknowns(
          R_0,
          N,
          Theta_plus_theta_initial,
          dtraction_f_ddotX,
          dtraction_f_ddotY,
          dtraction_f_ddotTheta,
          dtraction_f_ddotR_0_1,
          dtraction_f_ddotR_0_2);


        // Loop over drag, torque
        for (unsigned m = 0; m < n_internal; m++)
        {
          // Find the equation number
          local_eqn = internal_local_eqn(m, 0);

          /*IF it's not pinned*/
          if (local_eqn >= 0)
          {
            // Loop over the internal data
            for (unsigned m2 = 0; m2 < n_internal; m2++)
            {
              local_unknown = internal_local_eqn(m2, 0);

              // If at a non-zero degree of freedom add in the entry
              if (local_unknown >= 0)
              {
                // Note all entries are added a minus sign compared with
                // Maple's expressions

                // Drag in x direction
                if (m == 0)
                {
                  // Coeff of dX0/dt
                  if (m2 == 0)
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -dtraction_f_ddotX[0] * W;
                  }
                  // Coeff of dY0/dt
                  else if (m2 == 1)
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -dtraction_f_ddotY[0] * W;
                  }
                  // Coeff of dtheta/dt
                  else
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -dtraction_f_ddotTheta[0] * W;
                  }
                }
                // Drag in y direction
                else if (m == 1)
                {
                  // Coeff of dX0/dt
                  if (m2 == 0)
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -dtraction_f_ddotX[1] * W;
                  }
                  // Coeff of dY0/dt
                  else if (m2 == 1)
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -dtraction_f_ddotY[1] * W;
                  }
                  // Coeff of dtheta/dt
                  else
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -dtraction_f_ddotTheta[1] * W;
                  }
                }
                // Torque
                else
                {
                  // Coeff of dX0/dt
                  if (m2 == 0)
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -(((R[0] - steady_X) * dtraction_f_ddotX[1] -
                         (R[1] - steady_Y) * dtraction_f_ddotX[0]) *
                        W);
                  }
                  // Coeff of dY0/dt
                  else if (m2 == 1)
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -(((R[0] - steady_X) * dtraction_f_ddotY[1] -
                         (R[1] - steady_Y) * dtraction_f_ddotY[0]) *
                        W);
                  }
                  // Coeff of dtheta/dt
                  else
                  {
                    mass_matrix(local_eqn, local_unknown) +=
                      -(((R[0] - steady_X) * dtraction_f_ddotTheta[1] -
                         (R[1] - steady_Y) * dtraction_f_ddotTheta[0]) *
                        W);
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  // Assemble the entries of d(drag,torque)/dR_{ijk}

  // Loop over the drag, torque
  for (unsigned m = 0; m < n_internal; m++)
  {
    // Find the equation number
    local_eqn = internal_local_eqn(m, 0);

    /*IF it's not pinned*/
    if (local_eqn >= 0)
    {
      // Loop over two beam meshes (arms)
      for (unsigned g = 0; g < n_mesh; g++)
      {
        // Find number of elements in the mesh
        unsigned n_element = Beam_mesh_pt[g]->nelement();

        // Loop over the elements
        for (unsigned e = 0; e < n_element; e++)
        {
          // Upcast to the specific element type
          UnsteadyHaoHermiteBeamElement* elem_pt =
            dynamic_cast<UnsteadyHaoHermiteBeamElement*>(
              Beam_mesh_pt[g]->element_pt(e));

          // Set the number of integration points
          unsigned n_intpt = elem_pt->integral_pt()->nweight();

          // Find out how many nodes within this element
          const unsigned n_node = elem_pt->nnode();

          // Find out how many positional dofs within this element
          const unsigned n_position_type = elem_pt->nnodal_position_type();

          // Set the dimension of the global coordinates
          unsigned n_dim = 2;

          // # of nodes, # of positional dofs
          Shape psi(n_node, n_position_type);


          // Loop over the integration points
          for (unsigned ipt = 0; ipt < n_intpt; ipt++)
          {
            // Get the integral weight
            double w = elem_pt->integral_pt()->weight(ipt);

            // Return local coordinate s[j]  of i-th integration point.
            unsigned j = 0;
            s[j] = elem_pt->integral_pt()->knot(ipt, j);

            // Get shape functions
            elem_pt->shape(s, psi);

            // Get position vector to and non-unit tangent vector on wall:
            // dr/ds. NOTE: This is before we apply the rigid body motion!
            // so in terms of the write-up the position vector is R_0
            Vector<double> R_0(2);
            Vector<double> drds(2);
            elem_pt->get_non_unit_tangent(s, R_0, drds);

            // Jacobian of mapping between local and global coordinates
            double J = sqrt(drds[0] * drds[0] + drds[1] * drds[1]);

            // Premultiply the weights and the Jacobian
            double W = w * J;

            // Get the Eulerian position and the unit normal.
            // NOTE: This is before we apply the rigid body motion!
            // so in terms of the write-up the position vector is R_0 and
            // N_0
            Vector<double> N_0(2);
            elem_pt->get_normal(s, R_0, N_0);

            // Compute normal N which is after translation and rotation
            Vector<double> N(2);
            N[0] = cos(theta_eq + elem_pt->theta_initial()) * N_0[0] -
                   sin(theta_eq + elem_pt->theta_initial()) * N_0[1];
            N[1] = sin(theta_eq + elem_pt->theta_initial()) * N_0[0] +
                   cos(theta_eq + elem_pt->theta_initial()) * N_0[1];

            // herher!
            // Evalutate at t=0 since this is in steady state
            double t = 0.0;

            // Steady state for X and Y
            double steady_X = 0.5 * V * t * t + U0 * t + X0;
            double steady_Y = V * t + Y0;

            // herher!
            // Compute R after translation and rotation
            Vector<double> R(2);
            R[0] = cos(theta_eq + elem_pt->theta_initial()) * R_0[0] -
                   sin(theta_eq + elem_pt->theta_initial()) * R_0[1] + steady_X;
            R[1] = sin(theta_eq + elem_pt->theta_initial()) * R_0[0] +
                   cos(theta_eq + elem_pt->theta_initial()) * R_0[1] + steady_Y;

            // Orientaion plus theta_initial (first arm: theta_initial=0.0,
            // second arm: theta_initial=alpha)
            double Theta_plus_theta_initial =
              theta_eq + elem_pt->theta_initial();

            // The first derivative of fluid traction with respect to
            // unknowns (element level)
            Vector<double> dtraction_f_ddotX(2);
            Vector<double> dtraction_f_ddotY(2);
            Vector<double> dtraction_f_ddotTheta(2);
            Vector<double> dtraction_f_ddotR_0_1(2);
            Vector<double> dtraction_f_ddotR_0_2(2);

            // Compute the first derivative of fluid traction with respect to
            // unknowns (element level)
            // traction_f is the fulid traction computed from the slender body
            // theory
            // Input: R_0 is the deformed position before rotation and
            // translation, N is the normal to the deformed configuration,
            // Theta_add_theta_initial is the orientaion plus theta_initial
            // (first arm: theta_initial=0.0, second arm: theta_initial=alpha)
            // Output:
            // dtraction_f_ddotX,dtraction_f_ddotY,dtraction_f_ddotTheta,dtraction_f_ddotR_0_1,dtraction_f_ddotR_0_2
            Global_Physical_Variables::compute_dtraction_f_ddot_unknowns(
              R_0,
              N,
              Theta_plus_theta_initial,
              dtraction_f_ddotX,
              dtraction_f_ddotY,
              dtraction_f_ddotTheta,
              dtraction_f_ddotR_0_1,
              dtraction_f_ddotR_0_2);


            // Assemble the contributions to the mass matrix
            // Loop over the number of nodes
            for (unsigned n = 0; n < n_node; n++)
            {
              // Loop over the coordinate directions
              for (unsigned i = 0; i < n_dim; i++)
              {
                // Loop over the type of position
                for (unsigned k = 0; k < n_position_type; k++)
                {
                  // Initialise external data label
                  unsigned l = 0;

                  // Initialise label of values of external data
                  unsigned v = 0;

                  // Find the external data index l and its value index v in
                  // terms of mesh g, element e, local node n, position
                  // direction i, and position type k
                  external_data_index_and_value_index(g, e, n, i, k, l, v);

                  // Find the equation number
                  local_unknown = external_local_eqn(l, v);

                  // If at a non-zero degree of freedom add in the entry
                  if (local_unknown >= 0)
                  {
                    // Drag[1] (x direction)
                    if (m == 0)
                    {
                      // Coeff of dR_1nk/dt (x direction)
                      if (i == 0)
                      {
                        mass_matrix(local_eqn, local_unknown) +=
                          -(dtraction_f_ddotR_0_1[0] * psi(n, k) * W);
                      }
                      // Coeff of dR_2nk/dt (y direction)
                      else
                      {
                        mass_matrix(local_eqn, local_unknown) +=
                          -(dtraction_f_ddotR_0_2[0] * psi(n, k) * W);
                      }
                    }
                    // Drag[2] (y direction)
                    else if (m == 1)
                    { // Coeff of dR_1nk/dt (x direction)
                      if (i == 0)
                      {
                        mass_matrix(local_eqn, local_unknown) +=
                          -(dtraction_f_ddotR_0_1[1] * psi(n, k) * W);
                      }
                      // Coeff of dR_2nk/dt (y direction)
                      else
                      {
                        mass_matrix(local_eqn, local_unknown) +=
                          -(dtraction_f_ddotR_0_2[1] * psi(n, k) * W);
                      }
                    }
                    // Torque
                    else
                    { // Coeff of dR_1nk/dt (x direction)
                      if (i == 0)
                      {
                        mass_matrix(local_eqn, local_unknown) +=
                          -(((R[0] - steady_X) * dtraction_f_ddotR_0_1[1] *
                               psi(n, k) -
                             (R[1] - steady_Y) * dtraction_f_ddotR_0_1[0] *
                               psi(n, k)) *
                            W);
                      }
                      // Coeff of dR_2nk/dt (y direction)
                      else
                      {
                        mass_matrix(local_eqn, local_unknown) +=
                          -(((R[0] - steady_X) * dtraction_f_ddotR_0_2[1] *
                               psi(n, k) -
                             (R[1] - steady_Y) * dtraction_f_ddotR_0_2[0] *
                               psi(n, k)) *
                            W);
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }


  // Record the jacobian matrices of every step
  std::cout << "number of rows of the motion-element jacobian: "
            << jacobian.nrow() << std::endl;
  std::cout << "number of columns of the motion-element jacobian: "
            << jacobian.ncol() << std::endl;

  // Output jacobian
  std::cout << "motion-element jacobian: " << std::endl;
  for (unsigned i = 0; i < jacobian.nrow(); i++)
  {
    // Output the i-th row entries
    std::cout << "[";
    for (unsigned j = 0; j < jacobian.ncol(); j++)
    {
      std::cout << jacobian(i, j) << ", ";
    }
    std::cout << "]" << std::endl;
  }

  // Output mass matrix
  std::cout << "motion-element mass matrix: " << std::endl;
  for (unsigned i = 0; i < mass_matrix.nrow(); i++)
  {
    // Output the i-th row entries
    std::cout << "[";
    for (unsigned j = 0; j < mass_matrix.ncol(); j++)
    {
      std::cout << mass_matrix(i, j) << ", ";
    }
    std::cout << "]" << std::endl;
  }

} // end_of_fill_in_contribution_to_jacobian_and_mass_matrix


//======start_of_problem_class==========================================
/// Beam problem object
//======================================================================
class UnsteadyElasticBeamProblem : public Problem
{
public:
  /// Constructor: The arguments are the number of elements and the
  /// parameter to determine the length of the beam
  UnsteadyElasticBeamProblem(const unsigned& n_elem1,
                             const unsigned& n_elem2,
                             const std::string& restart_steady_file);

  /// No actions need to be performed after a solve
  void actions_after_newton_solve() {}

  /// No actions need to be performed before a solve
  void actions_before_newton_solve() {}

  /// Update the problem specs after solve (empty)
  void actions_after_implicit_timestep() {}

  /// Update the problem specs before next timestep:
  /// Set Dirchlet boundary conditions from exact solution.
  void actions_before_implicit_timestep() {}

  /// Set initial condition (incl previous timesteps) according
  /// to specified function.
  void set_initial_condition();

  /// Doc the solution
  void doc_solution(DocInfo& doc_info, ofstream& trace_file);

  /// Dump problem data to allow for later restart
  void dump_it(ofstream& dump_file)
  {
    // Current value of the FSI parameter
    dump_file << Global_Physical_Variables::I << " # FSI parameter"
              << std::endl;

    // hierher maybe add q and alpha but then issue warning
    // if it differs from the one specified on the command line

    // Dump the refinement pattern and the generic problem data
    Problem::dump(dump_file);
  }

  /// Read problem data for restart
  void restart(ifstream& restart_file)
  {
    // Read line up to termination sign
    string input_string;
    getline(restart_file, input_string, '#');

    // Ignore rest of line
    restart_file.ignore(80, '\n');

    // Read in FSI parameter
    Global_Physical_Variables::I = double(atof(input_string.c_str()));

    // Refine the mesh and read in the generic problem data
    Problem::read(restart_file);
  }

  /// Global error norm for adaptive time-stepping
  double global_temporal_error_norm();

  /// Solve the problem
  void solve(const unsigned& label);

  /// Doc the solution, pass the number of the case considered,
  /// so that output files can be distinguished.
  void doc_solution(const unsigned& label);


  /// Get the fully assembled residual vector and mass matrix
  /// in dense storage. The DoubleVector residuals returned will be
  /// non-distributed. If on calling this method the DoubleVector residuals is
  /// setup then it must be non-distributed and of the correct length.
  /// The matrix type DenseDoubleMatrix is not distributable and therefore
  /// the residual vector is also assumed to be non distributable.
  void get_mass_matrix(DoubleVector& residuals, DenseDoubleMatrix& mass_matrix);


private:
  /// Pointer to geometric object that represents the beam's
  /// undeformed shape
  GeomObject* Undef_beam_pt1;

  GeomObject* Undef_beam_pt2;

  /// Pointer to UnsteadyRigidBodyElement that actually contains the
  /// rigid body data
  UnsteadyRigidBodyElement* Rigid_body_element_pt;

  /// Pointer to beam mesh (first arm)
  OneDLagrangianMesh<UnsteadyHaoHermiteBeamElement>* Beam_mesh_first_arm_pt;

  /// Pointer to beam mesh (second arm)
  OneDLagrangianMesh<UnsteadyHaoHermiteBeamElement>* Beam_mesh_second_arm_pt;

  /// Pointer to mesh containing the rigid body element
  Mesh* Rigid_body_element_mesh_pt;

  /// Pointer to the SteadyElasticBeamProblem
  ElasticBeamProblem* SteadyElasticBeamProblem_pt;

}; // end of problem class


/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////


//=========================================================================
/// Steady, straight 1D line in 2D space
///  \f[ x = 0.0 \f]
///  \f[ y = \zeta \f]
//=========================================================================
class StraightLineVertical : public GeomObject
{
public:
  /// Constructor derives from GeomObject(1, 2)
  StraightLineVertical() : GeomObject(1, 2) {}

  /// Broken copy constructor
  StraightLineVertical(const StraightLineVertical& dummy) = delete;

  /// Broken assignment operator
  void operator=(const StraightLineVertical&) = delete;

  /// Position Vector at Lagrangian coordinate zeta
  void position(const Vector<double>& zeta, Vector<double>& r) const
  {
    r[0] = 0.0;
    r[1] = zeta[0];
  }


  /// Derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// Evaluated at current time.
  virtual void dposition(const Vector<double>& zeta,
                         DenseMatrix<double>& drdzeta) const
  {
    // Tangent vector
    drdzeta(0, 0) = 0.0;
    drdzeta(0, 1) = 1.0;
  }


  /// 2nd derivative of position Vector w.r.t. to coordinates:
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i). Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }


  /// Posn Vector and its  1st & 2nd derivatives
  /// w.r.t. to coordinates:
  /// \f$ \frac{dR_i}{d \zeta_\alpha}\f$ = drdzeta(alpha,i).
  /// \f$ \frac{d^2R_i}{d \zeta_\alpha d \zeta_\beta}\f$ =
  /// ddrdzeta(alpha,beta,i).
  /// Evaluated at current time.
  virtual void d2position(const Vector<double>& zeta,
                          Vector<double>& r,
                          DenseMatrix<double>& drdzeta,
                          RankThreeTensor<double>& ddrdzeta) const
  {
    // Position Vector
    r[0] = 0.0;
    r[1] = zeta[0];

    // Tangent vector
    drdzeta(0, 0) = 0.0;
    drdzeta(0, 1) = 1.0;

    // Derivative of tangent vector
    ddrdzeta(0, 0, 0) = 0.0;
    ddrdzeta(0, 0, 1) = 0.0;
  }
};


//=============start_of_constructor=====================================
/// Constructor for elastic beam problem
//======================================================================
UnsteadyElasticBeamProblem::UnsteadyElasticBeamProblem(
  const unsigned& n_elem1,
  const unsigned& n_elem2,
  const std::string& restart_steady_file)
{
  Problem::Always_take_one_newton_step = true;

  /////////////////////////////////////////////////////////////////////////////////////////
  // Get the steady state data
  /////////////////////////////////////////////////////////////////////////////////////////
  unsigned n1 = 20;
  unsigned n2 = 20;
  SteadyElasticBeamProblem_pt = new ElasticBeamProblem(n1, n2);

  RigidBodyElement* steady_rigid_body_element_pt = nullptr;
  SolidMesh* steady_beam_mesh_first_arm_pt = nullptr;
  SolidMesh* steady_beam_mesh_second_arm_pt = nullptr;

  // Get the pointer from the steady problem
  SteadyElasticBeamProblem_pt->get_steady_problem_rigid_body_element_pt(
    steady_rigid_body_element_pt);

  SteadyElasticBeamProblem_pt->get_steady_problem_beam_meshes(
    steady_beam_mesh_first_arm_pt, steady_beam_mesh_second_arm_pt);

  /*   // Do the restart?
    // This is to read the steady state data from the disk
    if
    (CommandLineArgs::command_line_flag_has_been_set("--restart_steady_file"))
    {
      // Open/read restart file
      std::ifstream file2;
      file2.open(restart_steady_file.c_str());
      SteadyElasticBeamProblem_pt->restart(file2);
      file2.close();
    } */

  // delete SteadyElasticBeamProblem_pt;
  /////////////////////////////////////////////////////////////////////////////////////////////

  // Create the eigen solver
  this->eigen_solver_pt() = new LAPACK_QZ();

  // Allocate the timestepper -- this constructs the Problem's
  // time object with a sufficient amount of storage to store the
  // previous timsteps.
  add_time_stepper_pt(new Steady<1>); // hierher  BDF<1>(true));

  // x position of clamped point
  double x = 0.0;

  // y position of clamped point
  double y = 0.0;

  // Initialize the orientaion
  double theta = 0.0;

  // Make the UnsteadyRigidBodyElement that stores the parameters for
  // the rigid body motion
  Rigid_body_element_pt =
    new UnsteadyRigidBodyElement(x, y, theta, time_stepper_pt());

  // Add the Unsteady rigid body element to its own mesh
  Rigid_body_element_mesh_pt = new Mesh;
  Rigid_body_element_mesh_pt->add_element_pt(Rigid_body_element_pt);

  // Still use the same expression of q as before to represent the
  // length of the two arms
  double* stretch_ratio_pt = &Global_Physical_Variables::Stretch_ratio;
  Undef_beam_pt1 = new StraightLineVertical_new(fabs(*stretch_ratio_pt + 0.5));
  Undef_beam_pt2 = new StraightLineVertical_new(fabs(*stretch_ratio_pt - 0.5));

  // Create the (Lagrangian!) mesh, using the StraightLineVertical
  // object Undef_beam_pt to specify the initial (Eulerian) position
  // of the nodes. (first arm) double length_1 = fabs(*q_pt + 0.5);
  double length_1 = 1.0;
  Beam_mesh_first_arm_pt =
    new OneDLagrangianMesh<UnsteadyHaoHermiteBeamElement>(
      n_elem1, length_1, Undef_beam_pt1, time_stepper_pt());

  // Create the (Lagrangian!) mesh, using the StraightLineVertical
  // object Undef_beam_pt to specify the initial (Eulerian) position
  // of the nodes. (second arm) double length_2 = fabs(*q_pt - 0.5);
  double length_2 = 1.0;
  Beam_mesh_second_arm_pt =
    new OneDLagrangianMesh<UnsteadyHaoHermiteBeamElement>(
      n_elem2, length_2, Undef_beam_pt2, time_stepper_pt());

  // Pass the pointer of the mesh to the UnsteadyRigidBodyElement
  // class so it can work out the drag and torque on the entire
  // structure
  Vector<SolidMesh*> Beam_mesh_pt(2);
  Beam_mesh_pt[0] = Beam_mesh_first_arm_pt;
  Beam_mesh_pt[1] = Beam_mesh_second_arm_pt;
  Rigid_body_element_pt->set_pointer_to_beam_meshes(Beam_mesh_pt);

  // Pass the steady state pointer of the mesh to the
  // UnsteadyRigidBodyElement class
  Vector<SolidMesh*> steady_beam_mesh_pt(2);
  steady_beam_mesh_pt[0] = steady_beam_mesh_first_arm_pt;
  steady_beam_mesh_pt[1] = steady_beam_mesh_second_arm_pt;
  Rigid_body_element_pt
    ->set_pointer_to_steady_beam_meshes_and_steady_rigid_body_element(
      steady_beam_mesh_pt, steady_rigid_body_element_pt);

  // Build the problem's global mesh
  // Meshes and solve of steady state are built on
  // SteadyElasticBeamProblem_pt, thus here it does not need to flush
  // the submeshes and rebuild the global mesh
  add_sub_mesh(Beam_mesh_first_arm_pt);
  add_sub_mesh(Beam_mesh_second_arm_pt);
  add_sub_mesh(Rigid_body_element_mesh_pt);
  build_global_mesh();

  // Set the boundary conditions: One end of the beam is clamped in
  // space Pin displacements in both x and y directions, and pin the
  // derivative of position Vector w.r.t. to coordinates in x
  // direction. (first arm)
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(0);
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(1);
  Beam_mesh_first_arm_pt->boundary_node_pt(0, 0)->pin_position(1, 0);

  /*   // Loop over nodes in mesh and pinned them
    n_node = Beam_mesh_first_arm_pt->nnode();
    for (unsigned j = 0; j < n_node; j++)
    {
      // Pin
      Beam_mesh_first_arm_pt->node_pt(j)->pin_position(0);
      Beam_mesh_first_arm_pt->node_pt(j)->pin_position(1);
      Beam_mesh_first_arm_pt->node_pt(j)->pin_position(1, 0);
      Beam_mesh_first_arm_pt->node_pt(j)->pin_position(1, 1);
    } */

  // Find number of elements in the mesh (first arm)
  unsigned n_element = Beam_mesh_first_arm_pt->nelement();

  // Loop over the elements to set physical parameters etc. (first
  // arm)
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast to the specific element type
    UnsteadyHaoHermiteBeamElement* elem_pt =
      dynamic_cast<UnsteadyHaoHermiteBeamElement*>(
        Beam_mesh_first_arm_pt->element_pt(e));

    // Pass the pointer of UnsteadyRigidBodyElement to the each
    // element so we can work out the rigid body motion
    elem_pt->set_pointer_to_rigid_body_element(Rigid_body_element_pt);

    // Set physical parameters for each element:
    elem_pt->h_pt() = &Global_Physical_Variables::H;
    elem_pt->i_pt() = &Global_Physical_Variables::I;
    elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

    // Pass the pointer of steady RigidBodyElement to the each element
    elem_pt->set_pointer_to_steady_rigid_body_element(
      steady_rigid_body_element_pt);

    // Note: no rotation!

    // Set the undeformed shape for each element
    elem_pt->undeformed_beam_pt() = Undef_beam_pt1;

  } // end of loop over elements

  // Set the boundary conditions: One end of the beam is clamped in
  // space Pin displacements in both x and y directions, and pin the
  // derivative of position Vector w.r.t. to coordinates in x
  // direction. (second arm)
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(0);
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(1);
  Beam_mesh_second_arm_pt->boundary_node_pt(0, 0)->pin_position(1, 0);

  /*   // Loop over nodes in mesh and pinned them
    n_node = Beam_mesh_second_arm_pt->nnode();
    for (unsigned j = 0; j < n_node; j++)
    {
      // Pin
      Beam_mesh_second_arm_pt->node_pt(j)->pin_position(0);
      Beam_mesh_second_arm_pt->node_pt(j)->pin_position(1);
      Beam_mesh_second_arm_pt->node_pt(j)->pin_position(1, 0);
      Beam_mesh_second_arm_pt->node_pt(j)->pin_position(1, 1);
    } */

  // Find number of elements in the mesh (second arm)
  n_element = Beam_mesh_second_arm_pt->nelement();

  // Loop over the elements to set physical parameters etc. (second
  // arm)
  for (unsigned e = 0; e < n_element; e++)
  {
    // Upcast to the specific element type
    UnsteadyHaoHermiteBeamElement* elem_pt =
      dynamic_cast<UnsteadyHaoHermiteBeamElement*>(
        Beam_mesh_second_arm_pt->element_pt(e));

    // Pass the pointer of UnsteadyRigidBodyElement to the each
    // element so we can work out the rigid body motion
    elem_pt->set_pointer_to_rigid_body_element(Rigid_body_element_pt);

    // Set physical parameters for each element:
    elem_pt->h_pt() = &Global_Physical_Variables::H;
    elem_pt->i_pt() = &Global_Physical_Variables::I;
    elem_pt->lambda_sq_pt() = &Global_Physical_Variables::Lambda_sq;

    // Pass the pointer of steady RigidBodyElement to the each element
    elem_pt->set_pointer_to_steady_rigid_body_element(
      steady_rigid_body_element_pt);

    // Rotate by opening angle
    elem_pt->theta_initial_pt(&Global_Physical_Variables::Alpha);

    // Set the undeformed shape for each element
    elem_pt->undeformed_beam_pt() = Undef_beam_pt2;

  } // end of loop over elements

  // Assign the global and local equation numbers
  cout << "# of dofs " << assign_eqn_numbers() << std::endl;


  // Create label for output
  DocInfo doc_info;

  // Set output directory -- this function checks if the output
  // directory exists and issues a warning if it doesn't.
  doc_info.set_directory("RESLT");

  // Output solution
  //-----------------
  char filename[100];

  // Output file stream used for writing results
  ofstream file;

  // Write the file name
  sprintf(filename, "RESLT/trace.dat");
  file.open(filename);

  // Counter to record the iterations for the while loop
  unsigned counter = 0;

  double ds = 0.0;

  while (Global_Physical_Variables::I >= 0.00)
  {
    if (counter == 0)
    {
      // Solve the system
      SteadyElasticBeamProblem_pt->newton_solve();
    }
    else
    {
      // To prevent large solution jumps in critical intervals for I, try to
      // reduce ds specifically in those areas for smoother and more precise
      // results

      // First I interval [Interval1_start,Interval1_end]
      if (Global_Physical_Variables::Ds_default >
            Global_Physical_Variables::Ds_interval1 &&
          Global_Physical_Variables::I >=
            Global_Physical_Variables::Interval1_start &&
          Global_Physical_Variables::I <=
            Global_Physical_Variables::Interval1_end)
      {
        ds = Global_Physical_Variables::Ds_interval1;
      }
      // Second I interval [Interval2_start,Interval2_end]
      else if (Global_Physical_Variables::Ds_default >
                 Global_Physical_Variables::Ds_interval2 &&
               Global_Physical_Variables::I >=
                 Global_Physical_Variables::Interval2_start &&
               Global_Physical_Variables::I <=
                 Global_Physical_Variables::Interval2_end)
      {
        ds = Global_Physical_Variables::Ds_interval2;
      }
      else
      {
        // Use the default one
        ds = Global_Physical_Variables::Ds_default;
      }

      /// Use the arclength solve
      ds = SteadyElasticBeamProblem_pt->arc_length_step_solve(
        &Global_Physical_Variables::I, ds);
    }

    // Document I
    file << Global_Physical_Variables::I << "  ";

    // Document the solution of Theta_eq, Theta_eq_orientation
    steady_rigid_body_element_pt->output(file);

    // Step label
    file << counter << "  ";

    // Loop over nodes in mesh and assign steady position values to this
    // unsteady mesh nodes
    unsigned n_node = Beam_mesh_first_arm_pt->nnode();
    for (unsigned j = 0; j < n_node; j++)
    {
      // Position
      Beam_mesh_first_arm_pt->node_pt(j)->x_gen(0, 0) =
        steady_beam_mesh_first_arm_pt->node_pt(j)->x_gen(0, 0);
      Beam_mesh_first_arm_pt->node_pt(j)->x_gen(0, 1) =
        steady_beam_mesh_first_arm_pt->node_pt(j)->x_gen(0, 1);

      // Slope
      Beam_mesh_first_arm_pt->node_pt(j)->x_gen(1, 0) =
        steady_beam_mesh_first_arm_pt->node_pt(j)->x_gen(1, 0);
      Beam_mesh_first_arm_pt->node_pt(j)->x_gen(1, 1) =
        steady_beam_mesh_first_arm_pt->node_pt(j)->x_gen(1, 1);
    }

    // Loop over nodes in mesh and assign steady position values to this
    // unsteady mesh nodes
    n_node = Beam_mesh_second_arm_pt->nnode();
    for (unsigned j = 0; j < n_node; j++)
    {
      // Position
      Beam_mesh_second_arm_pt->node_pt(j)->x_gen(0, 0) =
        steady_beam_mesh_second_arm_pt->node_pt(j)->x_gen(0, 0);
      Beam_mesh_second_arm_pt->node_pt(j)->x_gen(0, 1) =
        steady_beam_mesh_second_arm_pt->node_pt(j)->x_gen(0, 1);

      // Slope
      Beam_mesh_second_arm_pt->node_pt(j)->x_gen(1, 0) =
        steady_beam_mesh_second_arm_pt->node_pt(j)->x_gen(1, 0);
      Beam_mesh_second_arm_pt->node_pt(j)->x_gen(1, 1) =
        steady_beam_mesh_second_arm_pt->node_pt(j)->x_gen(1, 1);
    }

    // Copy the steady orientaion to unsteady problem
    double theta_eq =
      steady_rigid_body_element_pt->internal_data_pt(2)->value(0);
    Rigid_body_element_pt->internal_data_pt(2)->set_value(0, theta_eq);

    // Reinitilize X and Y
    Rigid_body_element_pt->internal_data_pt(0)->set_value(0, 0.0);
    Rigid_body_element_pt->internal_data_pt(1)->set_value(0, 0.0);

    // Output the steady configuration
    ofstream some_file1;
    sprintf(filename,
            "%s/beam_first_arm_%i.dat",
            doc_info.directory().c_str(),
            counter);
    some_file1.open(filename);
    Beam_mesh_first_arm_pt->output(some_file1);
    some_file1.close();

    ofstream some_file2;
    sprintf(filename,
            "%s/beam_second_arm_%i.dat",
            doc_info.directory().c_str(),
            counter);
    some_file2.open(filename);
    Beam_mesh_second_arm_pt->output(some_file2);
    some_file2.close();

    // Output the steady configuration
    ofstream some_file5;
    sprintf(filename,
            "%s/steady_beam_first_arm_%i.dat",
            doc_info.directory().c_str(),
            counter);
    some_file5.open(filename);
    steady_beam_mesh_first_arm_pt->output(some_file5);
    some_file5.close();

    ofstream some_file6;
    sprintf(filename,
            "%s/steady_beam_second_arm_%i.dat",
            doc_info.directory().c_str(),
            counter);
    some_file6.open(filename);
    steady_beam_mesh_second_arm_pt->output(some_file6);
    some_file6.close();

    // Translate rigid body parameters into meaningful variables
    double X = 0.0;
    double Y = 0.0;
    double Theta = 0.0;
    Rigid_body_element_pt->get_parameters(X, Y, Theta);

    // Output X
    file << X << "  ";

    // Output Y
    file << Y << "  ";

    // Output Theta
    file << Theta << "  ";

    solve(counter);

    // Translate rigid body parameters into meaningful variables
    Rigid_body_element_pt->get_parameters(X, Y, Theta);

    // Output X
    file << X << "  ";

    // Output Y
    file << Y << "  ";

    // Output Theta
    file << Theta << std::endl;

    // Bump counter for output
    counter = counter + 1;

    // exit(0);
  }

  file.close();

  DoubleVector residuals;
  DenseDoubleMatrix jacobian;
  DenseDoubleMatrix mass_matrix;
  get_jacobian(residuals, jacobian);
  get_mass_matrix(residuals, mass_matrix);

  // Output the full jacobian
  ofstream some_file3;
  sprintf(filename, "RESLT/full_jacobian.dat");
  some_file3.open(filename);
  jacobian.sparse_indexed_output(some_file3);
  some_file3.close();

  // Output the full mass matrix
  ofstream some_file4;
  sprintf(filename, "RESLT/full_mass_matrix.dat");
  some_file4.open(filename);
  mass_matrix.sparse_indexed_output(some_file4);
  some_file4.close();

} // end of constructor


//=====start_of_set_ic=====================================================
/// Setup initial conditions -- either restart from solution
/// specified via command line or impulsive start.
//=========================================================================
void UnsteadyElasticBeamProblem::set_initial_condition()
{
  // Assign impulsive start
  assign_initial_values_impulsive();

} // end of set_initial_condition


//=======start_of_doc_solution============================================
/// Doc the solution
//========================================================================
void UnsteadyElasticBeamProblem::doc_solution(DocInfo& doc_info,
                                              ofstream& trace_file)
{
  ofstream some_file1;
  char filename[100];

  // Number of plot points
  unsigned npts;
  npts = 5;


  cout << std::endl;
  cout << "=================================================" << std::endl;
  cout << "Docing solution for t=" << time_pt()->time() << std::endl;
  cout << "=================================================" << std::endl;


  // Output solution
  //-----------------
  sprintf(filename,
          "%s/soln_first_arm_%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file1.open(filename);
  Beam_mesh_first_arm_pt->output(some_file1, npts);
  some_file1.close();

  ofstream some_file2;
  sprintf(filename,
          "%s/soln_second_arm_%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file2.open(filename);
  Beam_mesh_second_arm_pt->output(some_file2, npts);
  some_file2.close();

  trace_file << Global_Physical_Variables::I << "  ";
  Rigid_body_element_pt->output(trace_file);
  trace_file << time_pt()->time() << std::endl;

  // Write restart file
  ofstream some_file3;
  sprintf(filename,
          "%s/restart%i.dat",
          doc_info.directory().c_str(),
          doc_info.number());
  some_file3.open(filename);
  dump_it(some_file3);
  some_file3.close();

} // end of doc_solution


//========start_of_global_temporal_error_norm==============================
/// Global error norm for adaptive timestepping: RMS error, based on
/// difference between predicted and actual value at all nodes.
//=========================================================================
double UnsteadyElasticBeamProblem::global_temporal_error_norm()
{
  double global_error = 0.0;

  // Rigid_body_element_mesh only contains one node
  unsigned n_node = 3;

  // Calculate the estimated error in the values

  // Get error in solution: Difference between predicted and actual
  // value for X0(t), Y0(t), orientation theta(t)
  for (unsigned i = 0; i < n_node; i++)
  {
    double error = Rigid_body_element_pt->internal_data_pt(i)
                     ->time_stepper_pt()
                     ->temporal_error_in_value(
                       Rigid_body_element_pt->internal_data_pt(i), 0);

    // Add the square of the individual error to the global error
    global_error += error * error;
  }

  // Divide by the number of nodes
  global_error /= double(n_node);

  // Return square root...
  return sqrt(global_error);

} // end of global_temporal_error_norm


//===start_of_doc=========================================================
/// Doc the solution in tecplot format. Label files with label.
//========================================================================
void UnsteadyElasticBeamProblem::doc_solution(const unsigned& label)
{
  ofstream some_file;
  char filename[100];

  // Output solution with specified number of plot points per element
  sprintf(filename, "RESLT/soln%i.dat", label);
  some_file.open(filename);
  Rigid_body_element_pt->output_eigenvectors(some_file);
  Beam_mesh_first_arm_pt->output(some_file);
  Beam_mesh_second_arm_pt->output(some_file);
  some_file.close();

} // end of doc

//=======================start_of_solve==============================
/// Solve the eigenproblem
//===================================================================
void UnsteadyElasticBeamProblem::solve(const unsigned& label)
{
  // Output the entire jacobian
  DoubleVector residuals;
  DenseDoubleMatrix jacobian;
  residuals.initialise(0.0);
  jacobian.initialise(0.0);
  get_jacobian(residuals, jacobian);
  std::cout << "number of rows of the entire jacobian: " << jacobian.nrow()
            << std::endl;
  std::cout << "number of columns of the entire jacobian: " << jacobian.ncol()
            << std::endl;

  std::cout << "entire_jacobian1: " << std::endl;
  for (unsigned i = 0; i < jacobian.nrow(); i++)
  {
    // Output the i-th row entries
    std::cout << "[";
    for (unsigned j = 0; j < jacobian.ncol(); j++)
    {
      std::cout << jacobian(i, j) << ", ";
    }
    std::cout << "]" << std::endl;
  }

  // Set external storage for the eigenvalues
  Vector<complex<double>> eigenvalue;
  // Set external storage for the eigenvectors
  Vector<DoubleVector> eigenvector_real;
  Vector<DoubleVector> eigenvector_imag;
  // Desired number eigenvalues
  unsigned n_eval = this->ndof();

  // Solve the eigenproblem
  this->solve_eigenproblem(
    n_eval, eigenvalue, eigenvector_real, eigenvector_imag);

  // We now need to sort the output based on the size of the real part
  // of the eigenvalues.
  // This is because the solver does not necessarily sort the
  // eigenvalues
  Vector<complex<double>> sorted_eigenvalue = eigenvalue;
  sort(
    sorted_eigenvalue.begin(), sorted_eigenvalue.end(), ComplexLess<double>());


  // Read out the largest eigenvalue
  complex<double> temp_evalue = sorted_eigenvalue[n_eval - 1];
  unsigned largest_index = 0;
  // Loop over the unsorted eigenvalues and find the entry that corresponds
  // to our largest eigenvalue.
  for (unsigned i = 0; i < eigenvalue.size(); i++)
  {
    // Note that equality tests for doubles are bad, but it was just
    // sorted data, so should be fine
    if (eigenvalue[i] == temp_evalue)
    {
      largest_index = i;
      break;
    }
  }

  // Normalise the eigenvector
  {
    // Get the dimension of the eigenvector
    unsigned dim = eigenvector_real[largest_index].nrow();
    double length = 0.0;
    // Loop over all the entries
    for (unsigned i = 0; i < dim; i++)
    {
      // Add the contribution to the length
      length += std::pow(eigenvector_real[largest_index][i], 2.0);
    }
    // Now take the magnitude
    length = sqrt(length);
    // Fix the sign
    if (eigenvector_real[largest_index][0] < 0)
    {
      length *= -1.0;
    }
    // Finally normalise
    for (unsigned i = 0; i < dim; i++)
    {
      eigenvector_real[largest_index][i] /= length;
    }
  }

  // The magnitude of the perturbation
  double beta = 1.0e-3;

  // Now assign the eigenvector to the dofs of the problem
  this->add_eigenvector_to_dofs(beta, eigenvector_real[largest_index]);
  // Output solution for this case (label output files with "j")
  this->doc_solution(label);

  char filename[100];
  sprintf(filename, "RESLT/eigenvalues%i.dat", label);

  // Open an output file for the sorted eigenvalues
  ofstream evalues(filename);
  for (unsigned i = 0; i < n_eval; i++)
  {
    // Print to screen
    cout << eigenvalue[i].real() << " " << eigenvalue[i].imag() << std::endl;
    // Send to file
    evalues << eigenvalue[i].real() << " " << eigenvalue[i].imag() << std::endl;
  }

  evalues.close();

  // Output the entire jacobian
  residuals.initialise(0.0);
  jacobian.initialise(0.0);
  get_jacobian(residuals, jacobian);
  std::cout << "number of rows of the entire jacobian: " << jacobian.nrow()
            << std::endl;
  std::cout << "number of columns of the entire jacobian: " << jacobian.ncol()
            << std::endl;

  std::cout << "entire_jacobian2: " << std::endl;
  for (unsigned i = 0; i < jacobian.nrow(); i++)
  {
    // Output the i-th row entries
    std::cout << "[";
    for (unsigned j = 0; j < jacobian.ncol(); j++)
    {
      std::cout << jacobian(i, j) << ", ";
    }
    std::cout << "]" << std::endl;
  }

} // end_of_solve


//=============================================================================
/// Get the fully assembled residual vector and mass matrix
/// in dense storage. The DoubleVector residuals returned will be
/// non-distributed. If on calling this method the DoubleVector residuals is
/// setup then it must be non-distributed and of the correct length.
/// The matrix type DenseDoubleMatrix is not distributable and therefore
/// the residual vector is also assumed to be non distributable.
//=============================================================================
void UnsteadyElasticBeamProblem::get_mass_matrix(DoubleVector& residuals,
                                                 DenseDoubleMatrix& mass_matrix)
{
  // get the number of degrees of freedom
  unsigned n_dof = ndof();

#ifdef PARANOID
  // PARANOID checks : if the distribution of residuals is setup then it must
  // must not be distributed, have the right number of rows, and the same
  // communicator as the problem
  if (residuals.built())
  {
    if (residuals.distribution_pt()->distributed())
    {
      std::ostringstream error_stream;
      error_stream << "If the DoubleVector residuals is setup then it must not "
                   << "be distributed.";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (residuals.distribution_pt()->nrow() != n_dof)
    {
      std::ostringstream error_stream;
      error_stream << "If the DoubleVector residuals is setup then it must have"
                   << " the correct number of rows";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
    if (!(*Communicator_pt == *residuals.distribution_pt()->communicator_pt()))
    {
      std::ostringstream error_stream;
      error_stream << "If the DoubleVector residuals is setup then it must have"
                   << " the same communicator as the problem.";
      throw OomphLibError(
        error_stream.str(), OOMPH_CURRENT_FUNCTION, OOMPH_EXCEPTION_LOCATION);
    }
  }
#endif

  // set the residuals distribution if it is not setup
  if (!residuals.built())
  {
    LinearAlgebraDistribution dist(Communicator_pt, n_dof, false);
    residuals.build(&dist, 0.0);
  }
  // else just zero the residuals
  else
  {
    residuals.initialise(0.0);
  }

  // Resize the matrices -- this cannot always be done externally
  // because get_mass_matrix exists in many different versions for
  // different storage formats -- resizing a CC or CR matrix doesn't
  // make sense.

  // resize the mass matrix
  mass_matrix.resize(n_dof, n_dof);
  mass_matrix.initialise(0.0);

  DenseDoubleMatrix jacobian;
  // resize the jacobian
  jacobian.resize(n_dof, n_dof);
  jacobian.initialise(0.0);

  // Locally cache pointer to assembly handler
  AssemblyHandler* const assembly_handler_pt = this->assembly_handler_pt();

  // Loop over all the elements
  unsigned long n_element = mesh_pt()->nelement();
  for (unsigned long e = 0; e < n_element; e++)
  {
    // Get the pointer to the element
    GeneralisedElement* elem_pt = mesh_pt()->element_pt(e);

    // Find number of dofs in the element
    unsigned n_element_dofs = assembly_handler_pt->ndof(elem_pt);

    // Set up an array
    Vector<double> element_residuals(n_element_dofs);

    // Set up a matrix
    DenseMatrix<double> element_jacobian(n_element_dofs);

    // Set up a matrix
    DenseMatrix<double> element_mass_matrix(n_element_dofs);

    // Fill the array
    elem_pt->get_jacobian_and_mass_matrix(
      element_residuals, element_jacobian, element_mass_matrix);

    // Now loop over the dofs and assign values to global Vector
    for (unsigned l = 0; l < n_element_dofs; l++)
    {
      unsigned long eqn_number = assembly_handler_pt->eqn_number(elem_pt, l);
      residuals[eqn_number] += element_residuals[l];

      for (unsigned l2 = 0; l2 < n_element_dofs; l2++)
      {
        mass_matrix(eqn_number, assembly_handler_pt->eqn_number(elem_pt, l2)) +=
          element_mass_matrix(l, l2);
      }
    }
  }
}


//========start_of_main================================================
/// Driver for beam (string under tension) test problem
//=====================================================================
int main(int argc, char** argv)
{
  // Store command line arguments
  CommandLineArgs::setup(argc, argv);

  // Stretch_ratio
  CommandLineArgs::specify_command_line_flag(
    "--q", &Global_Physical_Variables::Stretch_ratio);

  // Aspect ratio
  // CommandLineArgs::specify_command_line_flag("--q",
  //&Global_Physical_Variables::Q);

  // Opening angle in degrees
  double alpha_in_degrees = 45.0;
  CommandLineArgs::specify_command_line_flag("--alpha_in_degrees",
                                             &alpha_in_degrees);

  // Initial value for theta_eq in the Newton solve
  CommandLineArgs::specify_command_line_flag(
    "--Initial_value_for_theta_eq",
    &Global_Physical_Variables::Initial_value_for_theta_eq);

  // Initial value for theta_eq in the Newton solve
  CommandLineArgs::specify_command_line_flag(
    "--ds_default", &Global_Physical_Variables::Ds_default);

  // End point for the first I interval
  // [Interval1_start,Interval1_end]
  CommandLineArgs::specify_command_line_flag(
    "--interval1_start", &Global_Physical_Variables::Interval1_start);

  CommandLineArgs::specify_command_line_flag(
    "--interval1_end", &Global_Physical_Variables::Interval1_end);

  // Start and end points for the second I interval
  // [Interval2_start,Interval2_end]
  CommandLineArgs::specify_command_line_flag(
    "--interval2_start", &Global_Physical_Variables::Interval2_start);

  CommandLineArgs::specify_command_line_flag(
    "--interval2_end", &Global_Physical_Variables::Interval2_end);

  // Value of ds for first interval
  CommandLineArgs::specify_command_line_flag(
    "--ds_interval1", &Global_Physical_Variables::Ds_interval1);

  // Value of ds for second interval
  CommandLineArgs::specify_command_line_flag(
    "--ds_interval2", &Global_Physical_Variables::Ds_interval2);

  // Restart file
  std::string restart_file;
  CommandLineArgs::specify_command_line_flag("--restart_file", &restart_file);

  // Restart file
  std::string restart_steady_file;
  CommandLineArgs::specify_command_line_flag("--restart_steady_file",
                                             &restart_steady_file);

  // Parse command line
  CommandLineArgs::parse_and_assign();

  // Doc what has actually been specified on the command line
  CommandLineArgs::doc_specified_flags();

  // Now that we've read the opening angle in degrees, update the
  // value in radians
  Global_Physical_Variables::Alpha = 4.0 * atan(1.0) / 180.0 * alpha_in_degrees;

  // Set the non-dimensional thickness
  Global_Physical_Variables::H = 0.01;

  // Number of elements (choose an even number if you want the control
  // point to be located at the centre of the beam)
  unsigned n_element1 = 20;
  unsigned n_element2 = 20;

  // Assign the FSI parameter
  Global_Physical_Variables::I = 0.0;

  // Construct the problem
  UnsteadyElasticBeamProblem problem(
    n_element1, n_element2, restart_steady_file);

  // Setup labels for output
  // DocInfo doc_info;

  // Output directory
  // doc_info.set_directory("RESLT");

  // Output number
  // doc_info.number() = 0;

  /* // Open a trace file
  ofstream trace_file;
  char filename[100];
  sprintf(filename, "%s/trace.dat", doc_info.directory().c_str());
  trace_file.open(filename); */

  // Choose simulation interval
  // double t_max = 10.0;

  // Initial timestep: Use the one used when setting up the initial
  // condition
  double dt = 1.0;

  // Initialise timestep -- also sets the weights for all timesteppers
  // in the problem.
  problem.initialise_dt(dt);

  // Do the restart?
  if (CommandLineArgs::command_line_flag_has_been_set("--restart_file"))
  {
    // Open/read restart file
    std::ifstream file2;
    file2.open(restart_file.c_str());
    problem.restart(file2);
    file2.close();
  }
  else
  {
    // Set IC
    problem.set_initial_condition();
  }

  problem.describe_dofs();

  // exit(0);

  // Output initial condition
  // problem.doc_solution(doc_info, trace_file);

  // Increment counter for solutions
  // doc_info.number()++;

  // Target error for adaptive timestepping
  // double epsilon_t = 1.0e-4;

  // try
  //{
  // problem.steady_newton_solve();
  //}
  // catch (...)
  //{
  // Eigenvalues
  // Solve with LAPACK_QZ
  // problem.solve(1);
  //}

  // Eigenvalues
  // Solve with LAPACK_QZ
  // problem.solve(1);

  /*  // Timestepping loop: Don't know how many steps we're going to
   take
   // in advance
   while (problem.time_pt()->time() < t_max)
   {
     // Take an adaptive timestep -- the input dt is the suggested
   timestep.
     // The adaptive timestepper will adjust dt until the required
   error
     // tolerance is satisfied. The function returns a suggestion
     // for the timestep that should be taken for the next step. This
     // is either the actual timestep taken this time or a larger
     // value if the solution was found to be "too accurate".
     double dt_next = problem.adaptive_unsteady_newton_solve(dt,
   epsilon_t);

     // Use dt_next as suggestion for the next timestep
     dt = dt_next;

     // Output solution
     problem.doc_solution(doc_info, trace_file);

     // Increment counter for solutions
     doc_info.number()++;

   } // end of timestepping loop

   // Close trace file
   trace_file.close(); */


} // end of main
