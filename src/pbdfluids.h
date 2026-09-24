#pragma once

#include "LinearMath.h"
#include <cmath>
#include <algorithm>

using namespace LinearMath;

class CubicSpline
{
protected:
	static Real m_radius;
	static Real m_k;
	static Real m_l;
	static Real m_W_zero;
public:
	static Real getRadius() { return m_radius; }
	static void setRadius(Real val)
	{
		m_radius = val;
		const Real pi = static_cast<Real>(3.14159265358979323846);
		const Real h3 = m_radius * m_radius * m_radius;
		m_k = static_cast<Real>(8.0) / (pi * h3);
		m_l = static_cast<Real>(48.0) / (pi * h3);
		m_W_zero = W(Vector3r(0.0, 0.0, 0.0));
	}

	static Real W(const Vector3r &r)
	{
		Real res = 0.0;
		const Real rl = r.norm();
		const Real q = rl / m_radius;
		if (q <= static_cast<Real>(1.0))
		{
			if (q <= static_cast<Real>(0.5))
			{
				const Real q2 = q * q;
				const Real q3 = q2 * q;
				res = m_k * (static_cast<Real>(6.0) * q3 - static_cast<Real>(6.0) * q2 + static_cast<Real>(1.0));
			}
			else
			{
				const Real diff = static_cast<Real>(1.0) - q;
				res = m_k * (static_cast<Real>(2.0) * diff * diff * diff);
			}
		}
		return res;
	}

	static Real W_zero()
	{
		return m_W_zero;
	}

	static Vector3r gradW(const Vector3r &r)
	{
		Vector3r res;
		res.setZero();
		const Real rl = r.norm();
		const Real q = rl / m_radius;
		if (q <= static_cast<Real>(1.0))
		{
			if (rl > static_cast<Real>(1.0e-9))
			{
				const Vector3r gradq = r * (static_cast<Real>(1.0) / (rl * m_radius));
				if (q <= static_cast<Real>(0.5))
				{
					res = m_l * q * (static_cast<Real>(3.0) * q - static_cast<Real>(2.0)) * gradq;
				}
				else
				{
					const Real factor = static_cast<Real>(1.0) - q;
					res = m_l * (-factor * factor) * gradq;
				}
			}
		}
		return res;
	}
};

class pbdfluids
{
public:
	// -------------- Position Based Fluids  -----------------------------------------------------
		
	/** Perform an SPH computation of the density of a fluid particle: 
	* \f{equation*}{
	* \rho_i = \sum_j m_j W(\mathbf{x}_i-\mathbf{x}_j).		
	* \f}
	* An additional term is added for neighboring boundary particles 
	* according to \cite Akinci:2012 in order to perform boundary handling.\n\n
	* Remark: A neighboring particle with an index >= numberOfParticles is
	* handled as boundary particle.\n\n
	*
	* More information can be found in the following papers: \cite Macklin:2013:PBF, \cite BMOTM2014, \cite BMM2015
	*
	* @param particleIndex	index of current fluid particle
	* @param numberOfParticles	number of fluid particles
	* @param x	array of all particle positions
	* @param mass array of all particle masses
	* @param boundaryX array of all boundary particles
	* @param boundaryPsi array of all boundary psi values (see \cite Akinci:2012)
	* @param numNeighbors number of neighbors
	* @param neighbors array with indices of all neighbors (indices larger than numberOfParticles are boundary particles)
	* @param density0 rest density
	* @param boundaryHandling perform boundary handling (see \cite Akinci:2012)
	* @param density_err returns the clamped density error (can be used for enforcing a maximal global density error)
	* @param density return the density
	*/	
	static bool computePBFDensity(
		const unsigned int particleIndex,				// current fluid particle	
		const unsigned int numberOfParticles,			// number of fluid particles 
		const Vector3r x[],						// array of all particle positions
		const Real mass[],								// array of all particle masses
		const Vector3r boundaryX[],				// array of all boundary particles
		const Real boundaryPsi[],						// array of all boundary psi values (Akinci2012)
		const unsigned int numNeighbors,				// number of neighbors 
		const unsigned int neighbors[],					// array with indices of all neighbors (indices larger than numberOfParticles are boundary particles)
		const Real density0,							// rest density
		const bool boundaryHandling,					// perform boundary handling (Akinci2012)
		Real &density_err,								// returns the clamped density error (can be used for enforcing a maximal global density error)
		Real &density);								// return the density

	/**
		* Compute Lagrange multiplier \f$\lambda_i\f$ for a fluid particle which is required by
		* the solver step:
		* \f{equation*}{
		* \lambda_i = -\frac{C_i(\mathbf{x}_1,...,\mathbf{x}_n)}{\sum_k \|\nabla_{\mathbf{x}_k} C_i\|^2 + \varepsilon}
		* \f}
		* with the constraint gradient:
		* \f{equation*}{
		* \nabla_{\mathbf{x}_k}C_i = \frac{m_j}{\rho_0}
		* \begin{cases}
		* \sum\limits_j \nabla_{\mathbf{x}_k} W(\mathbf{x}_i-\mathbf{x}_j, h) & \text{if }  k = i \\
		* -\nabla_{\mathbf{x}_k} W(\mathbf{x}_i-\mathbf{x}_j, h)  & \text{if } k = j.
		* \end{cases}
		* \f}
		* \n
		* Remark: The computation of the gradient is extended for neighboring boundary particles
		* according to \cite Akinci:2012 to perform a boundary handling. A neighboring 
		* particle with an index >= numberOfParticles is handled as boundary particle.\n\n
		*
		* More information can be found in the following papers: \cite Macklin:2013:PBF, \cite BMOTM2014, \cite BMM2015
		*
		* @param particleIndex	index of current fluid particle
		* @param numberOfParticles	number of fluid particles
		* @param x	array of all particle positions
		* @param mass array of all particle masses
		* @param boundaryX array of all boundary particles
		* @param boundaryPsi array of all boundary psi values (see \cite Akinci:2012)
		* @param density density of current fluid particle
		* @param numNeighbors number of neighbors
		* @param neighbors array with indices of all neighbors (indices larger than numberOfParticles are boundary particles)
		* @param density0 rest density
		* @param boundaryHandling perform boundary handling (see \cite Akinci:2012)
		* @param lambda returns the Lagrange multiplier
		*/		
	static bool computePBFLagrangeMultiplier(
		const unsigned int particleIndex,				// current fluid particle	
		const unsigned int numberOfParticles,			// number of fluid particles 
		const Vector3r x[],						// array of all particle positions
		const Real mass[],								// array of all particle masses
		const Vector3r boundaryX[],				// array of all boundary particles
		const Real boundaryPsi[],						// array of all boundary psi values (Akinci2012)
		const Real density,							// density of current fluid particle
		const unsigned int numNeighbors,				// number of neighbors 
		const unsigned int neighbors[],					// array with indices of all neighbors
		const Real density0,							// rest density
		const bool boundaryHandling,					// perform boundary handling (Akinci2012)
		Real &lambda);									// returns the Lagrange multiplier

	/** Perform a solver step for a fluid particle:
	*
	* \f{equation*}{
	* \Delta\mathbf{x}_{i} = \frac{m_j}{\rho_0}\sum\limits_j{\left(\lambda_i + \lambda_j\right)\nabla W(\mathbf{x}_i-\mathbf{x}_j, h)},
	* \f}
	* where \f$h\f$ is the smoothing length of the kernel function \f$W\f$.\n
	*\n
	* Remark: The computation of the position correction is extended for neighboring boundary particles
	* according to \cite Akinci:2012 to perform a boundary handling. A neighboring
	* particle with an index >= numberOfParticles is handled as boundary particle.\n\n
	*
	* More information can be found in the following papers: \cite Macklin:2013:PBF, \cite BMOTM2014, \cite BMM2015
	*
	* @param particleIndex	index of current fluid particle
	* @param numberOfParticles	number of fluid particles
	* @param x	array of all particle positions
	* @param mass array of all particle masses
	* @param boundaryX array of all boundary particles
	* @param boundaryPsi array of all boundary psi values (see \cite Akinci:2012)
	* @param numNeighbors number of neighbors
	* @param neighbors array with indices of all neighbors (indices larger than numberOfParticles are boundary particles)
	* @param density0 rest density
	* @param boundaryHandling perform boundary handling (see \cite Akinci:2012)
	* @param lambda Lagrange multipliers
	* @param corr returns the position correction for the current fluid particle
	*/
	static bool solveDensityConstraint(
		const unsigned int particleIndex,				// current fluid particle	
		const unsigned int numberOfParticles,			// number of fluid particles 
		const Vector3r x[],						// array of all particle positions
		const Real mass[],								// array of all particle masses
		const Vector3r boundaryX[],				// array of all boundary particles
		const Real boundaryPsi[],						// array of all boundary psi values (Akinci2012)
		const unsigned int numNeighbors,				// number of neighbors 
		const unsigned int neighbors[],					// array with indices of all neighbors
		const Real density0,							// rest density
		const bool boundaryHandling,					// perform boundary handling (Akinci2012)
		const Real lambda[],							// Lagrange multiplier
		Vector3r &corr);							// returns the position correction for the current fluid particle
};