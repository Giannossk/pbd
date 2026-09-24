#pragma once

#include "LinearMath.h"
#include "rods.h"

using namespace LinearMath;

#include <vector>
#include <list>

// Implementation of "Position And Orientation Based Cosserat Rods" paper
// (https://animation.rwth-aachen.de/publication/0550/)
//
//	Implemented by:
//
//	Tassilo Kugelstadt
//	Computer Animation Group
//	RWTH Aachen University
//
//  kugelstadt[at] cs.rwth-aachen.de
//
class pbdcosseratrods
{
public:
	/** Determine the position and orientation corrections for the stretch and shear constraint constraint (eq. 37 in the paper). \n\n
	*
	* @param p0 position of first particle
	* @param invMass0 inverse mass of first particle
	* @param p1 position of second particle
	* @param invMass1 inverse mass of second particle
	* @param q0 Quaternionr at the center of the edge
	* @param invMassq0 inverse mass of the quaternion
	* @param stretchingAndShearingKs stiffness coefficients for stretching and shearing
	* @param restLength rest edge length
	* @param corr0 position correction of first particle
	* @param corr1 position correction of second particle
	* @param corrq0 orientation correction of quaternion
	*/
	static bool solve_StretchShearConstraint(
		const Vector3r& p0, Real invMass0,
		const Vector3r& p1, Real invMass1,
		const Quaternionr& q0, Real invMassq0,
		const Vector3r& stretchingAndShearingKs,
		const Real restLength,
		Vector3r& corr0, Vector3r& corr1, Quaternionr& corrq0);

	/** Determine the position corrections for the bending and torsion constraint constraint (eq. 40 in the paper). \n\n
	*
	* @param q0 first quaternion
	* @param invMassq0 inverse mass of the first quaternion
	* @param q1 second quaternion
	* @param invMassq1 inverse Mass of the second quaternion
	* @param bendingAndTwistingKs stiffness coefficients for stretching and shearing
	* @param restDarbouxVector rest Darboux vector
	* @param corrq0 position correction of first particle
	* @param corrq1 position correction of second particle
	*/
	static bool solve_BendTwistConstraint(
		const Quaternionr& q0, Real invMassq0,
		const Quaternionr& q1, Real invMassq1,
		const Vector3r& bendingAndTwistingKs,
		const Quaternionr& restDarbouxVector,
		Quaternionr& corrq0, Quaternionr& corrq1);
};