#include "pbdcosseratrods.h"
#include <cmath>

const Real eps = static_cast<Real>(1e-6);

bool pbdcosseratrods::solve_StretchShearConstraint(
	const Vector3r& p0, Real invMass0,
	const Vector3r& p1, Real invMass1,
	const Quaternionr& q0, Real invMassq0,
	const Vector3r& stretchingAndShearingKs,
	const Real restLength,
	Vector3r& corr0, Vector3r& corr1, Quaternionr& corrq0)
{
	Vector3r d3;	//third director d3 = q0 * e_3 * q0_conjugate
	d3[0] = static_cast<Real>(2.0) * (q0.x() * q0.z() + q0.w() * q0.y());
	d3[1] = static_cast<Real>(2.0) * (q0.y() * q0.z() - q0.w() * q0.x());
	d3[2] = q0.w() * q0.w() - q0.x() * q0.x() - q0.y() * q0.y() + q0.z() * q0.z();

	Vector3r gamma = (p1 - p0) / restLength - d3;
	gamma /= (invMass1 + invMass0) / restLength + invMassq0 * static_cast<Real>(4.0) * restLength + eps;

	if (std::abs(stretchingAndShearingKs[0] - stretchingAndShearingKs[1]) < eps && std::abs(stretchingAndShearingKs[0] - stretchingAndShearingKs[2]) < eps)	//all Ks are approx. equal
		for (int i = 0; i < 3; i++) gamma[i] *= stretchingAndShearingKs[i];
	else	//diffenent stretching and shearing Ks. Transform diag(Ks[0], Ks[1], Ks[2]) into world space using Ks_w = R(q0) * diag(Ks[0], Ks[1], Ks[2]) * R^T(q0) and multiply it with gamma
	{
		Matrix3r R = q0.toRotationMatrix();
		gamma = (R.transpose() * gamma).eval();
		for (int i = 0; i < 3; i++) gamma[i] *= stretchingAndShearingKs[i];
		gamma = (R * gamma).eval();
	}

	corr0 = invMass0 * gamma;
	corr1 = -invMass1 * gamma;

	Quaternionr q_e_3_bar(q0.z(), -q0.y(), q0.x(), -q0.w());	//compute q*e_3.conjugate (cheaper than quaternion product)
	corrq0 = Quaternionr(0.0, gamma.x(), gamma.y(), gamma.z()) * q_e_3_bar;
	corrq0.coeffs() *= static_cast<Real>(2.0) * invMassq0 * restLength;

	return true;
}

bool pbdcosseratrods::solve_BendTwistConstraint(
	const Quaternionr& q0, Real invMassq0,
	const Quaternionr& q1, Real invMassq1,
	const Vector3r& bendingAndTwistingKs,
	const Quaternionr& restDarbouxVector,
	Quaternionr& corrq0, Quaternionr& corrq1)
{
	Quaternionr omega = q0.conjugate() * q1;   //darboux vector

	Quaternionr omega_plus;
	omega_plus.coeffs() = omega.coeffs() + restDarbouxVector.coeffs();     //delta Omega with -Omega_0
	omega.coeffs() = omega.coeffs() - restDarbouxVector.coeffs();                 //delta Omega with + omega_0
	if (omega.squaredNorm() > omega_plus.squaredNorm()) omega = omega_plus;

	for (int i = 0; i < 3; i++) omega.coeffs()[i] *= bendingAndTwistingKs[i] / (invMassq0 + invMassq1 + static_cast<Real>(1.0e-6));
	omega.w() = 0.0;    //discrete Darboux vector does not have vanishing scalar part

	corrq0 = q1 * omega;
	corrq1 = q0 * omega;
	corrq0.coeffs() *= invMassq0;
	corrq1.coeffs() *= -invMassq1;
	return true;
}
