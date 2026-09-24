#include "pbdelasticrods.h"
#include "pbdstiffrods.h"

#define _USE_MATH_DEFINES
#include <cmath>

const Real eps = static_cast<Real>(1e-6);

const int permutation[3][3] = {
	0, 2, 1,
	1, 0, 2,
	2, 1, 0
};

void pbdstiffrods::initLists(int numberOfIntervals, std::list <Node*> * &forward, std::list <Node*> * &backward, Node* &root)
{
	if (forward != NULL)
		delete[] forward;
	if (backward != NULL)
		delete[] backward;
	if (root != NULL)
		delete[] root;
	forward = new std::list<Node*>[numberOfIntervals];
	backward = new std::list<Node*>[numberOfIntervals];
	root = new Node[numberOfIntervals];
}

bool pbdstiffrods::isSegmentInInterval(RodSegment *segment, int intervalIndex, Interval* intervals, std::vector<RodConstraint*> &rodConstraints, std::vector<RodSegment*> &rodSegments)
{
	for (int i = intervals[intervalIndex].start; i <= intervals[intervalIndex].end; i++)
	{
		if ((segment == rodSegments[rodConstraints[i]->segmentIndex(0)])
			|| (segment == rodSegments[rodConstraints[i]->segmentIndex(1)]))
			return true;
	}
	return false;
}

bool pbdstiffrods::isConstraintInInterval(RodConstraint *constraint, int intervalIndex, Interval* intervals, std::vector<RodConstraint*> &rodConstraints)
{
	for (int i = intervals[intervalIndex].start; i <= intervals[intervalIndex].end; i++)
	{
		if (constraint == rodConstraints[i])
			return true;
	}
	return false;
}

void pbdstiffrods::initSegmentNode(Node *n, int intervalIndex, std::vector<RodConstraint*> &rodConstraints, std::vector<RodSegment*> &rodSegments, std::vector <RodConstraint*> &markedConstraints, Interval* intervals)
{
	RodSegment *segment =
		(RodSegment*)n->object;

	std::vector<RodConstraint*> constraints;
	std::vector<int> constraintIndices;
	for (int j = 0; j < static_cast<int>(rodConstraints.size()); ++j)
	{
		RodConstraint* constraint(rodConstraints[j]);
		if (rodSegments[constraint->segmentIndex(0)] == segment
			|| rodSegments[constraint->segmentIndex(1)] == segment)
		{
			constraints.push_back(constraint);
			constraintIndices.push_back(j);
		}
	}

	for (unsigned int i = 0; i < constraints.size(); i++)
	{
		if (!isConstraintInInterval(constraints[i], intervalIndex, intervals, rodConstraints))
			continue;

		// Test whether the edge has been visited before
		bool marked = false;
		for (unsigned int j = 0; j < markedConstraints.size(); j++)
		{
			if (constraints[i] == markedConstraints[j])
			{
				marked = true;
				break;
			}
		}
		if (!marked)
		{
			Node *constraintNode = new Node();
			constraintNode->index = constraintIndices[i];
			constraintNode->object = constraints[i];
			constraintNode->isconstraint = true;
			constraintNode->parent = n;
			constraintNode->D.setZero();
			constraintNode->Dinv.setZero();
			constraintNode->J.setZero();
			constraintNode->soln.setZero();

			n->children.push_back(constraintNode);

			Node *segmentNode = new Node();
			segmentNode->isconstraint = false;
			segmentNode->parent = constraintNode;

			//	get other segment connected to constraint for new node
			if (rodSegments[constraints[i]->segmentIndex(0)] == segment)
			{
				segmentNode->object = rodSegments[constraints[i]->segmentIndex(1)];
				segmentNode->index = constraints[i]->segmentIndex(1);
			}
			else
			{
				segmentNode->object = rodSegments[constraints[i]->segmentIndex(0)];
				segmentNode->index = constraints[i]->segmentIndex(0);
			}

			segmentNode->D.setZero();
			segmentNode->Dinv.setZero();
			segmentNode->J.setZero();
			segmentNode->soln.setZero();

			constraintNode->children.push_back(segmentNode);

			// mark constraint
			markedConstraints.push_back(constraints[i]);

			initSegmentNode(segmentNode, intervalIndex, rodConstraints,
				rodSegments, markedConstraints, intervals);
		}
	}
}

void pbdstiffrods::orderMatrix(Node *n, int intervalIndex, std::list <Node*> * forward, std::list <Node*> * backward)
{
	for (unsigned int i = 0; i < n->children.size(); i++)
		orderMatrix(n->children[i], intervalIndex, forward, backward);
	forward[intervalIndex].push_back(n);
	backward[intervalIndex].push_front(n);
}

void pbdstiffrods::initNodes(int intervalIndex, std::vector<RodSegment*> &rodSegments, Node* &root, Interval* intervals, std::vector<RodConstraint*> &rodConstraints, std::list <Node*> * forward, std::list <Node*> * backward, std::vector <RodConstraint*> &markedConstraints)
{
	// find root
	for (int i = 0; i < (int)rodSegments.size(); i++)
	{
		RodSegment * rb(rodSegments[i]);
		if (!isSegmentInInterval(rb, intervalIndex, intervals, rodConstraints, rodSegments))
			continue;
		else
		{
			if (root[intervalIndex].object == NULL)
			{
				root[intervalIndex].object = rb;
				root[intervalIndex].index = i;
			}
		}

		if (!rb->isDynamic())
		{
			root[intervalIndex].object = rb;
			root[intervalIndex].index = i;
			break;
		}
	}
	root[intervalIndex].isconstraint = false;
	root[intervalIndex].parent = NULL;

	root[intervalIndex].D.setZero();
	root[intervalIndex].Dinv.setZero();
	root[intervalIndex].soln.setZero();

	initSegmentNode(&root[intervalIndex], intervalIndex, rodConstraints,
		rodSegments, markedConstraints, intervals);
	orderMatrix(&root[intervalIndex], intervalIndex, forward, backward);
}

void pbdstiffrods::initTree(std::vector<RodConstraint*> &rodConstraints, std::vector<RodSegment*> & rodSegments, Interval* &intervals, int &numberOfIntervals, std::list <Node*> * &forward, std::list <Node*> * &backward, Node* &root)
{
	numberOfIntervals = 1;
	intervals = new Interval[1];
	intervals[0].start = 0;
	intervals[0].end = (int)rodConstraints.size() - 1;
	initLists(numberOfIntervals, forward, backward, root);

	std::vector <RodConstraint*> markedConstraints;
	for (int i = 0; i < numberOfIntervals; i++)
	{
		initNodes(i, rodSegments, root, intervals, rodConstraints, forward, backward, markedConstraints);
		markedConstraints.clear();
	}
}

bool pbdstiffrods::computeDarbouxVector(const Quaternionr & q0, const Quaternionr & q1, const Real averageSegmentLength, Vector3r & darbouxVector)
{
	darbouxVector = 2. / averageSegmentLength * (q0.conjugate() * q1).vec();
	return true;
}

bool pbdstiffrods::computeBendingAndTorsionJacobians(const Quaternionr & q0, const Quaternionr & q1, const Real averageSegmentLength, Matrix<Real, 3, 4> & jOmega0, Matrix<Real, 3, 4> & jOmega1)
{
	jOmega0(0, 0) = -q1.w(); jOmega0(0, 1) = -q1.z(); jOmega0(0, 2) =  q1.y(); jOmega0(0, 3) =  q1.x();
	jOmega0(1, 0) =  q1.z(); jOmega0(1, 1) = -q1.w(); jOmega0(1, 2) = -q1.x(); jOmega0(1, 3) =  q1.y();
	jOmega0(2, 0) = -q1.y(); jOmega0(2, 1) =  q1.x(); jOmega0(2, 2) = -q1.w(); jOmega0(2, 3) =  q1.z();

	jOmega1(0, 0) =  q0.w(); jOmega1(0, 1) =  q0.z(); jOmega1(0, 2) = -q0.y(); jOmega1(0, 3) = -q0.x();
	jOmega1(1, 0) = -q0.z(); jOmega1(1, 1) =  q0.w(); jOmega1(1, 2) =  q0.x(); jOmega1(1, 3) = -q0.y();
	jOmega1(2, 0) =  q0.y(); jOmega1(2, 1) = -q0.x(); jOmega1(2, 2) =  q0.w(); jOmega1(2, 3) = -q0.z();

	jOmega0 *= static_cast<Real>(2.0) / averageSegmentLength;
	jOmega1 *= static_cast<Real>(2.0) / averageSegmentLength;
	return true;
}

bool pbdstiffrods::computeMatrixG(const Quaternionr & q, Matrix<Real, 4, 3> & G)
{
	// w component at index 3
	G(0, 0) =  static_cast<Real>(0.5)*q.w(); G(0, 1) =  static_cast<Real>(0.5)*q.z(); G(0, 2) = -static_cast<Real>(0.5)*q.y();
	G(1, 0) = -static_cast<Real>(0.5)*q.z(); G(1, 1) =  static_cast<Real>(0.5)*q.w(); G(1, 2) =  static_cast<Real>(0.5)*q.x();
	G(2, 0) =  static_cast<Real>(0.5)*q.y(); G(2, 1) = -static_cast<Real>(0.5)*q.x(); G(2, 2) =  static_cast<Real>(0.5)*q.w();
	G(3, 0) = -static_cast<Real>(0.5)*q.x(); G(3, 1) = -static_cast<Real>(0.5)*q.y(); G(3, 2) = -static_cast<Real>(0.5)*q.z();
	return true;
}

void pbdstiffrods::computeMatrixK(const Vector3r &connector, const Real invMass, const Vector3r &x, const Matrix3r &inertiaInverseW, Matrix3r &K)
{
	matrix::computeMatrixK(connector, invMass, x, inertiaInverseW, K);
}

void pbdstiffrods::getMassMatrix(RodSegment *segment, Matrix6r &M)
{
	if (!segment->isDynamic())
	{
		M = Matrix6r::Identity();
		return;
	}

	const Vector3r &inertiaLocal = segment->InertiaTensor();
	Matrix3r rotationMatrix(segment->Rotation().toRotationMatrix());
	Matrix3r scaledR;
	for (int i = 0; i < 3; ++i)
		for (int j = 0; j < 3; ++j)
			scaledR(i, j) = rotationMatrix(i, j) * inertiaLocal[j];
	Matrix3r inertia = scaledR * rotationMatrix.transpose();

	Real mass = segment->Mass();

	// Upper half
	for (int i = 0; i < 3; i++)
	for (int j = 0; j < 6; j++)
	if (i == j)
		M(i, j) = mass;
	else
		M(i, j) = 0.0;

	// lower left
	for (int i = 3; i < 6; i++)
	for (int j = 0; j < 3; j++)
		M(i, j) = 0.0;

	// lower right
	for (int i = 3; i < 6; i++)
	for (int j = 3; j < 6; j++)
		M(i, j) = inertia(i - 3, j - 3);
}

Real pbdstiffrods::factor(const int intervalIndex, const std::vector<RodConstraint*> &rodConstraints, std::vector<RodSegment*> & rodSegments, const Interval* &intervals, std::list <Node*> * forward, std::list <Node*> * backward, std::vector<Vector6r> & RHS, std::vector<Vector6r> & lambdaSums, std::vector<std::vector<Matrix3r>> & bendingAndTorsionJacobians)
{
	Real maxError(0.);
	
	// compute right hand side of linear equation system
	for (size_t currentConstraintIndex = 0; currentConstraintIndex < rodConstraints.size(); ++currentConstraintIndex)
	{
		RodConstraint* currentConstraint = rodConstraints[currentConstraintIndex];

		RodSegment* segment0 = rodSegments[currentConstraint->segmentIndex(0)];
		RodSegment* segment1 = rodSegments[currentConstraint->segmentIndex(1)];

		const Quaternionr &q0 = segment0->Rotation();
		const Quaternionr &q1 = segment1->Rotation();

		const Matrix<Real, 3, 4, DontAlign> &constraintInfo(currentConstraint->getConstraintInfo());
		Vector6r &rhs(RHS[currentConstraintIndex]);

		// Compute zero-stretch part of constraint violation
		const Vector3r &connector0 = constraintInfo.col(2);
		const Vector3r &connector1 = constraintInfo.col(3);
		Vector3r stretchViolation = connector0 - connector1;

		// compute Darboux vector (Equation (7))
		Vector3r omega;
		computeDarbouxVector(q0, q1, currentConstraint->getAverageSegmentLength(), omega);

		// Compute bending and torsion part of constraint violation
		Vector3r bendingAndTorsionViolation = omega - currentConstraint->getRestDarbouxVector();

		// fill right hand side of the linear equation system
		const Vector6r &lambdaSum(lambdaSums[currentConstraintIndex]);
		const Vector3r &stretchComp = currentConstraint->getStretchCompliance();
		rhs.block<3, 1>(0, 0) = -stretchViolation - Vector3r(
			stretchComp[0] * lambdaSum[0],
			stretchComp[1] * lambdaSum[1],
			stretchComp[2] * lambdaSum[2]);

		const Vector3r &btComp = currentConstraint->getBendingAndTorsionCompliance();
		rhs.block<3, 1>(3, 0) = -bendingAndTorsionViolation - Vector3r(
			btComp[0] * lambdaSum[3],
			btComp[1] * lambdaSum[4],
			btComp[2] * lambdaSum[5]);

		// compute max error
		for (unsigned char i(0); i < 6; ++i)
		{
			maxError = std::max(maxError, std::abs(rhs[i]));
		}

		// Compute a part of the Jacobian here, because the relationship
		// of the first and second segment to the constraint can be determined directly

		// compute G matrices
		Matrix<Real, 4, 3> G0, G1;
		computeMatrixG(q0, G0);
		computeMatrixG(q1, G1);

		// compute stretching bending Jacobians (Equation (10) and Equation (11))
		Matrix<Real, 3, 4> jOmega0, jOmega1;
		computeBendingAndTorsionJacobians(q0, q1, currentConstraint->getAverageSegmentLength(), jOmega0, jOmega1);

		bendingAndTorsionJacobians[currentConstraintIndex][0] = jOmega0*G0;
		bendingAndTorsionJacobians[currentConstraintIndex][1] = jOmega1*G1;
	}

	std::list<Node*>::iterator nodeIter;
	for (nodeIter = forward[intervalIndex].begin(); nodeIter != forward[intervalIndex].end(); nodeIter++)
	{
		Node *node = *nodeIter;
		// compute system matrix diagonal
		if (node->isconstraint)
		{
			RodConstraint* currentConstraint = (RodConstraint*)node->object;
			//insert compliance
			node->D.setZero();
			const Vector3r &stretchCompliance(currentConstraint->getStretchCompliance());

			node->D(0, 0) -= stretchCompliance[0];
			node->D(1, 1) -= stretchCompliance[1];
			node->D(2, 2) -= stretchCompliance[2];

			const Vector3r &bendingAndTorsionCompliance(currentConstraint->getBendingAndTorsionCompliance());
			node->D(3, 3) -= bendingAndTorsionCompliance[0];
			node->D(4, 4) -= bendingAndTorsionCompliance[1];
			node->D(5, 5) -= bendingAndTorsionCompliance[2];
		}
		else
		{
			getMassMatrix((RodSegment*)node->object, node->D);
		}

		// compute Jacobian
		if (node->parent != NULL)
		{
			if (node->isconstraint)
			{
				//compute J 
				RodConstraint *constraint = (RodConstraint*)node->object;
				RodSegment *segment = (RodSegment*)node->parent->object;

				Real sign = 1;
				int segmentIndex = 0;
				if (segment == rodSegments[constraint->segmentIndex(1)])
				{
					segmentIndex = 1;
					sign = -1;
				}

				const Matrix<Real, 3, 4, DontAlign> &constraintInfo(constraint->getConstraintInfo());
				const Vector3r r = constraintInfo.col(2 + segmentIndex) - segment->Position();
				Matrix3r r_cross;
				Real crossSign(-static_cast<Real>(1.0)*sign);
				matrix::crossProductMatrix(crossSign*r, r_cross);

				node->J.block<3, 3>(0, 0) = Matrix3r::Identity() * sign;

				Matrix3r lowerLeft(Matrix3r::Zero());
				node->J.block<3, 3>(3, 0) = lowerLeft;

				node->J.block<3, 3>(0, 3) = r_cross;

				Matrix3r &lowerRight(bendingAndTorsionJacobians[node->index][segmentIndex]);
				node->J.block<3, 3>(3, 3) = lowerRight;
			}
			else
			{
				//compute JT
				RodConstraint *constraint = (RodConstraint*)node->parent->object;
				RodSegment *segment = (RodSegment*)node->object;

				Real sign = 1;
				int segmentIndex = 0;
				if (segment == rodSegments[constraint->segmentIndex(1)])
				{
					segmentIndex = 1;
					sign = -1;
				}

				const Matrix<Real, 3, 4, DontAlign> &constraintInfo(constraint->getConstraintInfo());
				const Vector3r r = constraintInfo.col(2 + segmentIndex) - segment->Position();
				Matrix3r r_crossT;
				matrix::crossProductMatrix(sign*r, r_crossT);

				node->J.block<3, 3>(0, 0) = Matrix3r::Identity() * sign;

				node->J.block<3, 3>(3, 0) = r_crossT;

				Matrix3r upperRight(Matrix3r::Zero());
				node->J.block<3, 3>(0, 3) = upperRight;

				Matrix3r lowerRight(bendingAndTorsionJacobians[node->parent->index][segmentIndex].transpose());
				node->J.block<3, 3>(3, 3) = lowerRight;
			}
		}
	}

	for (nodeIter = forward[intervalIndex].begin(); nodeIter != forward[intervalIndex].end(); nodeIter++)
	{
		Node *node = *nodeIter;
		std::vector <Node*> children = node->children;
		for (size_t i = 0; i < children.size(); i++)
		{
			Matrix6r JT = (children[i]->J).transpose();
			Matrix6r &D = children[i]->D;
			Matrix6r &J = children[i]->J;
			Matrix6r JTDJ = ((JT * D) * J);
			node->D = node->D - JTDJ;
		}
		bool chk = false;
		if (!node->isconstraint)
		{
			RodSegment *segment = (RodSegment*)node->object;
			if (!segment->isDynamic())
			{
				node->Dinv.setZero();
				chk = true;
			}
		}

		node->DLDLT.compute(node->D); // result reused in solve()
		if (node->parent != NULL)
		{
			if (!chk)
			{
				node->J = node->DLDLT.solve(node->J);
			}
			else
			{
				node->J.setZero();
			}
		}
	}
	return maxError;
}

bool pbdstiffrods::solve(int intervalIndex, std::list <Node*> * forward, std::list <Node*> * backward, std::vector<Vector6r> & RHS, std::vector<Vector6r> & lambdaSums, std::vector<Vector3r> & corr_x, std::vector<Quaternionr> & corr_q)
{
	std::list<Node*>::iterator nodeIter;
	for (nodeIter = forward[intervalIndex].begin(); nodeIter != forward[intervalIndex].end(); nodeIter++)
	{
		Node *node = *nodeIter;
		if (node->isconstraint)
		{
			node->soln = -RHS[node->index];
		}
		else
		{
			node->soln.setZero();
		}
		std::vector <Node*> &children = node->children;
		for (size_t i = 0; i < children.size(); ++i)
		{
			Matrix6r cJT = children[i]->J.transpose();
			Vector6r &csoln = children[i]->soln;
			Vector6r v = cJT * csoln;
			node->soln = node->soln - v;
		}
	}

	for (nodeIter = backward[intervalIndex].begin(); nodeIter != backward[intervalIndex].end(); nodeIter++)
	{
		Node *node = *nodeIter;

		bool noZeroDinv(true);
		if (!node->isconstraint)
		{
			RodSegment *segment = (RodSegment*)node->object;
			noZeroDinv = segment->isDynamic();
		}
		if (noZeroDinv) // if DInv == 0 child value is 0 and node->soln is not altered
		{
			node->soln = node->DLDLT.solve(node->soln);

			if (node->parent != NULL)
			{
				node->soln -= node->J * node->parent->soln;
			}
		}
		else
		{
			node->soln.setZero(); // segment of node is not dynamic
		}

		if (node->isconstraint)
		{
			lambdaSums[node->index] += node->soln;
		}
	}

	// compute position and orientation updates
	for (nodeIter = forward[intervalIndex].begin(); nodeIter != forward[intervalIndex].end(); nodeIter++)
	{
		Node *node = *nodeIter;
		if (!node->isconstraint)
		{
			RodSegment *segment = (RodSegment *)node->object;
			if (!segment->isDynamic())
			{
				break;
			}

			const Vector6r & soln(node->soln);
			Vector3r deltaXSoln = Vector3r(-soln[0], -soln[1], -soln[2]);
			corr_x[node->index] = deltaXSoln;

			Matrix<Real, 4, 3> G;
			computeMatrixG(segment->Rotation(), G);
			Quaternionr deltaQSoln;
			deltaQSoln.coeffs() = G * Vector3r(-soln[3], -soln[4], -soln[5]);
			corr_q[node->index] = deltaQSoln;
		}
	}
	return true;
}

bool pbdstiffrods::init_DirectPositionBasedSolverForStiffRodsConstraint(
	std::vector<RodConstraint*> &rodConstraints, 
	std::vector<RodSegment*> & rodSegments, 
	Interval* &intervals, 
	int &numberOfIntervals, 
	std::list <Node*> * &forward, 
	std::list <Node*> * &backward, 
	Node* &root, 
	const std::vector<Vector3r> &constraintPositions,
	const std::vector<Real> &averageRadii,
	const std::vector<Real> &youngsModuli,
	const std::vector<Real> &torsionModuli, 
	std::vector<Vector6r> & RHS, 
	std::vector<Vector6r> & lambdaSums, 
	std::vector<std::vector<Matrix3r>> & bendingAndTorsionJacobians, 
	std::vector<Vector3r> & corr_x, 
	std::vector<Quaternionr> & corr_q
	)
{
	// init constraints
	for (size_t cIdx(0); cIdx < rodConstraints.size(); ++cIdx)
	{
		RodConstraint * constraint(rodConstraints[cIdx]);
		RodSegment * segment0(rodSegments[constraint->segmentIndex(0)]);
		RodSegment * segment1(rodSegments[constraint->segmentIndex(1)]);

		init_StretchBendingTwistingConstraint(
			segment0->Position(), segment0->Rotation(), segment1->Position(), segment1->Rotation(),
			constraintPositions[cIdx], averageRadii[cIdx], constraint->getAverageSegmentLength(),
			youngsModuli[cIdx], torsionModuli[cIdx], constraint->getConstraintInfo(),
			constraint->getStiffnessCoefficientK(), constraint->getRestDarbouxVector());
	}
	
	// compute tree data structure for direct solver
	initTree(rodConstraints, rodSegments, intervals, numberOfIntervals, forward, backward, root);

	RHS.resize(rodConstraints.size());
	std::fill(RHS.begin(), RHS.end(), Vector6r::Zero());

	lambdaSums.resize(rodConstraints.size());
	std::fill(lambdaSums.begin(), lambdaSums.end(), Vector6r::Zero());

	bendingAndTorsionJacobians.resize(rodConstraints.size());
	std::vector<Matrix3r> sampleJacobians(2);
	sampleJacobians[0].setZero();
	sampleJacobians[1].setZero();
	std::fill(bendingAndTorsionJacobians.begin(), bendingAndTorsionJacobians.end(), sampleJacobians);

	corr_x.resize(rodSegments.size());
	std::fill(corr_x.begin(), corr_x.end(), Vector3r::Zero());

	corr_q.resize(rodSegments.size());
	std::fill(corr_q.begin(), corr_q.end(), Quaternionr());

	return true;
}

bool pbdstiffrods::initBeforeProjection_DirectPositionBasedSolverForStiffRodsConstraint(
	const std::vector<RodConstraint*> &rodConstraints,
	const Real inverseTimeStepSize,
	std::vector<Vector6r> & lambdaSums
	)
{
	for (size_t cIdx(0); cIdx < rodConstraints.size(); ++cIdx)
	{
		RodConstraint * constraint(rodConstraints[cIdx]);

		initBeforeProjection_StretchBendingTwistingConstraint(
			constraint->getStiffnessCoefficientK(),
			inverseTimeStepSize,
			constraint->getAverageSegmentLength(),
			constraint->getStretchCompliance(),
			constraint->getBendingAndTorsionCompliance(),
			lambdaSums[cIdx]);
	}
	return true;
}



bool pbdstiffrods::update_DirectPositionBasedSolverForStiffRodsConstraint(
	const std::vector<RodConstraint*> &rodConstraints,
	const std::vector<RodSegment*> & rodSegments
	)
{
	// update rod constraints
	for (size_t cIdx(0); cIdx < rodConstraints.size(); ++cIdx)
	{
		RodConstraint * constraint(rodConstraints[cIdx]);
		RodSegment * segment0(rodSegments[constraint->segmentIndex(0)]);
		RodSegment * segment1(rodSegments[constraint->segmentIndex(1)]);

		update_StretchBendingTwistingConstraint(
			segment0->Position(), segment0->Rotation(), segment1->Position(), segment1->Rotation(),
			constraint->getConstraintInfo());
	}
	return true;
}



bool pbdstiffrods::solve_DirectPositionBasedSolverForStiffRodsConstraint(
	const std::vector<RodConstraint*> &rodConstraints, 
	std::vector<RodSegment*> & rodSegments, 
	const Interval* intervals, 
	const int &numberOfIntervals, 
	std::list <Node*> * forward, 
	std::list <Node*> * backward, 
	std::vector<Vector6r> & RHS, 
	std::vector<Vector6r> & lambdaSums, 
	std::vector<std::vector<Matrix3r>> & bendingAndTorsionJacobians, 
	std::vector<Vector3r> & corr_x, 
	std::vector<Quaternionr> & corr_q
	)
{
	for (int i = 0; i < numberOfIntervals; i++)
	{
		factor(i, rodConstraints, rodSegments, intervals,
			forward, backward, RHS, lambdaSums, bendingAndTorsionJacobians);
	}
	for (int i = 0; i < numberOfIntervals; i++)
	{
		solve(i, forward, backward, RHS, lambdaSums, corr_x, corr_q);
	}
	return true;
}

bool pbdstiffrods::init_StretchBendingTwistingConstraint(
	const Vector3r &x0, const Quaternionr &q0, 
	const Vector3r &x1, const Quaternionr &q1, 
	const Vector3r &constraintPosition, 
	const Real averageRadius, 
	const Real averageSegmentLength, 
	const Real youngsModulus, 
	const Real torsionModulus, 
	Matrix<Real, 3, 4, DontAlign> &constraintInfo, 
	Vector3r &stiffnessCoefficientK, 
	Vector3r &restDarbouxVector
	)
{
	// constraintInfo contains
	// 0:	connector in segment 0 (local)
	// 1:	connector in segment 1 (local)
	// 2:	connector in segment 0 (global)
	// 3:	connector in segment 1 (global)

	// transform in local coordinates
	const Matrix3r rot0T = q0.matrix().transpose();
	const Matrix3r rot1T = q1.matrix().transpose();

	constraintInfo.col(0) = rot0T * (constraintPosition - x0);
	constraintInfo.col(1) = rot1T * (constraintPosition - x1);
	constraintInfo.col(2) = constraintPosition;
	constraintInfo.col(3) = constraintPosition;

	// compute bending and torsion stiffness of the K matrix diagonal; assumption: the rod axis follows the y-axis of the local frame as with Blender's armatures
	Real secondMomentOfArea((static_cast<Real>(3.14159265358979323846) * static_cast<Real>(0.25)) * std::pow(averageRadius, static_cast<Real>(4.0)));
	Real bendingStiffness(youngsModulus * secondMomentOfArea);
	Real torsionStiffness(static_cast<Real>(2.0) * torsionModulus * secondMomentOfArea);
	stiffnessCoefficientK = Vector3r(bendingStiffness, torsionStiffness, bendingStiffness);

	// compute rest Darboux vector
	computeDarbouxVector(q0, q1, averageSegmentLength, restDarbouxVector);

	return true;
}

bool pbdstiffrods::initBeforeProjection_StretchBendingTwistingConstraint(
	const Vector3r &stiffnessCoefficientK,
	const Real inverseTimeStepSize,
	const Real averageSegmentLength,
	Vector3r &stretchCompliance,
	Vector3r &bendingAndTorsionCompliance,
	Vector6r &lambdaSum
	)
{

	Real inverseTSQuadratic(inverseTimeStepSize*inverseTimeStepSize);

	// compute compliance parameter of the stretch constraint part
	const Real stretchRegularizationParameter(static_cast<Real>(1.E-10));
	stretchCompliance = Vector3r(
		stretchRegularizationParameter * inverseTSQuadratic,
		stretchRegularizationParameter * inverseTSQuadratic,
		stretchRegularizationParameter * inverseTSQuadratic);

	// compute compliance parameter of the bending and torsion constraint part
	bendingAndTorsionCompliance = Vector3r(
		inverseTSQuadratic / stiffnessCoefficientK(0),
		inverseTSQuadratic / stiffnessCoefficientK(1),
		inverseTSQuadratic / stiffnessCoefficientK(2));
	bendingAndTorsionCompliance *= static_cast<Real>(1.0) / averageSegmentLength;

	// set sum of delta lambda values to zero
	lambdaSum.setZero();
	return true;
}

bool pbdstiffrods::update_StretchBendingTwistingConstraint(
	const Vector3r &x0, const Quaternionr &q0, 
	const Vector3r &x1, const Quaternionr &q1, 
	Matrix<Real, 3, 4, DontAlign> &constraintInfo
	)
{
	// constraintInfo contains
	// 0:	connector in segment 0 (local)
	// 1:	connector in segment 1 (local)
	// 2:	connector in segment 0 (global)
	// 3:	connector in segment 1 (global)

	// compute world space positions of connectors
	const Matrix3r rot0 = q0.matrix();
	const Matrix3r rot1 = q1.matrix();
	constraintInfo.col(2) = rot0 * constraintInfo.col(0) + x0;
	constraintInfo.col(3) = rot1 * constraintInfo.col(1) + x1;

	return true;
}

bool pbdstiffrods::solve_StretchBendingTwistingConstraint(
	const Real invMass0,
	const Vector3r &x0,
	const Matrix3r &inertiaInverseW0,
	const Quaternionr &q0,
	const Real invMass1,
	const Vector3r &x1,
	const Matrix3r &inertiaInverseW1,
	const Quaternionr &q1,
	const Vector3r &restDarbouxVector, 
	const Real averageSegmentLength, 
	const Vector3r &stretchCompliance, 
	const Vector3r &bendingAndTorsionCompliance, 
	const Matrix<Real, 3, 4, DontAlign> &constraintInfo,
	Vector3r &corr_x0, Quaternionr &corr_q0, 
	Vector3r &corr_x1, Quaternionr &corr_q1, 
	Vector6r &lambdaSum
	)
{
	// compute Darboux vector (Equation (7))
	Vector3r omega;
	computeDarbouxVector(q0, q1, averageSegmentLength, omega);

	// compute bending and torsion Jacobians (Equation (10) and Equation (11))
	Matrix<Real, 3, 4> jOmega0, jOmega1;
	computeBendingAndTorsionJacobians(q0, q1, averageSegmentLength, jOmega0, jOmega1);

	// compute G matrices (Equation (27))
	Matrix<Real, 4, 3> G0, G1;
	computeMatrixG(q0, G0);
	computeMatrixG(q1, G1);

	Matrix3r jOmegaG0(jOmega0*G0),
		jOmegaG1(jOmega1*G1);

	// Compute zero-stretch part of constraint violation (Equation (23))
	const Vector3r &connector0 = constraintInfo.col(2);
	const Vector3r &connector1 = constraintInfo.col(3);
	Vector3r stretchViolation = connector0 - connector1;

	// Compute bending and torsion part of constraint violation  (Equation (23))
	Vector3r bendingAndTorsionViolation = omega - restDarbouxVector;

	// fill right hand side of the linear equation system (Equation (19))
	Vector6r rhs;
	rhs.block<3, 1>(0, 0) = -stretchViolation - Vector3r(
		stretchCompliance[0] * lambdaSum[0],
		stretchCompliance[1] * lambdaSum[1],
		stretchCompliance[2] * lambdaSum[2]);

	rhs.block<3, 1>(3, 0) = -bendingAndTorsionViolation - Vector3r(
		bendingAndTorsionCompliance[0] * lambdaSum[3],
		bendingAndTorsionCompliance[1] * lambdaSum[4],
		bendingAndTorsionCompliance[2] * lambdaSum[5]);

	// compute matrix of the linear equation system (using Equations (25), (26), and (28) in Equation (19))
	Matrix6r JMJT(Matrix6r::Zero());

	// compute stretch block
	Matrix3r K1, K2;
	computeMatrixK(connector0, invMass0, x0, inertiaInverseW0, K1);
	computeMatrixK(connector1, invMass1, x1, inertiaInverseW1, K2);
	JMJT.block<3, 3>(0, 0) = K1 + K2;

	// compute coupling blocks
	const Vector3r ra = connector0 - x0;
	const Vector3r rb = connector1 - x1;

	Matrix3r ra_crossT, rb_crossT;
	matrix::crossProductMatrix(-ra, ra_crossT); // use -ra to get the transpose
	matrix::crossProductMatrix(-rb, rb_crossT); // use -rb to get the transpose

	Matrix3r offdiag(Matrix3r::Zero());
	if (invMass0 != 0.0)
	{
		offdiag = jOmegaG0 * inertiaInverseW0 * ra_crossT * (-1);
	}

	if (invMass1 != 0.0)
	{
		offdiag += jOmegaG1 * inertiaInverseW1 * rb_crossT;
	}
	JMJT.block<3, 3>(3, 0) = offdiag;
	JMJT.block<3, 3>(0, 3) = offdiag.transpose();

	// compute bending and torsion block
	Matrix3r MInvJT0(inertiaInverseW0 * jOmegaG0.transpose());
	Matrix3r MInvJT1(inertiaInverseW1 * jOmegaG1.transpose());

	Matrix3r JMJTOmega(Matrix3r::Zero());
	if (invMass0 != 0.0)
	{
		JMJTOmega = jOmegaG0*MInvJT0;
	}

	if (invMass1 != 0.0)
	{
		JMJTOmega += jOmegaG1*MInvJT1;
	}
	JMJT.block<3, 3>(3, 3) = JMJTOmega;

	// add compliance
	JMJT(0, 0) += stretchCompliance(0);
	JMJT(1, 1) += stretchCompliance(1);
	JMJT(2, 2) += stretchCompliance(2);
	JMJT(3, 3) += bendingAndTorsionCompliance(0);
	JMJT(4, 4) += bendingAndTorsionCompliance(1);
	JMJT(5, 5) += bendingAndTorsionCompliance(2);
	
	// solve linear equation system (Equation 19)
	auto decomposition(JMJT.ldlt());
	Vector6r deltaLambda(decomposition.solve(rhs));

	// update sum of delta lambda values for next Gauss-Seidel solver iteration step
	lambdaSum += deltaLambda;

	// compute position and orientation updates (using Equations (25), (26), and (28) in Equation (20))
	Vector3r deltaLambdaStretch(deltaLambda.block<3, 1>(0, 0)),
		deltaLambdaBendingAndTorsion(deltaLambda.block<3, 1>(3, 0));
	corr_x0.setZero();
	corr_x1.setZero();
	corr_q0.coeffs().setZero();
	corr_q1.coeffs().setZero();

	if (invMass0 != 0.)
	{
		corr_x0 += invMass0 * deltaLambdaStretch;
		corr_q0.coeffs() += G0 * (inertiaInverseW0 * ra_crossT * (-1 * deltaLambdaStretch) +
			MInvJT0 * deltaLambdaBendingAndTorsion);
	}

	if (invMass1 != 0.)
	{
		corr_x1 -= invMass1 * deltaLambdaStretch;
		corr_q1.coeffs() += G1 * (inertiaInverseW1 * rb_crossT * deltaLambdaStretch +
			MInvJT1 * deltaLambdaBendingAndTorsion);
	}

	return true;
}