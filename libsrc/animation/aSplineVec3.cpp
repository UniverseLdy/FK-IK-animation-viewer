#include "aSplineVec3.h"
#include <algorithm>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;
#pragma warning(disable:4018)
#pragma warning(disable:4244)

ASplineVec3::ASplineVec3() : mInterpolator(new ALinearInterpolatorVec3())
{
}

ASplineVec3::~ASplineVec3()
{
	delete mInterpolator;
}

void ASplineVec3::setFramerate(double fps)
{
	mInterpolator->setFramerate(fps);
}

double ASplineVec3::getFramerate() const
{
	return mInterpolator->getFramerate();
}

void ASplineVec3::setLooping(bool loop)
{
	mLooping = loop;
}

bool ASplineVec3::getLooping() const
{
	return mLooping;
}

void ASplineVec3::setInterpolationType(ASplineVec3::InterpolationType type)
{
	double fps = getFramerate();

	delete mInterpolator;
	switch (type)
	{
	case LINEAR: mInterpolator = new ALinearInterpolatorVec3();
		break;
	case CUBIC_BERNSTEIN: mInterpolator = new ABernsteinInterpolatorVec3();
		break;
	case CUBIC_CASTELJAU: mInterpolator = new ACasteljauInterpolatorVec3();
		break;
	case CUBIC_MATRIX: mInterpolator = new AMatrixInterpolatorVec3();
		break;
	case CUBIC_HERMITE: mInterpolator = new AHermiteInterpolatorVec3();
		break;
	case CUBIC_BSPLINE: mInterpolator = new ABSplineInterpolatorVec3();
		break;
	};

	mInterpolator->setFramerate(fps);
	computeControlPoints();
	cacheCurve();
}

ASplineVec3::InterpolationType ASplineVec3::getInterpolationType() const
{
	return mInterpolator->getType();
}

void ASplineVec3::editKey(int keyID, const vec3& value)
{
	assert(keyID >= 0 && keyID < mKeys.size());
	mKeys[keyID].second = value;
	computeControlPoints();
	cacheCurve();
}

void ASplineVec3::editControlPoint(int ID, const vec3& value)
{
	assert(ID >= 0 && ID < mCtrlPoints.size() + 2);
	if (ID == 0)
	{
		mStartPoint = value;
		computeControlPoints();
	}
	else if (ID == mCtrlPoints.size() + 1)
	{
		mEndPoint = value;
		computeControlPoints();
	}
	else mCtrlPoints[ID - 1] = value;
	cacheCurve();
}

void ASplineVec3::appendKey(double time, const vec3& value, bool updateCurve)
{
	mKeys.push_back(Key(time, value));

	if (mKeys.size() >= 2)
	{
		int totalPoints = mKeys.size();

		//If there are more than 1 interpolation point, set up the 2 end points to help determine the curve.
		//They lie on the tangent of the first and last interpolation points.
		vec3 tmp = mKeys[0].second - mKeys[1].second;
		double n = tmp.Length();
		mStartPoint = mKeys[0].second + (tmp / n) * n * 0.25;
		// distance to endpoint is 25% of distance between first 2 points

		tmp = mKeys[totalPoints - 1].second - mKeys[totalPoints - 2].second;
		n = tmp.Length();
		mEndPoint = mKeys[totalPoints - 1].second + (tmp / n) * n * 0.25;
	}

	if (updateCurve)
	{
		computeControlPoints();
		cacheCurve();
	}
}

void ASplineVec3::appendKey(const vec3& value, bool updateCurve)
{
	if (mKeys.size() == 0)
	{
		appendKey(0, value, updateCurve);
	}
	else
	{
		double lastT = mKeys[mKeys.size() - 1].first;
		appendKey(lastT + 1, value, updateCurve);
	}
}

void ASplineVec3::deleteKey(int keyID)
{
	assert(keyID >= 0 && keyID < mKeys.size());
	mKeys.erase(mKeys.begin() + keyID);
	computeControlPoints();
	cacheCurve();
}

vec3 ASplineVec3::getKey(int keyID)
{
	assert(keyID >= 0 && keyID < mKeys.size());
	return mKeys[keyID].second;
}

int ASplineVec3::getNumKeys() const
{
	return mKeys.size();
}

vec3 ASplineVec3::getControlPoint(int ID)
{
	assert(ID >= 0 && ID < mCtrlPoints.size() + 2);
	if (ID == 0) return mStartPoint;
	else if (ID == mCtrlPoints.size() + 1) return mEndPoint;
	else return mCtrlPoints[ID - 1];
}

int ASplineVec3::getNumControlPoints() const
{
	return mCtrlPoints.size() + 2; // include endpoints
}

void ASplineVec3::clear()
{
	mKeys.clear();
}

double ASplineVec3::getDuration() const
{
	return mKeys[mKeys.size() - 1].first;
}

double ASplineVec3::getNormalizedTime(double t) const
{
	return (t / getDuration());
}

vec3 ASplineVec3::getValue(double t)
{
	if (mCachedCurve.size() == 0) return vec3();

	double dt = mInterpolator->getDeltaTime();
	int rawi = (int)(t / dt); // assumes uniform spacing
	int i = rawi % mCachedCurve.size();
	double frac = t - rawi * dt;
	int inext = i + 1;
	if (!mLooping) inext = std::min<int>(inext, mCachedCurve.size() - 1);
	else inext = inext % mCachedCurve.size();

	vec3 v1 = mCachedCurve[i];
	vec3 v2 = mCachedCurve[inext];
	vec3 v = v1 * (1 - frac) + v2 * frac;
	return v;
}

void ASplineVec3::cacheCurve()
{
	mInterpolator->interpolate(mKeys, mCtrlPoints, mCachedCurve);
}

void ASplineVec3::computeControlPoints()
{
	mInterpolator->computeControlPoints(mKeys, mCtrlPoints, mStartPoint, mEndPoint);
}

int ASplineVec3::getNumCurveSegments() const
{
	return mCachedCurve.size();
}

vec3 ASplineVec3::getCurvePoint(int i) const
{
	return mCachedCurve[i];
}

//---------------------------------------------------------------------
AInterpolatorVec3::AInterpolatorVec3(ASplineVec3::InterpolationType t) : mDt(1.0 / 120.0), mType(t)
{
}

void AInterpolatorVec3::setFramerate(double fps)
{
	mDt = 1.0 / fps;
}

double AInterpolatorVec3::getFramerate() const
{
	return 1.0 / mDt;
}

double AInterpolatorVec3::getDeltaTime() const
{
	return mDt;
}

void AInterpolatorVec3::interpolate(const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints, std::vector<vec3>& curve)
{
	vec3 val = 0.0;
	double u = 0.0;

	curve.clear();

	int numSegments = keys.size() - 1;
	for (int segment = 0; segment < numSegments; segment++)
	{
		for (double t = keys[segment].first; t < keys[segment + 1].first - FLT_EPSILON; t += mDt)
		{
			// TODO: Compute u, fraction of duration between segment and segmentnext, for example,
			// u = 0.0 when t = keys[segment-1].first  
			// u = 1.0 when t = keys[segment].first
			u = ((t - keys[segment].first) / (keys[segment + 1].first - keys[segment].first));
			val = interpolateSegment(keys, ctrlPoints, segment, u);
			curve.push_back(val);
		}
	}
	/*add last point*/
	if (keys.size() > 1)
	{
		u = 1.0;
		val = interpolateSegment(keys, ctrlPoints, numSegments - 1, u);
		curve.push_back(val);
	}
}


vec3 ALinearInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 curveValue(0, 0, 0);
	vec3 key0 = keys[segment].second;
	vec3 key1 = keys[segment + 1].second;

	// TODO: 
	//Step 1: Create a Lerp helper function
	//Step 2: Linear i`nterpolate between key0 and key1 so that u = 0 returns key0 and u = 1 returns key1
	curveValue = key1 * u + key0 * (1 - u);
	return curveValue;
}

vec3 ABernsteinInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);
	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];
	// Step2: Compute the interpolated value f(u) point using  Bernstein polynomials
	curveValue = b0 * (1 - u)*(1 - u)*(1 - u) + b1 * 3 * u *(1 - u) *(1 - u) + b2 * 3 * u *u * (1 - u) + b3 * u*u*u;
	return curveValue;
}


vec3 ACasteljauInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  deCsteljau alogithm
	vec3 b01, b11, b21, b02, b12, b03;
	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];
	b01 = b0 * (1 - u) + b1 * u;
	b11 = b1 * (1 - u) + b2 * u;
	b21 = b2 * (1 - u) + b3 * u;
	b02 = b01 * (1 - u) + b11 * u;
	b12 = b11 * (1 - u) + b21 * u;
	b03 = b02 * (1 - u) + b12 * u;
	curveValue = b03;
	return curveValue;
}

vec3 AMatrixInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 b0;
	vec3 b1;
	vec3 b2;
	vec3 b3;
	vec3 curveValue(0, 0, 0);

	// TODO: 
	// Step1: Get the 4 control points, b0, b1, b2 and b3 from the ctrlPoints vector
	// Step2: Compute the interpolated value f(u) point using  matrix method f(u) = GMU
	b0 = ctrlPoints[4 * segment];
	b1 = ctrlPoints[4 * segment + 1];
	b2 = ctrlPoints[4 * segment + 2];
	b3 = ctrlPoints[4 * segment + 3];
	Eigen::Matrix4d G;
	G << b0[0], b1[0], b2[0], b3[0],
		b0[1], b1[1], b2[1], b3[1],
		b0[2], b1[2], b2[2], b3[2],
		0, 0, 0, 0;
	Eigen::Matrix4d M;
	M << 1, -3, 3, -1,
		0, 3, -6, 3,
		0, 0, 3, -3,
		0, 0, 0, 1;
	Eigen::Matrix4d U;
	U << 1, 0, 0, 0,
		u, 0, 0, 0,
		u*u, 0, 0, 0,
		u*u*u, 0, 0, 0;
	Eigen::Matrix4d F;
	F = G * M *U;
	curveValue[0] = F(0, 0);
	curveValue[1] = F(1, 0);
	curveValue[2] = F(2, 0);

	// Hint: Using Eigen::MatrixXd data representations for a matrix operations

	return curveValue;
}

vec3 AHermiteInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 p0 = keys[segment].second;
	vec3 p1 = keys[segment + 1].second;
	vec3 q0 = ctrlPoints[segment]; // slope at p0
	vec3 q1 = ctrlPoints[segment + 1]; // slope at p1
	vec3 curveValue(0, 0, 0);
	curveValue = (2 * u*u*u - 3 * u*u + 1)*p0 + (-2 * u*u*u + 3 * u*u)*p1 + (u*u*u - 2 * u*u + u)*q0 + (u*u*u - u * u)*q1;
	// TODO: Compute the interpolated value h(u) using a cubic Hermite polynomial  

	return curveValue;
}
vector<double> ComputeN(
	vector<double> knots,
	int n,
	int j,
	double t)
{
	vector<double> N;
	N.clear();
	Eigen::MatrixXd A(n + 1, n + 2);
	A.setZero(n + 1, n + 2);
	A(0, n) = 1;
	for (int a = 1; a <= n; a++)
	{
		for (int b = 0; b <= n; b++)
		{
			A(a, b) = A(a - 1, b)*(t - knots[j - 3 + b]) / (knots[j - 3 + a + b] - knots[j - 3 + b]) + A(a - 1, b + 1)*(knots[j - 3 + a + b + 1] - t) / (knots[j - 3 + a + b + 1] - knots[j - 3 + b + 1]);
		}
	}
	/*cout << A << endl;*/
	for (int i = 0; i <= n; i++)
	{
		N.push_back(A(n, i));
	}
	/*cout << A << endl;*/
	//cout << N[2] << endl;
	return N;
}
vec3 ABSplineInterpolatorVec3::interpolateSegment(
	const std::vector<ASplineVec3::Key>& keys,
	const std::vector<vec3>& ctrlPoints,
	int segment, double u)
{
	vec3 c0, c1, c2, c3;
	vec3 curveValue(0, 0, 0);
	c0 = ctrlPoints[segment];
	c1 = ctrlPoints[segment + 1];
	c2 = ctrlPoints[segment + 2];
	c3 = ctrlPoints[segment + 3];
	Eigen::Matrix4d G;
	G << c0[0], c1[0], c2[0], c3[0],
		c0[1], c1[1], c2[1], c3[1],
		c0[2], c1[2], c2[2], c3[2],
		0, 0, 0, 0;
	Eigen::Matrix4d M;
	M << 0.166667, -0.5, 0.5, -0.166667,
		0.666667, 0, -1, 0.5,
		0.166667, 0.5, 0.5, -0.5,
		0, 0, 0, 0.166667;
	/* cout << M << endl;*/
	Eigen::Matrix4d U;
	U << 1, 0, 0, 0,
		u, 0, 0, 0,
		u*u, 0, 0, 0,
		u*u*u, 0, 0, 0;
	Eigen::Matrix4d F;
	F = G * M *U;
	curveValue[0] = F(0, 0);
	curveValue[1] = F(1, 0);
	curveValue[2] = F(2, 0);

	// Hint: Create a recursive helper function N(knots,n,j,t) to calculate BSpline basis function values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n =3 for cubic)
	//     j = curve interval on knot vector in which to interpolate
	//     t = time value	

	// Step 1: determine the index j
	// Step 2: compute the n nonzero Bspline Basis functions N given j
	// Step 3: get the corresponding control points from the ctrlPoints vector
	// Step 4: compute the Bspline curveValue at time t


	return curveValue;
}

void ACubicInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;
	vec3 b0, b1, b2, b3;
	int num = keys.size();
	startPoint = keys[0].second;
	endPoint = keys[num - 1].second;
	if (keys.size() == 2)
	{
		b0 = startPoint;
		b3 = endPoint;
		b1 = startPoint + ((endPoint - startPoint) / 3);
		b2 = endPoint - ((endPoint - startPoint) / 3);
		ctrlPoints.push_back(b0);
		ctrlPoints.push_back(b1);
		ctrlPoints.push_back(b2);
		ctrlPoints.push_back(b3);
	}
	if (keys.size() >= 3)
	{
		ctrlPoints.clear();
		b0 = startPoint;
		b1 = startPoint + ((keys[1].second - startPoint) / 3);
		ctrlPoints.push_back(b0);
		ctrlPoints.push_back(b1);
		for (int i = 1; i < keys.size() - 1; i++)
		{
			b3 = keys[i].second;
			b2 = keys[i].second - ((keys[i + 1].second - keys[i - 1].second) / 6);
			ctrlPoints.push_back(b2);
			ctrlPoints.push_back(b3);
			b0 = keys[i].second;
			b1 = keys[i].second + ((keys[i + 1].second - keys[i - 1].second) / 6);
			// TODO: compute b0, b1, b2, b3
			ctrlPoints.push_back(b0);
			ctrlPoints.push_back(b1);
		}
		b3 = endPoint;
		b2 = endPoint - ((endPoint - keys[num - 1].second) / 3);
		ctrlPoints.push_back(b2);
		ctrlPoints.push_back(b3);
	}
}

void AHermiteInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPoint, vec3& endPoint)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;

	int numKeys = keys.size();

	Eigen::MatrixXd A(numKeys, numKeys);
	//A = MatrixXd::Zero(numKeys, numKeys);
	A.setZero(numKeys, numKeys);
	//for (int a = 0; a < numKeys;a++)
	//{
	   // for (int b = 0; b < numKeys; b++)
	   // {
		  //  A(a, b) = 0;
	   // }
	//}

	A(0, 0) = 2;
	A(0, 1) = 1;
	A(numKeys - 1, numKeys - 2) = 1;
	A(numKeys - 1, numKeys - 1) = 2;
	Eigen::MatrixXd D(3, numKeys);
	startPoint = keys[0].second;
	endPoint = keys[numKeys - 1].second;
	vec3 Startslope, Endslope;
	Startslope = 3 * (keys[1].second - startPoint);
	Endslope = 3 * (endPoint - keys[numKeys - 2].second);
	D(0, 0) = Startslope[0]; D(1, 0) = Startslope[1];  D(2, 0) = Startslope[2];
	D(0, numKeys - 1) = Endslope[0]; D(1, numKeys - 1) = Endslope[1]; D(2, numKeys - 1) = Endslope[2];
	for (int i = 1; i < numKeys - 1; i++)
	{
		A(i, i - 1) = 1;
		A(i, i) = 4;
		A(i, i + 1) = 1;
		vec3 slope = 3 * (keys[i + 1].second - keys[i - 1].second);
		D(0, i) = slope[0];
		D(1, i) = slope[1];
		D(2, i) = slope[2];
	}
	Eigen::MatrixXd C;
	C = D * (A.inverse());
	for (int j = 0; j < numKeys; j++)
	{
		vec3 Cp;
		Cp[0] = C(0, j);
		Cp[1] = C(1, j);
		Cp[2] = C(2, j);
		ctrlPoints.push_back(Cp);
	}
	// TODO: 
	// For each key point pi, compute the corresonding value of the slope pi_prime.
	// Hints: Using Eigen::MatrixXd for a matrix data structures, 
	// this can be accomplished by solving the system of equations AC=D for C.
	// Don't forget to save the values computed for C in ctrlPoints
	// For clamped endpoint conditions, set 1st derivative at first and last points (p0 and pm) to s0 and s1, respectively
	// For natural endpoints, set 2nd derivative at first and last points (p0 and pm) equal to 0

	// Step 1: Initialize A
	// Step 2: Initialize D
	// Step 3: Solve AC=D for C
	// Step 4: Save control points in ctrlPoints
}

void ComputedN(std::vector<double>& dN,
	vector<double> knots,
	int n,
	double t,
	int j,
	int l)
{
	Eigen::MatrixXd A(n + 1, n + 2);
	A.setZero(n + 1, n + 2);
	A(0, n) = 1;
	for (int a = 1; a <= (n - l); a++)
	{
		for (int b = 0; b <= n; b++)
		{
			A(a, b) = A(a - 1, b)*(t - knots[j - 3 + b]) / (knots[j - 3 + a + b] - knots[j - 3 + b]) + A(a - 1, b + 1)*(knots[j - 3 + a + b + 1] - t) / (knots[j - 3 + a + b + 1] - knots[j - 3 + b + 1]);
		}
	}
	for (int a = (n - l + 1); a <= n; a++)
		for (int b = 0; b <= n; b++)
		{
			A(a, b) = a * (A(a - 1, b) / (knots[j - 3 + b + a] - knots[j - 3 + b]) - A(a - 1, b + 1) / (knots[j - 3 + b + a + 1] - knots[j - 3 + b + 1]));
		}
	for (int i = 0; i <= n; i++)
	{
		dN.push_back(A(n, i));
	}
}
void ABSplineInterpolatorVec3::computeControlPoints(
	const std::vector<ASplineVec3::Key>& keys,
	std::vector<vec3>& ctrlPoints,
	vec3& startPt, vec3& endPt)
{
	ctrlPoints.clear();
	if (keys.size() <= 1) return;
	if (keys.size() >= 2)
	{
		int l = 2;
		int numkeys = keys.size();
		int m = numkeys - 1;
		startPt = keys[0].second;
		endPt = keys[numkeys - 1].second;
		double t0, tm;
		t0 = 0;
		tm = m;
		int n = 3;
		vector<double> knots;
		vector<double> dN;
		double ptmu;
		ptmu = 1;
		double ptm = 0;
		double t = t0;
		for (int i = n; i > 0; i--)
		{
			ptm = t0 - (i * ptmu);
			knots.push_back(ptm);
		}
		for (int i = 0; i <= m; i++)
		{
			knots.push_back(t);
			t += 1;
		}
		for (int i = 1; i <= n; i++)
		{
			ptm = tm + (i * ptmu);
			knots.push_back(ptm);
		}
		ComputedN(dN, knots, n, t0, 3, l);
		ComputedN(dN, knots, n, tm, m + 2, l);
		Eigen::MatrixXd A(m + 3, m + 3);
		A.setZero();
		for (int i = 0; i <= n; i++)
		{
			A(0, i) = dN[i];
		}
		for (int i = 0; i <= n; i++)
		{
			A(numkeys + 1, numkeys + 1 - n + i) = dN[i + n + 1];
		}
		vector<double> N;
		for (int i = 1; i < numkeys; i++)
		{
			N = ComputeN(knots, n, i + 2, knots[i + 2]);
			for (int a = i - 1; a <= i - 1 + n; a++)
			{
				A(i, a) = N[a - i + 1];
			}
			N.clear();
		}
		N = ComputeN(knots, n, m + 2, tm);
		for (int i = 0; i <= n; i++)
		{
			A(numkeys, numkeys - 2 + i) = N[i];
		}
		N.clear();
		/*cout << A << endl;*/
		Eigen::MatrixXd D(numkeys + 2, 3);
		D.setZero();

		for (int i = 1; i <= numkeys; i++)
		{
			for (int j = 0; j <= 2; j++)
			{
				D(i, j) = keys[i - 1].second[j];
			}
		}
		/*cout << D << endl;*/
		Eigen::MatrixXd C;
		C = (A.inverse())*D;
		/*cout << C << endl;*/
		for (int j = 0; j < numkeys + 2; j++)
		{
			vec3 Cp;
			Cp[0] = C(j, 0);
			Cp[1] = C(j, 1);
			Cp[2] = C(j, 2);
			ctrlPoints.push_back(Cp);
		}
	}
	// TODO: c
	// Hints: 

	// 1. use Eigen::MatrixXd to calculate the control points by solving the system of equations AC=D for C

	// 2. Create a recursive helper function dN(knots,n,t,l) to calculate derivative BSpline values at t, where
	//     knots = knot array
	//	   n = degree of the spline curves (n =3 for cubic)
	//     j = interval on knot vector in which to interpolate
	//     t = time value
	//     l = derivative (l = 1 => 1st derivative)

	// Step 1: Calculate knot vector using a uniform BSpline
	//         (assume knots are evenly spaced 1 apart and the start knot is at time = 0.0)

	// Step 2: Calculate A matrix  for a natural BSpline
	//         (Set 2nd derivative at t0 and tm to zero, where tm is the last point knot; m = #segments)

	// Step 3: Calculate  D matrix composed of our target points to interpolate

	// Step 4: Solve AC=D for C 

	// Step 5: save control points in ctrlPoints
}
