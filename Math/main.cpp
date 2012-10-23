/**************************************
* Math Test Program
* By Will Perone
\**************************************/

#include "vector3.h"
#include "matrix2.h"
#include "matrix3.h"
#include "matrix4.h"
#include "spline.h"
#include "line3.h"
#include "intersect.h"
#include "utility.h"
#include "rect.h"
#include "quaternion.h"
#include "random.h"
#include <iostream>
#include <iomanip>
using namespace std;


// outputs a 3x3 matrix
void outputmat3(matrix3 &m)
{
	int i, j;
	for (i= 0; i < 3; i++)
	{
	for (j= 0; j < 3; j++)
		cout << setw(8) << setprecision(4) << m[i][j] << " ";
	cout << endl;
	}
}

// outputs a 4x4 matrix
void outputmat4(matrix4 &m)
{
	int i, j;
	for (i= 0; i < 4; i++)
	{		
		for (j= 0; j < 4; j++)
		cout << setw(8) << setprecision(4) << m[i][j] << " ";
	cout << endl;
	}
}


#define outputv2(v) cout << #v ".x= " << v.x << " " #v ".y= " << v.y << endl;
#define outputv3(v) cout << #v ".x= " << v.x << " " #v ".y= " << v.y << " " #v ".z= " << v.z << endl;

void testvector()
{
	cout << "testing vector..." << endl;

	vector2f u, v;
	vector2i dim(10,20);	
	vector3f v3(1.5f,2.0f,1), v4(5,7,-3);
	vector3i v3i(10,20,40);
	
	u.x= 1;
	u.y= 0;
	
	v(0,1);
	v.reflect(u);
	outputv2(v);

	v3= -v3;
	float l= v3.length();
	v4.reflect(v3);
	float d= v3*v4;
	v3+=vector3f(6,6,6);
	v3= v4 + v3;
	outputv3(v3);

	u= u + v;
	v3*=5;
	
	v3i= vector3i(7,7,7);
	v3i[1]= 1337;
	outputv3(v3i);
	cout << endl;
}


void testmatrix()
{
	cout << "testing matrix..." << endl;

	matrix2 m2;
	vector2f u(2,3), v(1,5);
	m2.identity();	
	m2.transpose();
	m2*=5;
	float t= m2.trace();
	u= m2*v;
	m2.invert();
	
	matrix3 m3;
	m3.identity();
	t= m3.determinant();

	matrix3::rotate_euler(m3, (float)(90*Math::DEG2RAD), 0,0);
	cout << "x axis euler:" << endl;
	outputmat3(m3);
	matrix3::rotate_euler(m3, 0, (float)(90*Math::DEG2RAD), 0);
	cout << "y axis euler:" << endl;
	outputmat3(m3);
	matrix3::rotate_euler(m3, 0, 0, (float)(90*Math::DEG2RAD));
	cout << "z axis euler:" << endl;
	outputmat3(m3);

	matrix3 m3bernstein(1,0,0,
		                0,2,0,
						0,-1,1);
	m3= m3bernstein;
	m3bernstein.invert();

	m3*=m3bernstein;

	cout << "3x3 bernstein test" << endl;
	outputmat3(m3);

	matrix4 m4bernstein(1,0,0,0,
						0,3,0,0,
						0,-3,3,0,
						0,1,-3,1);

	matrix4 m4= m4bernstein;

	m4bernstein.invert();
	m4*=m4bernstein;
	cout << "4x4 bernstein test" << endl;
	outputmat4(m4);
}


void testquaternion()
{
	cout << "testing quaternion..." << endl;
	
	matrix3 m3= quaternion::from_axis_angle(vector3f(1,0,0), (float)(45*Math::DEG2RAD));
	cout << "axis angle to quaternion to matrix:" << endl;
	outputmat3(m3);

	// quaternion commutivity test
	quaternion p(1, 2, 3, 4), q(2, 3, 4, 5);
	quaternion ans= p*q;
	cout << "p*q= " << ans.s << " " << ans.v.x << " " << ans.v.y << " " << ans.v.z << endl;  
	ans= q*p;
	cout << "q*p= " << ans.s << " " << ans.v.x << " " << ans.v.y << " " << ans.v.z << endl;  	

	// quaternion log test
	p= quaternion(0,0,0,1);
	q= quaternion(0,0,1,0);
	q.invert();
	ans= (p*q).log();
	cout << "log(p*q)= " << ans.s << " " << ans.v.x << " " << ans.v.y << " " << ans.v.z << endl;  

	// summation test
	quaternion a1(1,0,0,0), a2(0,1,0,0), a3(0,0,1,0), a4(0,0,0,1);
	quaternion x(6,0,0,4);
	quaternion b1(0,1,0,0), b2(0,0,1,0), b3(0,1,0,0), b4(0,0,1,0);
	ans= a1*x*b1 + a2*x*b2 + a3*x*b3 + a4*x*b4;
	cout << "sum(axb)= " << ans.s << " " << ans.v.x << " " << ans.v.y << " " << ans.v.z << endl;  
	cout << endl;
}


void testrect()
{
	cout << "testing rect..." << endl;

	// rect test
	Rectf  r(0,0, 20,20);
	r.Contains(10, 10);
	r.Contains(-1, 10);
	r.Contains(20, 10);	
}


void testfloat()
{
	cout << "testing float..." << endl;
	cout << "float equal test: " << Math::equal(1.5, 1.5) << endl;

	float array[4]= {1, 2, 3, 4};
	float x= 10;	
	cout << "horner's rule eval: " << Math::HornerEval<float>(x, array, 4) << endl;	
	cout << "normal eval: " << x*x*x*array[0] + x*x*array[1] + x*array[2] + array[3] << endl;

	cout << "1.2 rounds to: " << Math::round(1.2f) << endl;
	cout << "10.7 rounds to: " << Math::round(10.7f) << endl;
}


void testrand()
{
	cout << "testing rand..." << endl;

	random *r= new random(55);
	
	cout << "max val= " << r->max() << endl;
	for (int i= 0; i < 10; i++)
	{
		r->get();
		cout << r->current() << "  " << r->currentReal() << endl;
	}
	delete r;
}


int main()
{		
	testvector();
	testmatrix();
	testquaternion();	
	testfloat();
	testrand();		
	return 0;
}