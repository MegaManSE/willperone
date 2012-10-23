/****************************************
 * Quaternion class
 * By Will Perone (will.perone@gmail.com)
 * Original: 12-09-2003  
 * Revised:  27-09-2003
 *           22-11-2003
 *           10-12-2003
 *           15-01-2004
 *           16-04-2004
 *           07-29-2011    added corrections from website
 *           22-12-2011    added correction to *= operator, thanks Steve Rogers
 *           22-10-2012    fixed ctor from euler angles & added non windows platform fixes, thanks to Art Golf
 *   
 * Dependancies: My 4x4 matrix class
 * 
 * © 2003, This code is provided "as is" and you can use it freely as long as 
 * credit is given to Will Perone in the application it is used in
 *
 * Notes:  
 * if |q|=1 then q is a unit quaternion
 * if q=(0,v) then q is a pure quaternion 
 * if |q|=1 then q conjugate = q inverse
 * if |q|=1 then q= [cos(angle), u*sin(angle)] where u is a unit vector 
 * q and -q represent the same rotation 
 * q*q.conjugate = (q.length_squared, 0) 
 * ln(cos(theta),sin(theta)*v)= ln(e^(theta*v))= (0, theta*v)
 ****************************************/

#pragma once


#include "matrix4.h"
#include "assert.h"


struct quaternion
{
	#ifdef WIN32
	union {
		struct {
			float    s; //!< the real component
			vector3f v; //!< the imaginary components
		};

		struct { float elem[4]; }; //! the raw elements of the quaternion
	};
	#else
	float    s; //!< the real component
	vector3f v; //!< the imaginary components
	#endif
	

	//! ctors
	quaternion() {}
	quaternion(float real, float x, float y, float z): s(real), v(x,y,z) {}
	quaternion(float real, const vector3f &i): s(real), v(i) {}

	//! from 3 euler angles
	quaternion(float theta_z, float theta_y, float theta_x)//float heading, float attitude, float bank) 
	{
		float cos_z_2 = cosf(0.5f*theta_z);
		float cos_y_2 = cosf(0.5f*theta_y);
		float cos_x_2 = cosf(0.5f*theta_x);

		float sin_z_2 = sinf(0.5f*theta_z);
		float sin_y_2 = sinf(0.5f*theta_y);
		float sin_x_2 = sinf(0.5f*theta_x);

		// and now compute quaternion
		s   = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;
		v.x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
		v.y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
		v.z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;
	}
	
	//! from 3 euler angles 
	quaternion(const vector3f &angles)
	{	
		float cos_z_2 = cosf(0.5f*angles.z);
		float cos_y_2 = cosf(0.5f*angles.y);
		float cos_x_2 = cosf(0.5f*angles.x);

		float sin_z_2 = sinf(0.5f*angles.z);
		float sin_y_2 = sinf(0.5f*angles.y);
		float sin_x_2 = sinf(0.5f*angles.x);

		// and now compute quaternion
		s   = cos_z_2*cos_y_2*cos_x_2 + sin_z_2*sin_y_2*sin_x_2;
		v.x = cos_z_2*cos_y_2*sin_x_2 - sin_z_2*sin_y_2*cos_x_2;
		v.y = cos_z_2*sin_y_2*cos_x_2 + sin_z_2*cos_y_2*sin_x_2;
		v.z = sin_z_2*cos_y_2*cos_x_2 - cos_z_2*sin_y_2*sin_x_2;
	} 
		
	//! basic operations
	quaternion &operator =(const quaternion &q)		
	{	s= q.s; v= q.v;	return *this;		}

	const quaternion operator +(const quaternion &q) const	
	{	return quaternion(s+q.s, v+q.v);	}

	const quaternion operator -(const quaternion &q) const	
	{	return quaternion(s-q.s, v-q.v);	}

	const quaternion operator *(const quaternion &q) const	
	{	return quaternion(s*q.s - v*q.v,
						  v.y*q.v.z - v.z*q.v.y + s*q.v.x + v.x*q.s,
						  v.z*q.v.x - v.x*q.v.z + s*q.v.y + v.y*q.s,
						  v.x*q.v.y - v.y*q.v.x + s*q.v.z + v.z*q.s);
	}

	/*const quaternion operator /(const quaternion &q) const
	{
		float denominator = q.length();
       
		// unsure if this is correct 
        return quaternion((s*q.s + v*q.v)/denominator,  
                          (-s*q.v.x + v.x*q.s - v.y*q.v.z + v.z*q.v.y)/denominator,  
                          (-s*q.v.y + v.x*q.v.z + v.y*q.s - v.z*q.v.x)/denominator,  
                          (-s*q.v.z - v.x*q.v.y + v.y*q.v.x + v.z*q.s)/denominator);  
	}*/
	const quaternion operator /(const quaternion &q) const	
	{	
			quaternion p(q); 
			p.invert(); 
			return *this * p;
	}

	const quaternion operator *(float scale) const			
		{	return quaternion(s*scale,v*scale);		}

	const quaternion operator /(float scale) const
		{	return quaternion(s/scale,v/scale);		}

	const quaternion operator -() const						
		{	return quaternion(-s, -v);				}
	
	const quaternion &operator +=(const quaternion &q)		
		{	v+=q.v; s+=q.s; return *this;			}

	const quaternion &operator -=(const quaternion &q)		
		{	v-=q.v; s-=q.s; return *this;			}

	const quaternion &operator *=(const quaternion &q)		
	{	
			float x= v.x, y= v.y, z= v.z, sn= s*q.s - v*q.v;
			v.x= y*q.v.z - z*q.v.y + s*q.v.x + x*q.s;
			v.y= z*q.v.x - x*q.v.z + s*q.v.y + y*q.s;
			v.z= x*q.v.y - y*q.v.x + s*q.v.z + z*q.s;
			s= sn;
			return *this;
	}
	
	const quaternion &operator *= (float scale)			
		{	v*=scale; s*=scale; return *this;		}

	const quaternion &operator /= (float scale)			
		{	v/=scale; s/=scale; return *this;		}
	

	//! gets the length of this quaternion
	float length() const
	{	return (float)sqrt(s*s + v*v);   }

	//! gets the squared length of this quaternion
	float length_squared() const
	{	return (float)(s*s + v*v);   }

	//! normalizes this quaternion
	void normalize()
	{	*this/=length();	}

	//! returns the normalized version of this quaternion
	quaternion normalized() const
	{   return  *this/length();  }

	//! computes the conjugate of this quaternion
	void conjugate()
	{	v=-v;   }

	//! inverts this quaternion
	void invert()
	{	conjugate(); *this/=length_squared(); 	}
	
	//! returns the logarithm of a quaternion = v*a where q = [cos(a),v*sin(a)]
	quaternion log() const
	{
		float a = (float)acos(s);
		float sina = (float)sin(a);
		quaternion ret;

		ret.s = 0;
		if (sina > 0)
		{
			ret.v.x = a*v.x/sina;
			ret.v.y = a*v.y/sina;
			ret.v.z = a*v.z/sina;
		} else {
			ret.v.x= ret.v.y= ret.v.z= 0;
		}
		return ret;
	}

	//! returns e^quaternion = exp(v*a) = [cos(a),vsin(a)]
	quaternion exp() const
	{
		float a = (float)v.length();
		float sina = (float)sin(a);
		float cosa = (float)cos(a);
		quaternion ret;

		ret.s = cosa;
		if (a > 0)
		{
			ret.v.x = sina * v.x / a;
			ret.v.y = sina * v.y / a;
			ret.v.z = sina * v.z / a;
		} else {
			ret.v.x = ret.v.y = ret.v.z = 0;
		}
		return ret;
	}

	//! casting to a 4x4 isomorphic matrix for right multiplication with vector
	operator matrix4() const
	{			
		return matrix4(s,  -v.x, -v.y,-v.z,
			           v.x,  s,  -v.z, v.y,
					   v.y, v.z,    s,-v.x,
					   v.z,-v.y,  v.x,   s);
	}

	//! casting to a 4x4 isomorphic matrix for left multiplication with vector
	/*operator matrix4() const
	{			
		return matrix4(   s, v.x, v.y, v.z,
			           -v.x,   s,-v.z, v.y,
					   -v.y, v.z,   s,-v.x,
					   -v.z,-v.y, v.x,   s);
	}*/

	//! casting to 3x3 rotation matrix
	operator matrix3() const
	{
		Assert(length() > 0.9999 && length() < 1.0001, "quaternion is not normalized");		
		return matrix3(1-2*(v.y*v.y+v.z*v.z), 2*(v.x*v.y-s*v.z),   2*(v.x*v.z+s*v.y),   
					   2*(v.x*v.y+s*v.z),   1-2*(v.x*v.x+v.z*v.z), 2*(v.y*v.z-s*v.x),   
					   2*(v.x*v.z-s*v.y),   2*(v.y*v.z+s*v.x),   1-2*(v.x*v.x+v.y*v.y));
	}

	//! computes the dot product of 2 quaternions
	static inline float dot(const quaternion &q1, const quaternion &q2) 
	{   return q1.v*q2.v + q1.s*q2.s;  }

	//! linear quaternion interpolation
	static quaternion lerp(const quaternion &q1, const quaternion &q2, float t) 
	{	return (q1*(1-t) + q2*t).normalized();	}

	//! spherical linear interpolation
	static quaternion slerp(const quaternion &q1, const quaternion &q2, float t) 
	{
		quaternion q3;
		float dot = quaternion::dot(q1, q2);

		/*	dot = cos(theta)
			if (dot < 0), q1 and q2 are more than 90 degrees apart,
			so we can invert one to reduce spinning	*/
		if (dot < 0)
		{
			dot = -dot;
			q3 = -q2;
		} else q3 = q2;
		
		if (dot < 0.95f)
		{
			float angle = acosf(dot);
			return (q1*sinf(angle*(1-t)) + q3*sinf(angle*t))/sinf(angle);
		} else // if the angle is small, use linear interpolation								
			return lerp(q1,q3,t);		
	}

	//! This version of slerp, used by squad, does not check for theta > 90.
	static quaternion slerpNoInvert(const quaternion &q1, const quaternion &q2, float t) 
	{
		float dot = quaternion::dot(q1, q2);

		if (dot > -0.95f && dot < 0.95f)
		{
			float angle = acosf(dot);			
			return (q1*sinf(angle*(1-t)) + q2*sinf(angle*t))/sinf(angle);
		} else  // if the angle is small, use linear interpolation								
			return lerp(q1,q2,t);			
	}

	//! spherical cubic interpolation
	static quaternion squad(const quaternion &q1,const quaternion &q2,const quaternion &a,const quaternion &b,float t)
	{
		quaternion c= slerpNoInvert(q1,q2,t),
			       d= slerpNoInvert(a,b,t);		
		return slerpNoInvert(c,d,2*t*(1-t));
	}
	
	//! calculate inner quad point for spherical cubic interpolation 2
	static quaternion innerQuadPoint(const quaternion &q0,const quaternion &q1,const quaternion &q2)
	{
		quaternion vTemp, vTemp1, vTemp2;
		quaternion qResult;
		
		vTemp = q1;
		vTemp.invert();
		
		vTemp1 = (vTemp * q0);
		vTemp1 = vTemp1.log();

		vTemp2 = (vTemp * q2);
		vTemp2 = vTemp2.log();
		
		vTemp = (vTemp1 + vTemp2) * (-0.25f);
		
		qResult = q1 * vTemp.exp();
		
		return qResult;
	}

	//! spherical cubic interpolation alternate approach
	//! this method is smoother than the original squad but only works when the angles are less than 90 degrees apart
	static quaternion squad2(const quaternion &q1,const quaternion &q2,const quaternion &q3,const quaternion &q4,float t)
	{
		quaternion a= innerQuadPoint(q1,q2,q3),
			       b= innerQuadPoint(q2,q3,q4);
		return slerpNoInvert(slerpNoInvert(q2,q3,t),slerpNoInvert(a,b,t),2.0f*t*(1.0f-t));
	}

	//! Shoemake-Bezier interpolation using De Castlejau algorithm
	static quaternion bezier(const quaternion &q1,const quaternion &q2,const quaternion &a,const quaternion &b,float t)
	{
		// level 1
		quaternion q11= slerpNoInvert(q1,a,t),
			       q12= slerpNoInvert(a,b,t),
			       q13= slerpNoInvert(b,q2,t);		
		// level 2 and 3
		return slerpNoInvert(slerpNoInvert(q11,q12,t), slerpNoInvert(q12,q13,t), t);
	}

	//! Given 3 quaternions, qn-1,qn and qn+1, calculate a control point to be used in spline interpolation
	static quaternion spline(const quaternion &qnm1,const quaternion &qn,const quaternion &qnp1)
	{
		quaternion qni(qn.s, -qn.v);		
		
		return qn * (( (qni*qnm1).log()+(qni*qnp1).log() )/-4).exp();
	}

	//! converts from a normalized axis - angle pair rotation to a quaternion
	static inline quaternion from_axis_angle(const vector3f &axis, float angle)
	{	return quaternion(cosf(angle/2), axis*sinf(angle/2)); 	}

	//! returns the axis and angle of this unit quaternion
	void to_axis_angle(vector3f &axis, float &angle) const
	{
		angle = acosf(s);

		// pre-compute to save time
		float sinf_theta_inv = 1.0/sinf(angle);

		// now the vector
		axis.x = v.x*sinf_theta_inv;
		axis.y = v.y*sinf_theta_inv;
		axis.z = v.z*sinf_theta_inv;

		// multiply by 2
		angle*=2;
	}

	//! rotates v by this quaternion (quaternion must be unit)
	vector3f rotate(const vector3f &v)
	{   
		quaternion V(0, v);
	    quaternion conjugate(*this);
		conjugate.conjugate();
		return (*this * V * conjugate).v;
	}

	//! returns the euler angles from a rotation quaternion
	vector3f euler_angles(bool homogenous=true) const
	{
		float sqw = s*s;    
		float sqx = v.x*v.x;    
		float sqy = v.y*v.y;    
		float sqz = v.z*v.z;    

		vector3f euler;
		if (homogenous) {
			euler.x = atan2f(2.f * (v.x*v.y + v.z*s), sqx - sqy - sqz + sqw);    		
			euler.y = asinf(-2.f * (v.x*v.z - v.y*s));
			euler.z = atan2f(2.f * (v.y*v.z + v.x*s), -sqx - sqy + sqz + sqw);    
		} else {
			euler.x = atan2f(2.f * (v.z*v.y + v.x*s), 1 - 2*(sqx + sqy));
			euler.y = asinf(-2.f * (v.x*v.z - v.y*s));
			euler.z = atan2f(2.f * (v.x*v.y + v.z*s), 1 - 2*(sqy + sqz));
		}
		return euler;
	}
};


