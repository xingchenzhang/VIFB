/*! \file

\verbatim

Copyright (c) 2004, Sylvain Paris and Francois Sillion
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
    notice, this list of conditions and the following disclaimer.
    
    * Redistributions in binary form must reproduce the above
    copyright notice, this list of conditions and the following
    disclaimer in the documentation and/or other materials provided
    with the distribution.

    * Neither the name of ARTIS, GRAVIR-IMAG nor the names of its
    contributors may be used to endorse or promote products derived
    from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

\endverbatim


 *  This file contains code made by Sylvain Paris under supervision of
 * François Sillion for his PhD work with <a
 * href="http://www-artis.imag.fr">ARTIS project</a>. ARTIS is a
 * research project in the GRAVIR/IMAG laboratory, a joint unit of
 * CNRS, INPG, INRIA and UJF.
 *
 *  Use <a href="http://www.stack.nl/~dimitri/doxygen/">Doxygen</a>
 * with DISTRIBUTE_GROUP_DOC option to produce an nice html
 * documentation.
 *
 *  The file defines several common mathematical functions.
 */

#ifndef __MATH_TOOLS__
#define __MATH_TOOLS__




#include <cmath>
#include <cstdlib>
#include <ctime>

#include <algorithm>
#include <limits>
#include <numeric>
#include <vector>

#ifndef WITHOUT_LIMITS
#include <limits>
#endif

#include "msg_stream.h"


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace Math_tools{

  // #########################################################
  
  template<typename Vec2,typename Real>
  inline void barycentric_coordinates(const Vec2& A,
				      const Vec2& B,
				      const Vec2& C,
				      const Vec2& M,
				      Real* const a,
				      Real* const b,
				      Real* const c);
  
  // #########################################################
  
  inline void init_random() {std::srand(std::time(NULL));}

  template<typename Real>
  inline Real random(const Real min_value_included,
		     const Real max_value_included);

  // #########################################################
  
  template<typename ValueIterator>
  inline double entropy(ValueIterator begin,
			ValueIterator end);

  // #########################################################
  
  template<typename Real,typename RealDummy1,typename RealDummy2>
  inline Real clamp(const RealDummy1 min_value,
		    const RealDummy2 max_value,
		    const Real x);

  // #########################################################
  
  template<typename Real>
  inline Real square(const Real x);


  // #########################################################
  
  template<unsigned int N,typename Real>
  inline Real power(const Real x);

  // #########################################################
  
  template<unsigned int N,typename Real>
  inline Real power2(Real x);
  
  // #########################################################
  
  template<typename Real,typename RealDummy1,typename RealDummy2>
  inline Real smooth_step(const RealDummy1 min_value,
			  const RealDummy2 max_value,
			  const Real x);

  // #########################################################
  
  template<typename Array,typename Real>
  inline
  typename Array::value_type
  bilinear_interpolation(const Array& array,
			 const Real x,
			 const Real y);

  // #########################################################
  
  template<typename Array,typename Real>
  inline
  typename Array::value_type
  trilinear_interpolation(const Array& array,
			  const Real x,
			  const Real y,
			  const Real z);

  // #########################################################
  
  template<typename Array,typename Real>
  inline
  typename Array::value_type
  bicubic_interpolation(const Array& array,
			const Real x,
			const Real y);

  // #########################################################
  
  template<typename Array,typename Real>
  inline
  typename Array::value_type
  tricubic_interpolation(const Array& array,
			 const Real x,
			 const Real y,
			 const Real z);

  // #########################################################  

  inline double degree_to_radian(const double d);
  
  inline double radian_to_degree(const double r);


  // #########################################################  
  
  template<typename Real>
  inline Real sign(const Real r);

  // #########################################################  

  template<typename Iterator>
  inline typename std::iterator_traits<Iterator>::value_type
  mean(Iterator begin,Iterator end);

  template<typename Iterator>
  inline double standard_deviation(Iterator begin,Iterator end);

  template<typename Iterator>
  inline double standard_deviation(Iterator begin,Iterator end,
				   typename std::iterator_traits<Iterator>::value_type* const mean_value);

  template<typename Iterator>
  inline typename std::iterator_traits<Iterator>::value_type
  median(Iterator begin,Iterator end,
	 const float position = 0.5);

   // #########################################################

  template<typename BoolArray,typename ScalarArray>
  inline void
  compute_distance_field(const BoolArray&   input,
			 ScalarArray* const result,
			 const typename BoolArray::value_type out_value = typename BoolArray::value_type(),
			 const int          window = 5);
  
  
  // #########################################################  

#ifndef WITHOUT_LIMITS    

  template<typename T> inline bool is_quiet_NaN(T x);

  template<typename T> inline bool is_signaling_NaN(T x);

  template<typename T> inline bool is_NaN(T x);
  
#endif


  
/*
  
  #############################################
  #############################################
  #############################################
  ######                                 ######
  ######   I M P L E M E N T A T I O N   ######
  ######                                 ######
  #############################################
  #############################################
  #############################################
  
*/


  
  
  template<typename Vec2,typename Real>
  void barycentric_coordinates(const Vec2& A,
			       const Vec2& B,
			       const Vec2& C,
			       const Vec2& M,
			       Real* const a,
			       Real* const b,
			       Real* const c){

    const Real Xa = A[0] - C[0];
    const Real Ya = A[1] - C[1];

    const Real Xb = B[0] - C[0];
    const Real Yb = B[1] - C[1];

    const Real Xm = M[0] - C[0];
    const Real Ym = M[1] - C[1];

    const Real det = Xa*Yb - Ya*Xb;

    if (det != 0){
      
      const Real inv_det = 1.0/det;

      *a = inv_det * (Yb*Xm  - Xb*Ym);
      *b = inv_det * (-Ya*Xm + Ym*Xa);
      *c = 1.0 - (*a) - (*b);
    }
    else{

      Message::error
	<<"barycentric_coordinates: colinear points not yet implemented"
	<<Message::done;
      
    }
  }


  

  template<typename Real>
  Real random(const Real min_value_included,
	      const Real max_value_included){
    
    const double delta = static_cast<double>(max_value_included - min_value_included);
    
    return min_value_included + static_cast<Real>(delta*rand()/(RAND_MAX));
  }



  
  template<typename ValueIterator>
  double entropy(ValueIterator begin,
		 ValueIterator end){
    
    double sum = 0;
    for(ValueIterator i=begin;i!=end;i++){
      if (*i>=0){
	sum+=*i;
      }
      else{
	Message::error<<"entropy: non-positive value"<<Message::done;
      }
    }

    if (sum==0){
      Message::error<<"entropy: sum==0"<<Message::done;
    }

    double result = 0;
    for(ValueIterator i=begin;i!=end;i++){
      const double p = (*i)/sum;
      result -= (p==0) ? 0 : p*log(p);
    }

    return result;
  }

  

  template<typename Real,typename RealDummy1,typename RealDummy2>
  inline Real clamp(const RealDummy1 min_value,
		    const RealDummy2 max_value,
		    const Real x){
    
    return std::max(std::min(x,static_cast<Real>(max_value)),static_cast<Real>(min_value));
  }

  
  template<typename Real>
  inline Real square(const Real x){
    
    return x*x;
  }

  
  template<unsigned int N,typename Real>
  inline Real power(const Real x){

    Real y = x;
    for(unsigned int i=1;i<N;i++){
      y *= x;
    }

    return y;
  }

  template<unsigned int N,typename Real>
  inline Real power2(Real x){

    for(unsigned int i=1;i<N;i++){
      x *= x;
    }

    return x;
  }

  
  template<typename Real,typename RealDummy1,typename RealDummy2>
  inline Real smooth_step(const RealDummy1 min_value,
			  const RealDummy2 max_value,
			  const Real x){
    
    const Real rm = static_cast<Real>(min_value);
    const Real rM = static_cast<Real>(max_value);

    if (x<=rm){
      return 0;
    }

    if (x>=rM){
      return 1;
    }
    
    const Real delta = rM - rm;

    const Real alpha = 1.0 - (x-rm) / delta;

    const Real tmp = 1.0 - alpha * alpha;
    
    return tmp*tmp;
  }




  template<typename Array,typename Real>
  inline
  typename Array::value_type
  bilinear_interpolation(const Array& array,
			 const Real x,
			 const Real y){

    typedef unsigned int size_type;
    typedef float        real_type;

    const size_type x_size = array.x_size();
    const size_type y_size = array.y_size();
    
    const size_type x_index  = clamp(0,x_size-1,static_cast<size_type>(x));
    const size_type xx_index = clamp(0,x_size-1,x_index+1);
    
    const size_type y_index  = clamp(0,y_size-1,static_cast<size_type>(y));
    const size_type yy_index = clamp(0,y_size-1,y_index+1);

    const real_type x_alpha = x - x_index;
    const real_type y_alpha = y - y_index;

    return
      (1.0f-x_alpha) * (1.0f-y_alpha) * array(x_index, y_index) +
      x_alpha        * (1.0f-y_alpha) * array(xx_index,y_index) +
      (1.0f-x_alpha) * y_alpha        * array(x_index, yy_index) +
      x_alpha        * y_alpha        * array(xx_index,yy_index);
      
  }

  template<typename Array,typename Real>
  inline
  typename Array::value_type
  trilinear_interpolation(const Array& array,
			  const Real x,
			  const Real y,
			  const Real z){

    typedef unsigned int size_type;
    typedef float        real_type;

    const size_type x_size = array.x_size();
    const size_type y_size = array.y_size();
    const size_type z_size = array.z_size();
    
    const size_type x_index  = clamp(0,x_size-1,static_cast<size_type>(x));
    const size_type xx_index = clamp(0,x_size-1,x_index+1);
    
    const size_type y_index  = clamp(0,y_size-1,static_cast<size_type>(y));
    const size_type yy_index = clamp(0,y_size-1,y_index+1);

    const size_type z_index  = clamp(0,z_size-1,static_cast<size_type>(z));
    const size_type zz_index = clamp(0,z_size-1,z_index+1);

    const real_type x_alpha = x - x_index;
    const real_type y_alpha = y - y_index;
    const real_type z_alpha = z - z_index;

    return
      (1.0f-x_alpha) * (1.0f-y_alpha) * (1.0f-z_alpha) * array(x_index, y_index, z_index) +
      x_alpha        * (1.0f-y_alpha) * (1.0f-z_alpha) * array(xx_index,y_index, z_index) +
      (1.0f-x_alpha) * y_alpha        * (1.0f-z_alpha) * array(x_index, yy_index,z_index) +
      x_alpha        * y_alpha        * (1.0f-z_alpha) * array(xx_index,yy_index,z_index) +
      (1.0f-x_alpha) * (1.0f-y_alpha) * z_alpha        * array(x_index, y_index, zz_index) +
      x_alpha        * (1.0f-y_alpha) * z_alpha        * array(xx_index,y_index, zz_index) +
      (1.0f-x_alpha) * y_alpha        * z_alpha        * array(x_index, yy_index,zz_index) +
      x_alpha        * y_alpha        * z_alpha        * array(xx_index,yy_index,zz_index);
  }


  namespace Cubic_Interpolation{

    inline double cubic_positive_part(const double x){
      const double p = std::max(x,0.0);
      return p*p*p;
    }

    inline double cubic_polynom(const double x){
      return 1.0 / 6.0 * (cubic_positive_part(x+2.0)
			  - 4.0 * cubic_positive_part(x+1.0)
			  + 6.0 * cubic_positive_part(x)
			  - 4.0 * cubic_positive_part(x-1.0));
    }
    
  }
  
  
  template<typename Array,typename Real>
  inline
  typename Array::value_type
  bicubic_interpolation(const Array& array,
			const Real x,
			const Real y){

    using namespace std;
    using namespace Cubic_Interpolation;
    
    typedef unsigned int size_type;
    typedef double       real_type;

    const size_type x_size = array.x_size();
    const size_type y_size = array.y_size();

    vector<size_type> x_index(4);
    vector<size_type> y_index(4);
    
    const size_type x_floor = clamp(0,x_size-1,static_cast<size_type>(x));
    const size_type y_floor = clamp(0,y_size-1,static_cast<size_type>(y));

    const real_type x_alpha = x - x_floor;
    const real_type y_alpha = y - y_floor;
    
    for(int i=-1;i<=2;i++){
      x_index[i+1] = clamp<int>(0,x_size-1,x_floor+i);
      y_index[i+1] = clamp<int>(0,y_size-1,y_floor+i);
    }

    real_type res = 0;
    
    for(int i=-1;i<=2;i++){
      for(int j=-1;j<=2;j++){

	res += array(x_index[i+1],y_index[j+1]) * cubic_polynom(i-x_alpha) * cubic_polynom(j-y_alpha);
      }
    }

    return res;
  }


    
  template<typename Array,typename Real>
  inline
  typename Array::value_type
  tricubic_interpolation(const Array& array,
			 const Real x,
			 const Real y,
			 const Real z){

    using namespace std;
    using namespace Cubic_Interpolation;
    
    typedef unsigned int size_type;
    typedef double       real_type;

    const size_type x_size = array.x_size();
    const size_type y_size = array.y_size();
    const size_type z_size = array.z_size();

    vector<size_type> x_index(4);
    vector<size_type> y_index(4);
    vector<size_type> z_index(4);
    
    const size_type x_floor = clamp(0,x_size-1,static_cast<size_type>(x));
    const size_type y_floor = clamp(0,y_size-1,static_cast<size_type>(y));
    const size_type z_floor = clamp(0,z_size-1,static_cast<size_type>(z));

    const real_type x_alpha = x - x_floor;
    const real_type y_alpha = y - y_floor;
    const real_type z_alpha = z - z_floor;
    
    for(int i=-1;i<=2;i++){
      x_index[i+1] = clamp<int>(0,x_size-1,x_floor+i);
      y_index[i+1] = clamp<int>(0,y_size-1,y_floor+i);
      z_index[i+1] = clamp<int>(0,z_size-1,z_floor+i);
    }

    real_type res = 0;
    
    for(int i=-1;i<=2;i++){
      for(int j=-1;j<=2;j++){
	for(int k=-1;k<=2;k++){

	  res += array(x_index[i+1],y_index[j+1],z_index[k+1]) * cubic_polynom(i-x_alpha) * cubic_polynom(j-y_alpha) * cubic_polynom(k-z_alpha);
	}
      }
    }

    return res;
  }



  
  inline double degree_to_radian(const double d){
    return d*M_PI/180;
  }
  

  inline double radian_to_degree(const double r){
    return r*180/M_PI;
  }


  template<typename Real>
  inline Real sign(const Real r){
    return (r==0) ? 0 : ((r>0) ? 1 : -1);
  }

  template<typename Iterator>
  inline typename std::iterator_traits<Iterator>::value_type
  mean(Iterator begin,Iterator end){

    typedef typename std::iterator_traits<Iterator>::value_type value_type;
    
    return std::accumulate(begin,end,value_type()) / std::distance(begin,end);
  }

  template<typename Iterator>
  double standard_deviation(Iterator begin,Iterator end,
			   typename std::iterator_traits<Iterator>::value_type* const mean_value){

    typedef typename std::iterator_traits<Iterator>::value_type value_type;

    unsigned int size = std::distance(begin,end);
    
    std::vector<double> dist_to_mean(size);

    const value_type m = mean(begin,end);
    
    Iterator i=begin;
    for(std::vector<double>::iterator d=dist_to_mean.begin();i!=end;i++,d++){
      
      const double dist = std::abs(*i-m);
      *d = dist*dist;
    }

    *mean_value = m;
    return std::sqrt(mean(dist_to_mean.begin(),dist_to_mean.end()));
  }

  
  template<typename Iterator>
  double standard_deviation(Iterator begin,Iterator end){
    typename std::iterator_traits<Iterator>::value_type proxy;

    return standard_deviation(begin,end,&proxy);
  }

  
  template<typename Iterator>
  typename std::iterator_traits<Iterator>::value_type
  median(Iterator begin,Iterator end,
	 const float position){

    std::vector<typename std::iterator_traits<Iterator>::value_type> proxy(begin,end);
    std::sort(proxy.begin(),proxy.end());
    return proxy[static_cast<unsigned int>((proxy.size()-1)*position)];
  }


  
  template<typename BoolArray,typename ScalarArray>
  void compute_distance_field(const BoolArray&   input,
			      ScalarArray* const result,
			      const typename BoolArray::value_type out_value,
			      const int          window){

    using namespace std;
    
    typedef typename ScalarArray::value_type real_type;
    typedef unsigned int size_type;
    
    const size_type width  = input.x_size();
    const size_type height = input.y_size();

    ScalarArray& distance = *result;
    
    distance.resize(width,height);

    for(size_type x=0;x<width;x++){
      for(size_type y=0;y<height;y++){

	distance(x,y) = (input(x,y) != out_value) ? 0 : numeric_limits<real_type>::max();
	
      }
    }

    for(int x=window;x+window<static_cast<int>(width);x++){
      for(int y=window;y+window<static_cast<int>(height);y++){

	const real_type d = distance(x,y);

	for(int dx=0;dx<=window;dx++){
	  for(int dy=((dx==0)?1:-window);dy<=window;dy++){
	    
	    const real_type delta = sqrt(1.0*dx*dx+1.0*dy*dy);

	    real_type& r = distance(x+dx,y+dy);

	    r = min(d+delta,r);
	  }
	}
	
      }
    }

    
    for(int x=width-window-1;x>=window;x--){
      for(int y=height-window-1;y>=window;y--){

	const real_type d = distance(x,y);

	for(int dx=0;dx<=window;dx++){
	  for(int dy=((dx==0)?1:-window);dy<=window;dy++){
	    
	    const real_type delta = sqrt(1.0*dx*dx+1.0*dy*dy);

	    real_type& r = distance(x-dx,y-dy);

	    r = min(d+delta,r);
	  }
	}
	
      }
    }    
  }


  
#ifndef WITHOUT_LIMITS    

  template<typename T> bool is_quiet_NaN(T x){

    const unsigned int size = sizeof(T)/sizeof(int);
    
    T ref = std::numeric_limits<T>::quiet_NaN();
    
    int* index_ref = reinterpret_cast<int*>(&ref);
    int* index_cmp = reinterpret_cast<int*>(&x);

    for(unsigned int i=0;i<size;i++){
      if (index_ref[i]!=index_cmp[i]){
	return false;
      }
    }

    return true;
  }

  
  template<typename T> bool is_signaling_NaN(T x){
    
    const unsigned int size = sizeof(T)/sizeof(int);
    
    T ref = std::numeric_limits<T>::signaling_NaN();
    
    int* index_ref = reinterpret_cast<int*>(&ref);
    int* index_cmp = reinterpret_cast<int*>(&x);

    for(unsigned int i=0;i<size;i++){
      if (index_ref[i]!=index_cmp[i]){
	return false;
      }
    }

    return true;
  }


  template<typename T> bool is_NaN(T x){
    return (is_quiet_NaN(x)||is_signaling_NaN(x));
  }
  
#endif

  
} // end of namespace


#endif
