
//work around for bug in old eigen version... -.-

#ifndef EIGENWRAPPER_HPP
#define EIGENWRAPPER_HPP

#include <eigen3/Eigen/Dense>


template<unsigned int dataDim, typename dataT>
class myAbs{
private:
  dataT X;
public:
  typedef typename Eigen::MatrixBase<dataT>::RealScalar real_type;

  myAbs(dataT Y):X(Y){}

  operator real_type(){
    return X.norm();
  }
};

template<typename dataT>
class myAbs<1u,dataT>{
private:
  dataT X;
public:
  myAbs(dataT Y):X(Y){}

  typedef typename Eigen::MatrixBase<dataT>::RealScalar real_type;

  operator real_type(){
    return abs(X);
  }
};




#endif
