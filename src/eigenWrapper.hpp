
//work around for bug in old eigen version... -.-

#ifndef EIGENWRAPPER_HPP
#define EIGENWRAPPER_HPP

#include <eigen3/Eigen/Dense>


template<unsigned int dataDim, typename dataT>
class myAbs{
private:
  dataT X;
public:
  myAbs(dataT Y):X(Y){}

  operator typename Eigen::NumTraits<dataT>::Real(){
    return X.norm();
  }
};

template<typename dataT>
class myAbs<1u,dataT>{
private:
  dataT X;
public:
  myAbs(dataT Y):X(Y){}

  operator typename Eigen::NumTraits<dataT>::Real(){
    return abs(X);
  }
};



#endif
