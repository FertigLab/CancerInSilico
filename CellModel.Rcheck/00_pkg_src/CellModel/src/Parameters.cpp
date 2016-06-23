#include "Parameters.hpp"

#include <iostream>

void Parameters::InitializeRadiusSolver() {

  double theta;  

  for (int i = 0; i <= 3140; i++) {

    theta = (double) i / 1000;  
    m_radius_key.push_back(m_max_radius * pow(M_PI / (2 * M_PI - theta + sin(theta)),0.5));

  }

}

double Parameters::GetTheta(double rad) {

  int imin = 0, imax = 3140;
  int imid;

  while (imin + 1 < imax) {

    imid = (int) ((imin + imax) / 2);

    if (m_radius_key[imid] <= rad) {

      imin = imid;

    } else {

      imax = imid;

    }

  }

  if (rad - m_radius_key[imin] < m_radius_key[imax] - rad) {

    return (double) imin / 1000;

  } else {

    return (double) imax / 1000;

  }  

}

