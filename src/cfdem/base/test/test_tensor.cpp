#include <iostream>
#include "../tensor/tensor.H"
using std::cout;
using std::endl;

int main() {
  int nrow = 3;
  int ncol = 5;
  auto shape2 = base::makeShape2(nrow, ncol);
  auto shape1 = shape2.flatTo1D();
  cout << shape2 << endl;
  cout << shape1 << endl;
  double *pnum = new double[3 * 5];
  std::fill(pnum, pnum + nrow * ncol, 0);
  base::Tensor<1, base::cpu, double> tensor1(pnum, shape1);
  base::CDTensor1 tensor2(pnum, shape1);
  base::CDTensor2 tensor3(pnum, shape2);
  free(pnum);
  pnum = nullptr;
  return 0;
}