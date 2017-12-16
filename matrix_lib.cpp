#include <mpi.h>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace std;

typedef unsigned int uint;

class matrix {
private:
  uint n;
  uint m;
  vector<vector<float> > v;

public:
  matrix();
  matrix(uint N, uint M);
  float operator()(uint i, uint j);
  void set(uint i, uint j, float val);
  int size_row();
  int size_col();
  friend ostream &operator<<(ostream &stream, matrix ob);

  int lgt();
  float max();

  static matrix arange(uint n, uint m);
  static matrix eye(uint w);
  static matrix ones(uint n, uint m);
  static matrix random_int(uint n, uint m, int min = 0, int max = 100);
  static matrix random_real(uint n, uint m, float min = 0.0, float max = 1.0);
  static matrix zeros(uint n, uint m);
};

matrix::matrix() {
  n = 0;
  m = 0;
}

matrix::matrix(uint N, uint M) {
  n = N;
  m = M;
  v.resize(n);
  for (int i = 0; i < n; i++) {
    v[i].resize(m);
  }
}

float matrix::operator()(uint i, uint j) {
  if (i >= n || j >= m) {
    throw "index out of bounds";
  }
  return v[i][j];
}

void matrix::set(uint i, uint j, float val) {
  if (i >= n || j >= m) {
    throw "index out of bounds";
  }
  v[i][j] = val;
}

int matrix::size_row() {
  return n;
}

int matrix::size_col() {
  return m;
}

ostream &operator<<(ostream &stream, matrix ob) {
  uint i, j;

  int lg = ob.lgt() + 1;

  stream << '[' << endl;
  for (i = 0; i < ob.n; i++) {
    stream << "\t[";
    for (j = 0; j < ob.m; j++) {
      stream << setw(lg) << ob.v[i][j];
    }
    stream << " ]" << endl;
  }
  stream << ']' << endl;
  return stream;
}

int matrix::lgt() {
  float mx = max();
  return int(log10(mx)) + 1;
}

float matrix::max() {
  float m = v[0][0];
  int i, j;

  for (i = 0; i < v.size(); i++) {
    for (j = 0; j < v[i].size(); j++) {
      if (m < v[i][j]) m = v[i][j];
    }
  }
  return m;
}

matrix matrix::arange(uint n, uint m) {
  matrix mat(n, m);
  for (int i = 0; i < n * m; ++i) {
    mat.set(int(i / m), int(i % m), float(i));
  }
  return mat;
}

matrix matrix::eye(uint w) {
  matrix mat(w, w);
  for (int i = 0; i < w; ++i) {
    mat.set(i, i, 1.0);
  }
  return mat;
}

matrix matrix::ones(uint n, uint m) {
  matrix mat(n, m);
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) mat.set(i, j, 1);
  return mat;
}

matrix matrix::random_int(uint n, uint m, int min, int max) {
  matrix mat(n, m);
  int i, j;

  random_device rd;
  mt19937 mt(rd());
  uniform_int_distribution<> rnd(min, max);

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) mat.set(i, j, float(rnd(mt)));
  return mat;
}

matrix matrix::random_real(uint n, uint m, float min, float max) {
  matrix mat(n, m);
  int i, j;

  random_device rd;
  mt19937 mt(rd());
  uniform_real_distribution<> rnd(min, max);

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) mat.set(i, j, rnd(mt));
  return mat;
}

matrix matrix::zeros(uint n, uint m) {
  matrix mat(n, m);
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < m; j++) mat.set(i, j, 0);
  return mat;
}

matrix dot(matrix a, matrix b) {
  matrix ans = matrix(a.size_row(), b.size_col());
  int i, j, k, kmax;

  if (a.size_col() != b.size_row()) {
    throw "cannot dot calc";
  }

  kmax = a.size_col();

  for (i = 0; i < ans.size_row(); i++) {
    for (j = 0; j < ans.size_col(); j++) {
      float s = 0.0;
      for (k = 0; k < kmax; k++) s += a(i, k) * b(k, j);
      ans.set(i, j, s);
    }
  }

  return ans;
}

float inner_product(float a[], float b[], uint length) {
  float sum = 0.0;
  for (int i = 0; i < length; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

int min(int x, int y) {
  if (x < y) {
    return x;
  }
  return y;
}

int main(int argc, char *argv[]) {
  int g = 300;
  int size, rank;
  matrix a, b, c;
  double startTime, endTime, elapsedTime;

#if 1  // initialize
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  a = matrix::arange(g, g);
  b = matrix::eye(g);

  MPI_Barrier(MPI_COMM_WORLD);
#endif

#if 0  // Pre Prarell
  startTime = MPI_Wtime();
  c = dot(a, b);
  endTime = MPI_Wtime();
#endif

#if 1  // Do Pararell
  startTime = MPI_Wtime();
  int i, j, k;

  int xBegin = int(g / size) * rank + min(g % size, rank);
  int xSize = int(g / size) + int(g % size > rank);

  c = matrix(a.size_row(), b.size_col());

  for (i = xBegin; i < xBegin + xSize; i++) {
    for (j = 0; j < g; j++) {
      float sum = 0.0;
      for (k = 0; k < g; k++) {
        sum += a(i, k) * b(k, j);
      }
      c.set(i, j, sum);
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  endTime = MPI_Wtime();


#if 0 // show results
  float result[xSize][g];
  float results[g + size][g];

  for (i = 0; i < xSize; ++i) {
    for (j = 0; j < g; ++j) {
      result[i][j] = c(xBegin + i, j);
    }
  }

  MPI_Gather(&result, xSize * g, MPI_FLOAT, &results, xSize * g, MPI_FLOAT, 0,
      MPI_COMM_WORLD);

  if (rank == 0) {
    cout << "result:" <<endl;
    for (k = 0; k < size; k++) {
      for (i = 0; i < xSize; ++i) {
        if (k != 0 && g / k != xSize && i == xSize) {
          continue;
        }
        for (j = 0; j < g; ++j) {
          cout << results[k * xSize + i][j] << " ";
        }
        cout << endl;
      }
    }
    cout << endl << endl;
  }
#endif // end show result

#endif

#if 1  // finalize
  elapsedTime = endTime - startTime;
  MPI_Finalize();
  if (rank == 0) {
    cout << "elapsed time: " << elapsedTime << endl;
  }
#endif

  return 0;
}
