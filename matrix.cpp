#include <mpi.h>
#include <iostream>
#define g 1500

using namespace std;

int main(int argc, char *argv[]) {
  int i, j, k;
  int rank, size;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  int xBegin[size];
  int xSize[size];

  int h;

  for (i = 0; i < size; i++) {
    xBegin[i] = int(g / size) * i + min(g % size, i);
    xSize[i] = int(g / size) + int(g % size > i);
  }

  if (rank == 0) {
    h = g;
  } else {
    h = xSize[rank];
  }

  float matA[h][g];
  float matB[g][g];
  float matC[h][g];

  double startTime, endTime, elapsedTime;
  double startTimeNet, endTimeNet, elapsedTimeNet;

  elapsedTime = 0.0;
  elapsedTimeNet = 0.0;

  startTime = MPI_Wtime();
  // Generate
  if (rank == 0) {
    for (i = 0; i < g; i++) {
      for (j = 0; j < g; j++) {
        matA[i][j] = i * g + j;
        if (i == j) {
          matB[i][j] = 1;
        } else {
          matB[i][j] = 0;
        }
      }
    }
  }

  // Send matA
  startTimeNet = MPI_Wtime();
  if (rank == 0) {
    for (i = 1; i < size; i++) {
      MPI_Send(&matA[xBegin[i]][0], xSize[i] * g, MPI_FLOAT, i, 0,
               MPI_COMM_WORLD);
    }
  } else {
    MPI_Recv(&matA, xSize[rank] * g, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
  }

  // Send matB
  MPI_Bcast(&matB, g * g, MPI_FLOAT, 0, MPI_COMM_WORLD);
  endTimeNet = MPI_Wtime();
  elapsedTimeNet += endTimeNet - startTimeNet;

#if 0  // Send Check
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank == 0) {
    for (i = 0; i < xSize[0]; i++) {
      for (j = 0; j < g; j++) {
        cout << matA[i][j] << " ";
      }
      cout << endl;
    }

    for (i = 0; i < g; i++) {
      for (j = 0; j < g; j++) {
        cout << matB[i][j] << " ";
      }
      cout << endl;
    }
  }
  MPI_Barrier(MPI_COMM_WORLD);
  if (rank != 0) {
    for (i = 0; i < h; i++) {
      for (j = 0; j < g; j++) {
        cout << matA[i][j] << " ";
      }
      cout << endl;
    }

    for (i = 0; i < g; i++) {
      for (j = 0; j < g; j++) {
        cout << matB[i][j] << " ";
      }
      cout << endl;
    }
  }
#endif

  // startTime = MPI_Wtime();
  for (i = 0; i < xSize[rank]; i++) {
    for (j = 0; j < g; j++) {
      float sum = 0.0;
      for (k = 0; k < g; k++) {
        sum += matA[i][k] * matB[k][j];
      }
      matC[i][j] = sum;
    }
  }
  // endTime = MPI_Wtime();

#if 0  // debug
  if (rank == 2) {
    for (i = 0; i < h; i++) {
      for (j = 0; j < g; j++) {
        cout << matC[i][j] << " ";
      }
      cout << endl;
    }
  }
#endif

  float results[xSize[0] * size][g];
  startTimeNet = MPI_Wtime();
  MPI_Gather(&matC, xSize[0] * g, MPI_FLOAT, &results, xSize[0] * g, MPI_FLOAT,
             0, MPI_COMM_WORLD);
  endTimeNet = MPI_Wtime();
  elapsedTimeNet += endTimeNet - startTimeNet;

  if (rank == 0) {
    for (k = 1; k < size; k++) {
      for (i = 0; i < xSize[k]; ++i) {
        for (j = 0; j < g; ++j) {
          matC[xBegin[k] + i][j] = results[xSize[0] * k][j];
        }
      }
    }

#if 0  // debug
    for (i = 0; i < g; ++i) {
      for (j = 0; j < g; ++j) {
        cout << matC[i][j] << " ";
      }
      cout << endl;
    }
#endif
  }
  endTime = MPI_Wtime();

  double maxElapsedTime;
  double maxElapsedTimeNet;
  elapsedTime = endTime - startTime;
  MPI_Reduce(&elapsedTime, &maxElapsedTime, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);
  MPI_Reduce(&elapsedTimeNet, &maxElapsedTimeNet, 1, MPI_DOUBLE, MPI_MAX, 0,
             MPI_COMM_WORLD);

  if (rank == 0) {
    cout << maxElapsedTime << endl;
    cout << maxElapsedTimeNet << endl;
  }

  MPI_Finalize();
  return 0;
}
