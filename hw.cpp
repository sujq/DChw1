#include <cstdio>
#include <cmath>
#include <mpi.h>
#include "vector"
#include "iostream"
using namespace std;

//array of numbers
int num = 100;//total number
bool is_prime[512000];
int pace[1432];

int main (int argc, char **argv)
{
    for (int i = 0; i < 512000; i++) is_prime[i] = true;
    is_prime[0] = false;
    is_prime[1] = false;
    //for (int i = 0; i < 716; i++) pace[i] = -1;
    int length = 2*ceil(sqrt(num));
    int all = 0;// total numbers of prime
    int cnt = 0;// numbers of prime in one process
    double total_time = 0;
    double time = 0;
    double start_time = 0;
    double end_time = 0;
    //operation
    int rank, size;
    
    MPI_Init (&argc, &argv);  /* starts MPI */
    
    MPI_Comm_rank (MPI_COMM_WORLD, &rank);    /* get current process id */
    MPI_Comm_size (MPI_COMM_WORLD, &size);    /* get number of processes */
    
    start_time = MPI_Wtime();
    // each part is [left, right)
    int left = rank*length;
    if (left < 2) left = 2;
    
    int right = (rank+1)*length;
    if (right > num+1) right = num;
    
    int k = 0;//index of pace
    if (rank == 0) {
        for (int i = 2; i < length; i++) {
            if (is_prime[i]) {
                for (int j = i+i; j < right; j = j+i) {
                    is_prime[j] = false;
                    //cout << j << " is disabled" << endl;
                }
                pace[k] = i;
                MPI_Bcast(&pace[k], 1, MPI_INT, 0, MPI_COMM_WORLD);
                k++;
            }
        }
        pace[k] = 0;// stop
        MPI_Bcast(&pace[k], 1, MPI_INT, 0, MPI_COMM_WORLD);
    }
    else {
        int start = 0;// real start position
        while (1) {
            MPI_Bcast(&pace[k], 1, MPI_INT, 0, MPI_COMM_WORLD);
            if (pace[k] == 0) break;
            start = left+(pace[k]-left%pace[k])%pace[k];
            if (start > right) break;
            for (int i = start; i < right; i = i+pace[k]) {
                is_prime[i] = false;
                //cout << i << " is disabled" << endl;
            }
            k++;
        }
    }
    for (int i = left; i < right; i++) {
        if (is_prime[i]) cnt++;
    }
    end_time = MPI_Wtime();
    time = end_time - start_time;
    //cout << "prime number in this process is: " << cnt << endl;
    MPI_Reduce(&cnt, &all, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&time, &total_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    if (rank == 0) {
        cout << "total prime number is: " << all << endl;
        cout << "maximum time is: " << total_time << endl;
    }
    MPI_Finalize();
    
    return 0;
}
