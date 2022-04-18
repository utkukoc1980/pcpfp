
// #include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>

#include "mpi.h"

// #include "../templates.h"


using namespace std;

int main(int argc,char* argv[]){
// 	cout << __LINE__ << endl;
// 	cout << "zikkim" << endl;
// 	cout.flush();
	int myid = 1, numprocs= 1;
	int source,count;
	MPI_Status status;
	MPI_Request request =MPI_REQUEST_NULL ;

	MPI::Init_thread(argc,argv,MPI_THREAD_MULTIPLE);
	myid = MPI::COMM_WORLD.Get_rank();
	numprocs = MPI::COMM_WORLD.Get_size();
	source=0;
	count=4;
// 	cout << "numprocs " << numprocs<< endl;
	int TAG = 123;
	int sleeptime = 9800000;
	if(myid == source){
		cout << "processor " << myid << " will sleep for " << sleeptime/(1000000.0) << " seconds " << endl;
		usleep(sleeptime);
		int d = 1;
		
		for(int i=1;i<numprocs;i++){
			cout << "sending data " <<endl;
			MPI_Send(&d,1,MPI_INT,i,TAG+1,MPI_COMM_WORLD);
			cout << "sending completed" << endl;
		}
	}
	else {
		bool recv = false;
		int d2 = 0;
		while (!recv){
			usleep(1000000);
// 			cout << " Irecv:" << 
// 			<<endl;
// 			cout.flush();
			cout << "d2 = " << d2 << endl;
			int flag = 99;
			int t =MPI_Iprobe(0,TAG,MPI_COMM_WORLD, &flag, &status);
			cout << "flag: " << flag << endl;
            cout << "TAG m : " << TAG << endl;
			if (flag==1){
				cout << "__ " << __LINE__ << endl;
				MPI_Recv(&d2, 1, MPI_INT,0,TAG,MPI_COMM_WORLD,&status);
                cout << "TAG s : " << TAG << endl;
                
			}
			if(d2 == 1)
				recv = true;
			
			
		}
		
	}
	
	cout << __LINE__ << endl;
	cout << "process " << myid << " terminates" << endl;
		
	MPI_Finalize();
	return 0;
}