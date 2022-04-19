CFLAGS = -Wall -O3 -c 
LFLAGS = -Wall -O3  -lm 
LLFLAGS = -lpthread -lcplex -L"/opt/ibm/ILOG/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic" -L"/home/kocu/3501g/Cbc-2.9.8/lib" -L"/home/kocu/3501g/OSI-CONIC/lib" -lClp -lOsi -lCgl -lCoinUtils -lClpSolver -lOsiClp -lOsiCpx -lOsiIpopt -lOsiConic -L"/home/kocu/3501g/OsiIpopt/lib" 
GCFLAGS = -g -c 
GLFLAGS = -Wall -g  -lm -lpthread -lcplex -L"/opt/ibm/ILOG/CPLEX_Studio201/cplex/lib/x86-64_linux/static_pic"
INC = -I"/home/kocu/3501g/Cbc-2.9.8/include/coin" -I"/home/kocu/3501g/Cbc-2.9.8/Osi/src/OsiCpx" -I"/opt/ibm/ILOG/CPLEX_Studio201/cplex/include" -I"/opt/ibm/ILOG/CPLEX_Studio201/concert/include" -I"/home/kocu/3501g/Ipopt/src/Interfaces/" -I"/home/kocu/3501g/Ipopt/src/Common/" -I"/home/kocu/3501g/OsiIpopt/src" -I"/home/kocu/3501g/OSI-CONIC/src/" -I"/home/kocu/3501g/Ipopt/src/LinAlg/"
MCC     = mpic++
CC     = g++

pcpfp: my_solution.o my_objective_function.o  myserializablerowcut.o  others.o p3arallel_cpfp.o parallel_cpfp_stuff.o cpfp_executable.o md5.o 
	$(CC) -o pcpfp my_solution.o my_objective_function.o parallel_cpfp_stuff.o myserializablerowcut.o others.o p3arallel_cpfp.o cpfp_executable.o  md5.o $(LFLAGS) $(LLFLAGS) 

.cpp.o:
	$(CC) $(GCFLAGS) $(INC) $<
clean:
	rm pcpfp my_solution.o my_objective_function.o  myserializablerowcut.o others.o p3arallel_cpfp.o parallel_cpfp_stuff.o cpfp_executable.o  md5.o 
 


