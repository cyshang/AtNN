objects = main.o Parameter.o NeuralNetwork.o NetworkInfo.o Group.o SymFunction.o FuncGroup.o FuncInfo.o Molecule.o 
EIGEN = -I ../Tools/Eigen3/
LIB = -L /work1/soft/intel2015/composer_xe_2015.2.164/mkl/lib/intel64 -lmkl_rt
MKL = -I /work1/soft/intel2015/mkl/include
options = -O2

TARGET = AtNN.exe

$(TARGET) : $(objects)
	icpc -o $(TARGET) $(objects) -lm $(LIB) $(options)

main.o : main.cpp global.h NeuralNetwork.h SymFunction.h Molecule.h
	icpc $(EIGEN) $(MKL) -c main.cpp $(options)

Parameter.o : Parameter.cpp Parameter.h
	icpc $(EIGEN) $(MKL) -c Parameter.cpp $(options)	
	
NeuralNetwork.o : NeuralNetwork.cpp NeuralNetwork.h SymFunction.h NetworkInfo.h Group.h global.h
	icpc $(EIGEN) $(MKL) -c NeuralNetwork.cpp $(options)

NetworkInfo.o : NetworkInfo.cpp NetworkInfo.h SymFunction.h FuncGroup.h global.h
	icpc $(EIGEN) $(MKL) -c NetworkInfo.cpp $(options)

Group.o : Group.cpp Group.h NeuralNetwork.h NetworkInfo.h global.h
	icpc $(EIGEN) $(MKL) -c Group.cpp $(options)

SymFunction.o : SymFunction.cpp SymFunction.h FuncGroup.h FuncInfo.h Molecule.h NeuralNetwork.h global.h
	icpc $(EIGEN) $(MKL) -c SymFunction.cpp $(options)

FuncGroup.o : FuncGroup.cpp FuncGroup.h FuncInfo.h global.h
	icpc $(EIGEN) $(MKL) -c FuncGroup.cpp $(options)

FuncInfo.o : FuncInfo.cpp FuncInfo.h global.h
	icpc $(EIGEN) $(MKL) -c FuncInfo.cpp $(options)

Molecule.o : Molecule.cpp Molecule.h SymFunction.h FuncGroup.h FuncInfo.h global.h
	icpc $(EIGEN) $(MKL) -c Molecule.cpp $(options)
	
	
.PHONY : clean
clean :
	-rm $(TARGET) $(objects)
