// Translate this file with
//
// g++ -O3 assignment-code.cpp -o assignment-code
//
// Run it with
//
// ./demo-code
//
// There should be a result.pvd file that you can open with Paraview.
// Sometimes, Paraview requires to select the representation "Point Gaussian"
// to see something meaningful.
//
// (C) 2018-2020 Tobias Weinzierl

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <math.h>
#include <limits>
#include <iomanip>



#include <cmath>
#include <vector>


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;


#define CACHE_LINE 64
struct alignas(CACHE_LINE) AlignedTriple {
  double _[3];
};


class VectorArray {
  private:
    size_t length;
    AlignedTriple ** data;
  public:
    VectorArray(size_t length){
        this->length = length;
        data = new AlignedTriple*[length];
        for(size_t i = 0; i<length; i++){
          data[i] = new AlignedTriple();
        }
    }

    ~VectorArray(){
      for(size_t i=0; i<length; i++){
        delete data[i];
      }

      delete[] data;
    }

    double & operator()(int x, int y){
      return (*data[x])._[y];
    };
};



/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
VectorArray* x;
VectorArray* x_half;

/**
 * Equivalent to x storing the velocities.
 */
VectorArray* v;
VectorArray* v_half;

/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;
double   halfTimeStepSize = 0.0;

/**
 * Maximum velocity of all particles.
 */
double   maxV;

/**
 * Minimum distance between two elements.
 */
double   minDx;


/**
 * Force experienced by a particle.
 */
VectorArray* force;

bool * potentialCollision;
int * pcMap;



#define X(a,b) (*x)(a,b)
#define V(a,b) (*v)(a,b)
#define FORCE(a,b) (*force)(a,b)

#define X_HALF(a,b) (*x_half)(a,b)
#define V_HALF(a,b) (*v_half)(a,b)



/**
 * Set up scenario from the command line.
 *
 * If you need additional helper data structures, you can
 * initialise them here. Alternatively, you can introduce a
 * totally new function to initialise additional data fields and
 * call this new function from main after setUp(). Either way is
 * fine.
 *
 * This operation's semantics is not to be changed in the assignment.
 */
void setUp(int argc, char** argv) {
  NumberOfBodies = (argc-4) / 7;


  mass = new double [NumberOfBodies];
  potentialCollision = new bool [NumberOfBodies];
  pcMap = new int [NumberOfBodies];

  x = new VectorArray(NumberOfBodies);
  v = new VectorArray(NumberOfBodies);
  force = new VectorArray(NumberOfBodies);

  x_half = new VectorArray(NumberOfBodies);
  v_half = new VectorArray(NumberOfBodies);

  

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;
  halfTimeStepSize = timeStepSize * 0.5;

  for (int i=0; i<NumberOfBodies; i++) {
    //x[i] = new AlignedTriple();
    //v[i] = new AlignedTriple();
    //force[i] = new AlignedTriple();

    X(i, 0) = std::stof(argv[readArgument]); readArgument++;
    X(i, 1) = std::stof(argv[readArgument]); readArgument++;
    X(i, 2) = std::stof(argv[readArgument]); readArgument++;
    
    V(i, 0) = std::stof(argv[readArgument]); readArgument++;
    V(i, 1) = std::stof(argv[readArgument]); readArgument++;
    V(i, 2) = std::stof(argv[readArgument]); readArgument++;

    mass[i] = std::stof(argv[readArgument]); readArgument++;

    if (mass[i]<=0.0 ) {
      std::cerr << "invalid mass for body " << i << std::endl;
      exit(-2);
    }
  }

  std::cout << "created setup with " << NumberOfBodies << " bodies" << std::endl;
  
  if (tPlotDelta<=0.0) {
    std::cout << "plotting switched off" << std::endl;
    tPlot = tFinal + 1.0;
  }
  else {
    std::cout << "plot initial setup plus every " << tPlotDelta << " time units" << std::endl;
    tPlot = 0.0;
  }
}


std::ofstream videoFile;


/**
 * This operation is not to be changed in the assignment.
 */
void openParaviewVideoFile() {
  videoFile.open( "result.pvd" );
  videoFile << "<?xml version=\"1.0\"?>" << std::endl
            << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">" << std::endl
            << "<Collection>";
}





/**
 * This operation is not to be changed in the assignment.
 */
void closeParaviewVideoFile() {
  videoFile << "</Collection>"
            << "</VTKFile>" << std::endl;
}


/**
 * The file format is documented at http://www.vtk.org/wp-content/uploads/2015/04/file-formats.pdf
 *
 * This operation is not to be changed in the assignment.
 */
void printParaviewSnapshot() {
  static int counter = -1;
  counter++;
  std::stringstream filename;
  filename << "result-" << counter <<  ".vtp";
  std::ofstream out( filename.str().c_str() );
  out << "<VTKFile type=\"PolyData\" >" << std::endl
      << "<PolyData>" << std::endl
      << " <Piece NumberOfPoints=\"" << NumberOfBodies << "\">" << std::endl
      << "  <Points>" << std::endl
      << "   <DataArray type=\"Float64\" NumberOfComponents=\"3\" format=\"ascii\">";
//      << "   <DataArray type=\"Float32\" NumberOfComponents=\"3\" format=\"ascii\">";

  for (int i=0; i<NumberOfBodies; i++) {
    out << X(i, 0)
        << " "
        << X(i, 1)
        << " "
        << X(i, 2)
        << " ";
  }

  out << "   </DataArray>" << std::endl
      << "  </Points>" << std::endl
      << " </Piece>" << std::endl
      << "</PolyData>" << std::endl
      << "</VTKFile>"  << std::endl;

  videoFile << "<DataSet timestep=\"" << counter << "\" group=\"\" part=\"0\" file=\"" << filename.str() << "\"/>" << std::endl;
}



#define SQUARED(e) (e)*(e)
#define PRINT_PARTICAL(i) printf("Partical[%d]: x(%f, %f, %f) v(%f, %f, %f) mass: %f\n", i, x[i][0], x[i][1], x[i][2], v[i][0], v[i][1], v[i][2], mass[i])


inline void filterMerge(double& minDxSquared){

  for(int i=0; i<NumberOfBodies; i++){
    potentialCollision[i] = false;
  }

  const double C = 10e-2;
  double loopMinDX = std::numeric_limits<double>::max();

  #pragma omp parallel for reduction(min:loopMinDX)
  for(int i=0; i<NumberOfBodies; i++){
    for (int j=0; j<NumberOfBodies; j++) {
      if(i==j) continue;

      /// Calculate i,j distance
      const double distanceSquared = SQUARED(X(i, 0)-X(j,0)) + SQUARED(X(i, 1)-X(j,1)) + SQUARED(X(i,2)-X(j,2));
      loopMinDX = std::min( loopMinDX, distanceSquared );

      if(distanceSquared<=SQUARED(C*(mass[j]+mass[i]))){
        potentialCollision[i] = true;
        potentialCollision[j] = true;
      }

    }
  }

  minDxSquared = std::min(minDxSquared, loopMinDX);

}

bool isIn(int*list, int&length, int item, int&returnV){
  for(int i=item; i<length; i++){
    if(item==list[i]){
      returnV = i;
      return true;
    }
  }

  return false;
} 

inline void mergeParticales(double & maxVSquared, double &minDxSquared){
  
  int pcCount = 0;
   for(int i=0; i<NumberOfBodies; i++){
     if(potentialCollision[i]){
       pcMap[pcCount] = i;
       pcCount ++;
     }
  }
  
  
  
  // Check and merge
  int i=0;
  const double C = 10e-2;
  while(i<pcCount){
    int j = i+1;
    int s = pcMap[i];
    int t = pcMap[j];
    
    bool merged = false;
    // Check to end of list or until merged
    while( j<pcCount && !merged){
      
      const double distanceSquared = SQUARED(X(s, 0)-X(t,0)) + SQUARED(X(s, 1)-X(t,1)) + SQUARED(X(s,2)-X(t,2));
      minDxSquared = std::min( minDxSquared, distanceSquared );

      if(distanceSquared<=SQUARED(C*(mass[t]+mass[s]))){
        merged = true;
        break;
      }else{
        j++;
      }
    }

    if(merged){

      // Merge particales s and t into s
      V(s, 0) = (mass[s]*V(s, 0)+mass[t]*V(t, 0))/(mass[s]+mass[t]);
      V(s, 1) = (mass[s]*V(s, 1)+mass[t]*V(t, 1))/(mass[s]+mass[t]);
      V(s, 2) = (mass[s]*V(s, 2)+mass[t]*V(t, 2))/(mass[s]+mass[t]);
      
      X(s, 0) = (mass[s]*X(s, 0)+mass[t]*X(t, 0))/(mass[s]+mass[t]);
      X(s, 1) = (mass[s]*X(s, 1)+mass[t]*X(t, 1))/(mass[s]+mass[t]);
      X(s, 2) = (mass[s]*X(s, 2)+mass[t]*X(t, 2))/(mass[s]+mass[t]);

      mass[s] = mass[s]+mass[t];

      maxVSquared = std::max( maxVSquared, ( SQUARED(V(s, 0)) + SQUARED(V(s, 1)) + SQUARED(V(s, 2))));
      
      // Remove t by swapping in the last particle
      int lastBody = NumberOfBodies - 1;
      
      V(t, 0) = V(lastBody, 0);
      V(t, 1) = V(lastBody, 1);
      V(t, 2) = V(lastBody, 2);

      X(t, 0) = X(lastBody, 0);
      X(t, 1) = X(lastBody, 1);
      X(t, 2) = X(lastBody, 2);
      
      mass[t] = mass[lastBody];

      NumberOfBodies--;

      // Remove t from map
      pcMap[j] = pcMap[pcCount-1];
      pcCount--; 

      // If last body is in map update reference
      int matchI;
      if(isIn(pcMap, pcCount, lastBody, matchI)){
        pcMap[matchI] = j;
      }


    }else{
      i++;
    }
  }
}

inline void takeStep(VectorArray& in_x, VectorArray& in_v, VectorArray& out_x, VectorArray& out_v, double timeStep){
    // Zero forces array
  #pragma omp simd
  for(int i=0; i<NumberOfBodies; i++){
    for(int dim =0; dim<3; dim++){
     FORCE(i, dim) = 0.0;
    }
  }

  #pragma omp parallel for
  for(int i=0; i<NumberOfBodies; i++){
    for (int j=0; j<NumberOfBodies; j++) {
      if(i==j) continue;

      /// Calculate i,j distance
      const double distance = sqrt(
        SQUARED(in_x(i, 0)-in_x(j, 0)) +
        SQUARED(in_x(i, 1)-in_x(j ,1)) +
        SQUARED(in_x(i, 2)-in_x(j, 2))
      );

      // Calculate Forces
      const double invarient = (mass[j])/(distance*distance*distance);

      #pragma omp simd aligned(x:CACHE_LINE) aligned(v:CACHE_LINE) aligned(force:CACHE_LINE)
      for(int dim=0; dim<3; dim++){
        double f = (in_x(j, dim)-in_x(i, dim)) * invarient;
        FORCE(i, dim) += f;
      }
      
    }

     // Incremet x and v
    #pragma omp simd aligned(x:CACHE_LINE) aligned(v:CACHE_LINE) aligned(force:CACHE_LINE)
    for(int dim=0; dim<3; dim++){
      out_x(i, dim) = X(i, dim) + timeStep * in_v(i, dim);
      out_v(i, dim) = V(i, dim) + timeStep * FORCE(i, dim);
    }  

  }

}

inline void setMaxV(double & maxVSquared){
  for(int i =0; i<NumberOfBodies; i++){
    maxVSquared = std::max( maxVSquared, ( SQUARED(V(i, 0)) + SQUARED(V(i, 1)) + SQUARED(V(i, 2))));
  }
}

void printVectorArray(VectorArray* array){
  for(int i=0; i<NumberOfBodies; i++){
    std::cout << "(" << (*array)(i, 0) << ", " << (*array)(i, 1) << ", " << (*array)(i, 2) << ")";
  }
  std::cout << std::endl;
}

/**
 * This is the main operation you should change in the assignment. You might
 * want to add a few more variables or helper functions, but this is where the
 * magic happens.
 */
void updateBody() {
  double maxVSquared   = 0.0;
  double minDxSquared  = std::numeric_limits<double>::max();

  //printVectorArray(x);
  //printVectorArray(v);
  takeStep(*x, *v, *x_half, *v_half, halfTimeStepSize);
  //printVectorArray(x_half);
  //printVectorArray(v_half);
  takeStep(*x_half, *v_half, *x, *v, timeStepSize);
  setMaxV(maxVSquared);
  filterMerge(minDxSquared);
  mergeParticales(maxVSquared, minDxSquared);
  

  t += timeStepSize;
  maxV = std::sqrt(maxVSquared);
  minDx = minDxSquared==std::numeric_limits<double>::max()?std::numeric_limits<double>::max():std::sqrt(minDxSquared);

}


/**
 * Main routine.
 *
 * No major changes in assignment. You can add a few initialisation
 * or stuff if you feel the need to do so. But keep in mind that you
 * may not alter what the program plots to the terminal.
 */
int main(int argc, char** argv) {
  if (argc==1) {
    std::cerr << "usage: " + std::string(argv[0]) + " snapshot final-time dt objects" << std::endl
              << "  snapshot        interval after how many time units to plot. Use 0 to switch off plotting" << std::endl
              << "  final-time      simulated time (greater 0)" << std::endl
              << "  dt              time step size (greater 0)" << std::endl
              << std::endl
              << "Examples:" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0 \t One body moving form the coordinate system's centre along x axis with speed 1" << std::endl
              << "0.01  100.0  0.001    0.0 0.0 0.0  1.0 0.0 0.0  1.0     0.0 1.0 0.0  1.0 0.0 0.0  1.0  \t One spiralling around the other one" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0 \t Three body setup from first lecture" << std::endl
              << "0.01  100.0  0.001    3.0 0.0 0.0  0.0 1.0 0.0  0.4     0.0 0.0 0.0  0.0 0.0 0.0  0.2     2.0 0.0 0.0  0.0 0.0 0.0  1.0     2.0 1.0 0.0  0.0 0.0 0.0  1.0     2.0 0.0 1.0  0.0 0.0 0.0  1.0 \t Five body setup" << std::endl
              << std::endl
              << "In this naive code, only the first body moves" << std::endl;

    return -1;
  }
  else if ( (argc-4)%7!=0 ) {
    std::cerr << "error in arguments: each planet is given by seven entries (position, velocity, mass)" << std::endl;
    std::cerr << "got " << argc << " arguments (three of them are reserved)" << std::endl;
    std::cerr << "run without arguments for usage instruction" << std::endl;
    return -2;
  }

  std::cout << std::setprecision(15);

  setUp(argc,argv);

  openParaviewVideoFile();

  int snapshotCounter = 0;
  if (t > tPlot) {
    printParaviewSnapshot();
    std::cout << "plotted initial setup" << std::endl;
    tPlot = tPlotDelta;
  }

  int timeStepCounter = 0;
  while (t<=tFinal) {
    updateBody();
    timeStepCounter++;
    if (t >= tPlot) {
      printParaviewSnapshot();
      std::cout << "plot next snapshot"
    		    << ",\t time step=" << timeStepCounter
    		    << ",\t t="         << t
				<< ",\t dt="        << timeStepSize
				<< ",\t v_max="     << maxV
				<< ",\t dx_min="    << minDx
				<< std::endl;

      tPlot += tPlotDelta;
    }
  }

  std::cout << "Number of remaining objects: " << NumberOfBodies << std::endl;
  std::cout << "Position of first remaining object: " << X(0, 0) << ", " << X(0, 1) << ", " << X(0, 2) << std::endl;

  closeParaviewVideoFile();

  return 0;
}
