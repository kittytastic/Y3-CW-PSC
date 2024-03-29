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


double t          = 0;
double tFinal     = 0;
double tPlot      = 0;
double tPlotDelta = 0;

int NumberOfBodies = 0;


#define CACHE_LINE 64
struct alignas(CACHE_LINE) AlignedTriple {
  double _[3];
};

// Future me: Turns out an array of aligned structs could be an anti-patten and not creeate the desired alignment
// This data-structure was reworked in step-4 
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
/**
 * Equivalent to x storing the velocities.
 */
VectorArray* v;
/**
 * One mass entry per molecule/particle.
 */
double*  mass;

/**
 * Global time step size used.
 */
double   timeStepSize = 0.0;

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

double C;

#define X(a,b) (*x)(a,b)
#define V(a,b) (*v)(a,b)
#define FORCE(a,b) (*force)(a,b)


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
  C = 10e-2/NumberOfBodies;

  mass = new double [NumberOfBodies];
  x = new VectorArray(NumberOfBodies);
  v = new VectorArray(NumberOfBodies);
  force = new VectorArray(NumberOfBodies);

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {

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


/**
 * This is the main operation you should change in the assignment. You might
 * want to add a few more variables or helper functions, but this is where the
 * magic happens.
 */
void updateBody() {
  double maxVSquared   = 0.0;
  double minDxSquared  = std::numeric_limits<double>::max();

  // Zero forces array
  #pragma omp simd
  for(int i=0; i<NumberOfBodies; i++){
    for(int dim =0; dim<3; dim++){
     FORCE(i, dim) = 0.0;
    }
  }

  // Euler Step
  for(int i=0; i<NumberOfBodies; i++){
    for (int j=i+1; j<NumberOfBodies; j++) {

      /// Calculate i,j distance
      // SIMD reductions appear to be a detriment
      const double distance = sqrt(
        SQUARED(X(i, 0)-X(j, 0)) +
        SQUARED(X(i, 1)-X(j ,1)) +
        SQUARED(X(i, 2)-X(j, 2))
      );

      // Calculate Forces
      const double invarient = (mass[j]*mass[i])/(distance*distance*distance);

      #pragma omp simd aligned(x:CACHE_LINE) aligned(v:CACHE_LINE) aligned(force:CACHE_LINE)
      for(int dim=0; dim<3; dim++){
        double f = (X(j, dim)-X(i, dim)) * invarient;
        FORCE(i, dim) += f;
        FORCE(j, dim) -= f;
      }
      
    }

    // Incremet x and v
    #pragma omp simd aligned(x:CACHE_LINE) aligned(v:CACHE_LINE) aligned(force:CACHE_LINE)
    for(int dim=0; dim<3; dim++){
      X(i, dim) += timeStepSize * V(i,dim);
      V(i, dim) += timeStepSize * FORCE(i, dim) / mass[i];
    }
  }

  int i=0;
  
  while(i<NumberOfBodies){
    int j = i+1;
    bool merged = false;
    while( j<NumberOfBodies && !merged){
      
      // SIMD reductions appear to be a detriment
      const double distanceSquared = SQUARED(X(i, 0)-X(j,0)) + SQUARED(X(i, 1)-X(j,1)) + SQUARED(X(i,2)-X(j,2));
      
      if(distanceSquared<=SQUARED(C*(mass[j]+mass[i]))){
        merged = true;
        break;
      }else{
        minDxSquared = std::min( minDxSquared, distanceSquared );
        j++;
      }
    }

    if(merged){
 
      // Merge i and j into i
      // Very cold code, no point SIMD, undecided if this or SIMD looks better...
      V(i, 0) = (mass[i]*V(i, 0)+mass[j]*V(j, 0))/(mass[i]+mass[j]);
      V(i, 1) = (mass[i]*V(i, 1)+mass[j]*V(j, 1))/(mass[i]+mass[j]);
      V(i, 2) = (mass[i]*V(i, 2)+mass[j]*V(j, 2))/(mass[i]+mass[j]);
      
      X(i, 0) = (mass[i]*X(i, 0)+mass[j]*X(j, 0))/(mass[i]+mass[j]);
      X(i, 1) = (mass[i]*X(i, 1)+mass[j]*X(j, 1))/(mass[i]+mass[j]);
      X(i, 2) = (mass[i]*X(i, 2)+mass[j]*X(j, 2))/(mass[i]+mass[j]);

      mass[i] = mass[i]+mass[j];

      // Move last body into j and decrement body count
      int lastBody = NumberOfBodies - 1;
      V(j, 0) = V(lastBody, 0);
      V(j, 1) = V(lastBody, 1);
      V(j, 2) = V(lastBody, 2);

      X(j, 0) = X(lastBody, 0);
      X(j, 1) = X(lastBody, 1);
      X(j, 2) = X(lastBody, 2);
      
      mass[j] = mass[lastBody];

      NumberOfBodies--;

    }else{
      i++;
    }
  }
  
  for(int i=0; i<NumberOfBodies; i++){
    maxVSquared = std::max( maxVSquared, ( SQUARED(V(i, 0)) + SQUARED(V(i, 1)) + SQUARED(V(i, 2))));
  }

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
