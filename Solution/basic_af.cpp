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

/**
 * Pointer to pointers. Each pointer in turn points to three coordinates, i.e.
 * each pointer represents one molecule/particle/body.
 */
double** x;

/**
 * Equivalent to x storing the velocities.
 */
double** v;

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

  x    = new double*[NumberOfBodies];
  v    = new double*[NumberOfBodies];
  mass = new double [NumberOfBodies];

  int readArgument = 1;

  tPlotDelta   = std::stof(argv[readArgument]); readArgument++;
  tFinal       = std::stof(argv[readArgument]); readArgument++;
  timeStepSize = std::stof(argv[readArgument]); readArgument++;

  for (int i=0; i<NumberOfBodies; i++) {
    x[i] = new double[3];
    v[i] = new double[3];

    x[i][0] = std::stof(argv[readArgument]); readArgument++;
    x[i][1] = std::stof(argv[readArgument]); readArgument++;
    x[i][2] = std::stof(argv[readArgument]); readArgument++;

    v[i][0] = std::stof(argv[readArgument]); readArgument++;
    v[i][1] = std::stof(argv[readArgument]); readArgument++;
    v[i][2] = std::stof(argv[readArgument]); readArgument++;

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
    out << x[i][0]
        << " "
        << x[i][1]
        << " "
        << x[i][2]
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
  maxV   = 0.0;
  minDx  = std::numeric_limits<double>::max();

  // force0 = force along x direction
  // force1 = force along y direction
  // force2 = force along z direction
  double* force0 = new double[NumberOfBodies];
  double* force1 = new double[NumberOfBodies];
  double* force2 = new double[NumberOfBodies];

  for(int i=0; i<NumberOfBodies; i++){

    //PRINT_PARTICAL(i);
    force0[i] = 0.0;
    force1[i] = 0.0;
    force2[i] = 0.0;
    for (int j=0; j<NumberOfBodies; j++) {
      if(i==j){ continue;}

      const double distance = sqrt(
        SQUARED(x[i][0]-x[j][0]) +
        SQUARED(x[i][1]-x[j][1]) +
        SQUARED(x[i][2]-x[j][2])
      );

      // x,y,z forces acting on particle 0
      force0[i] += (x[j][0]-x[i][0]) * mass[j]*mass[i] / distance / distance / distance ;
      force1[i] += (x[j][1]-x[i][1]) * mass[j]*mass[i] / distance / distance / distance ;
      force2[i] += (x[j][2]-x[i][2]) * mass[j]*mass[i] / distance / distance / distance ;

      
    }

    x[i][0] = x[i][0] + timeStepSize * v[i][0];
    x[i][1] = x[i][1] + timeStepSize * v[i][1];
    x[i][2] = x[i][2] + timeStepSize * v[i][2];

    v[i][0] = v[i][0] + timeStepSize * force0[i] / mass[i];
    v[i][1] = v[i][1] + timeStepSize * force1[i] / mass[i];
    v[i][2] = v[i][2] + timeStepSize * force2[i] / mass[i];

    //PRINT_PARTICAL(i);
    maxV = std::max( maxV, std::sqrt( SQUARED(v[i][0]) + SQUARED(v[i][1]) + SQUARED(v[i][2]) ));
  }

  int i=0;
  double C = 10e-2; 
  while(i<NumberOfBodies){
    int j = i+1;
    bool merged = false;
    while( j<NumberOfBodies && !merged){
       const double distance = sqrt(
        SQUARED(x[i][0]-x[j][0]) +
        SQUARED(x[i][1]-x[j][1]) +
        SQUARED(x[i][2]-x[j][2])
      );

      minDx = std::min( minDx, distance );

      //printf("i: %d  j: %d   dx: %f\n", i, j, distance);
      if(distance<=C*(mass[j]+mass[i])){
        merged = true;
        break;
      }else{
        j++;
      }
    }

    if(merged){
      std::cout << "Merging 2 particals" <<std::endl;
      PRINT_PARTICAL(i);
      PRINT_PARTICAL(j);

      // Merge i and j into i
      v[i][0] = (mass[i]*v[i][0]+mass[j]*v[j][0])/(mass[i]+mass[j]);
      v[i][1] = (mass[i]*v[i][1]+mass[j]*v[j][1])/(mass[i]+mass[j]);
      v[i][2] = (mass[i]*v[i][2]+mass[j]*v[j][2])/(mass[i]+mass[j]);
      
      x[i][0] = (mass[i]*x[i][0]+mass[j]*x[j][0])/(mass[i]+mass[j]);
      x[i][1] = (mass[i]*x[i][1]+mass[j]*x[j][1])/(mass[i]+mass[j]);
      x[i][2] = (mass[i]*x[i][2]+mass[j]*x[j][2])/(mass[i]+mass[j]);

      mass[i] = mass[i]+mass[j];

      maxV = std::max( maxV, std::sqrt( SQUARED(v[i][0]) + SQUARED(v[i][1]) + SQUARED(v[i][2]) ));

      //printf("Merged Partical\n");
      //PRINT_PARTICAL(i);

      // Move last body into j and decrement body count
      int lastBody = NumberOfBodies - 1;
      v[j][0] = v[lastBody][0];
      v[j][1] = v[lastBody][1];
      v[j][2] = v[lastBody][2];
      
      x[j][0] = x[lastBody][0];
      x[j][1] = x[lastBody][1];
      x[j][2] = x[lastBody][2];
      
      mass[j] = mass[lastBody];

      NumberOfBodies--;

    }else{
      i++;
    }


  }
  

  t += timeStepSize;

  delete[] force0;
  delete[] force1;
  delete[] force2;
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
  std::cout << "Position of first remaining object: " << x[0][0] << ", " << x[0][1] << ", " << x[0][2] << std::endl;

  closeParaviewVideoFile();

  return 0;
}
