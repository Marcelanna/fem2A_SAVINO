#include <iostream>
#include <string>
#include <vector>

#include "src/fem.h"
#include "src/mesh.h"
#include "src/solver.h"
#include "src/tests.h"
#include "src/simu.h"

/* Global variables */
std::vector< std::string > arguments;

/* To parse command line arguments */
bool flag_is_used(
    const std::string& flag,
    const std::vector< std::string >& arguments )
{
    for( int i = 0; i < arguments.size(); ++i ) {
        if( flag == arguments[i] ) {
            return true;
        }
    }
    return false;
}

using namespace FEM2A;

void run_tests()
//to test the desired function, give the value "true" to the associated bool//

{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quadrature = false;
    
    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    //change the varible depending on the element tested
    if( t_quadrature ) Tests::test_quadrature(4, false);
    
    
    
    ////Tests ElementMapping////
    const bool t_elmt_map = false;
    const bool t_transf = false;
    const bool t_jacob = false;

    //change the mesh name, the element tested (true = edge, false = triangle) and its index
    if( t_elmt_map ) Tests::test_ElementMapping("data/square.mesh", false, 4);
    //change the coordinates depending the studied point
    if( t_transf ) Tests::test_transf(0.2, 0.4);
    //change the coordinates depending the studied point
    if( t_jacob ) Tests::test_jacob(0.2, 0.4);
    
    
    
    ////Tests ShapeFunction////
    const bool t_shpfunct = false;
    const bool t_evaluate = false;
    
    //change the dimension and the order depending on the studied case
    if( t_shpfunct ) Tests::test_shpfct(2,1);
    //change the dimension, the order and the coordinates depending on the evaluation wanted
    if( t_evaluate ) Tests::test_evaluate(2, 2, 0.2, 0.4);
    
    
    
    ////Test Finite Element functions////
    const bool t_ass_el_matrix = false;
    const bool t_loc_to_mat = false;
    const bool t_dirichlet_boundary_cdt = false;
    const bool t_vector_Fe = false;
    const bool t_loc_to_vect = false;
    
    //change the mesh name, the element tested (true = edge, false = triangle), the order and the dimensionof the fucntion and the element index
    if( t_ass_el_matrix ) Tests::test_ass_el_matrix("data/square.mesh", false, 2, 2, 4);
    //change the mesh name and the triangle index
    if( t_loc_to_mat ) Tests::test_local_to_global_m("data/square.mesh", 4);
    //cgange the mesh name
    if( t_dirichlet_boundary_cdt ) Tests::test_dirichlet_boundary_cdt("data/square.mesh");
    //change the mesh name, the element tested (true = edge, false = triangle), the order and the dimensionof the fucntion and the element index
    if( t_vector_Fe ) Tests::test_ass_el_vector("data/square.mesh", true, 2, 1, 4);
    //change the mesh name, the element tested (true = edge, false = triangle), the order and the dimensionof the fucntion and the element index
    if (t_loc_to_vect) Tests::test_local_to_global_v("data/square.mesh", false, 2, 2, 4);
}

void run_simu()
//to test the desired function, give the value "true" to the associated bool//

{
    const bool simu_pure_dirichlet = false;
    const bool simu_source_dirichlet = false;
    const bool simu_sin_bump_dirichlet = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    //change in the mesh used for the simulation in the calling up of the function//
    
    if( simu_pure_dirichlet ) {
        Simu::pure_dirichlet_pb("data/square_fine.mesh", verbose );
    }
    
    if( simu_source_dirichlet ) {
        Simu::source_dirichlet_pb("data/square_fine.mesh", verbose );
    }
    
    if( simu_sin_bump_dirichlet ) {
        Simu::sinus_bump_dirichlet_pb("data/square_fine.mesh", verbose );
    }
}

int main( int argc, const char * argv[] )
{
    /* Command line parsing */
    for( int i = 1; i < argc; ++i ) {
        arguments.push_back( std::string(argv[i]) );
    }

    /* Show usage if asked or no arguments */
    if( arguments.size() == 0 || flag_is_used("-h", arguments)
        || flag_is_used("--help", arguments) ) {
        std::cout << "Usage: ./fem2a [options]" << std::endl
            << "Options: " << std::endl;
        std::cout << " -h, --help:        show usage" << std::endl;
        std::cout << " -t, --run-tests:   run the tests" << std::endl;
        std::cout << " -s, --run-simu:    run the simulations" << std::endl;
        std::cout << " -v, --verbose:     print lots of details" << std::endl;
        return 0;
    }

    /* Run the tests if asked */
    if( flag_is_used("-t", arguments)
        || flag_is_used("--run-tests", arguments) ) {
        run_tests();
    }

    /* Run the simulation if asked */
    if( flag_is_used("-s", arguments)
        || flag_is_used("--run-simu", arguments) ) {
        run_simu();
    }

    return 0;
}
