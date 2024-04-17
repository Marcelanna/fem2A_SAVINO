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
{
    const bool t_opennl = false;
    const bool t_lmesh = false;
    const bool t_io = false;
    const bool t_quadrature = false;
    
    if( t_opennl ) test_opennl();
    if( t_lmesh ) Tests::test_load_mesh();
    if( t_io ) Tests::test_load_save_mesh();
    if( t_quadrature ) Tests::test_quadrature(4, false);
    
    //Tests ElementMapping
    const bool t_elmt_map = false;
    const bool t_transf = false;
    const bool t_jacob = false;

    if( t_elmt_map ) Tests::test_ElementMapping("data/square.mesh", false, 4);
    if( t_transf ) Tests::test_transf(0.2, 0.4);
    if( t_jacob ) Tests::test_jacob(0.2, 0.4);
    
    //Tests ShapeFunction
    const bool t_shpfunct = false;
    const bool t_evaluate = false;
    
    if( t_shpfunct ) Tests::test_shpfct(2,1);
    if( t_evaluate ) Tests::test_evaluate(2, 2, 0.2, 0.4);
    
    //Test Finite Element functions
    const bool t_ass_el_matrix = false;
    const bool t_loc_to_mat = false;
    const bool t_dirichlet_boundary_cdt = false;
    
    if( t_ass_el_matrix ) Tests::test_ass_el_matrix("data/square.mesh", false, 2, 2, 4);
    if( t_loc_to_mat ) Tests::test_local_to_global("data/square.mesh", 4);
    if( t_dirichlet_boundary_cdt ) Tests::test_dirichlet_boundary_cdt("data/square.mesh");
}

void run_simu()
{

    const bool simu_pure_dirichlet = true;

    const bool verbose = flag_is_used( "-v", arguments )
        || flag_is_used( "--verbose", arguments );

    if( simu_pure_dirichlet ) {
        Simu::pure_dirichlet_pb("data/square.mesh", verbose);
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
