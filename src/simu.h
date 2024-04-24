#pragma once

#include "mesh.h"
#include "fem.h"
#include <math.h>
#include <cmath>
#include <iostream>

namespace FEM2A {
    namespace Simu {

        //#################################
        //  Useful functions
        //#################################

        double unit_fct( vertex v )
        {
            return 1.;
        }

        double zero_fct( vertex v )
        {
            return 0.;
        }

        double xy_fct( vertex v )
        {
            return v.x + v.y;
        }

        //#################################
        //  Simulations
        //#################################

        void pure_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a pure Dirichlet problem" << std::endl;           
            if ( verbose ) 
            {
                std::cout << " run simulation " << std::endl << "The limit conditions are imposed to all the edges of the mesh following the Dirichlet condition." << std::endl;
            }
            
            Mesh mesh;
            mesh.load(mesh_filename);
            SparseMatrix K(mesh.nb_vertices());  
                      
            for(int t=0; t < mesh.nb_triangles(); t++)
            {
            	ElementMapping elt_mapping = ElementMapping( mesh, false, t );
            	ShapeFunctions reference_functions  = ShapeFunctions( 2, 1 );
    		Quadrature quadrature = Quadrature::get_quadrature(2, false);
    		DenseMatrix Ke;
    		assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, unit_fct, Ke);

    		local_to_global_matrix( mesh, t, Ke, K );
    		
            	
            }
            std::vector< bool > attribute_is_dirichlet(1, true);
            /*std::vector< bool > attribute_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;*/
            	
            mesh.set_attribute(unit_fct, 1, true );

            std::vector< double > values;
            for( int i=0; i<mesh.nb_vertices(); i++)
            {
            	values.push_back(mesh.get_vertex(i).x + mesh.get_vertex(i).y);
            }
            std::vector< double > F(mesh.nb_vertices(), 0);
            
            apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F );
            
            std::vector<double> U(mesh.nb_vertices());
            solve(K, F, U);
            std::string export_name = "pure_dirichlet";
            mesh.save(export_name + ".mesh");
            save_solution(U, export_name + ".bb");
            
            
        }

    }

}
