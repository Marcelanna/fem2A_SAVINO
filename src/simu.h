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

	double sinus_bump( vertex v )
	{
		return 2 * pow(3.1415, 2) * sin(3.1415* v.x) * sin(3.1415 * v.y);
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
        
        
        
        
        void source_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source term Dirichlet problem" << std::endl;           
            if ( verbose ) 
            {
                std::cout << " run simulation " << std::endl << "A null limit Dirichlet condition  is imposed to all the edges of the mesh with a source term equal to one." << std::endl;
            }
            
            Mesh mesh;
            mesh.load(mesh_filename);
            ShapeFunctions reference_functions  = ShapeFunctions( 2, 1 );
            Quadrature quadrature = Quadrature::get_quadrature(2, false);
            SparseMatrix K(mesh.nb_vertices());  
            
            
            //K initialization          
            for(int t=0; t < mesh.nb_triangles(); t++)
            {
            	ElementMapping elt_mapping = ElementMapping( mesh, false, t );
    		DenseMatrix Ke;
    		assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, unit_fct, Ke);

    		local_to_global_matrix( mesh, t, Ke, K );
            }
            
            //attribut_is_dirichlet initialization
            std::vector< bool > attribute_is_dirichlet(1, true);
            /*std::vector< bool > attribute_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;*/
            	
            mesh.set_attribute(unit_fct, 1, true );

            //Value initialization
            std::vector< double > values;
            for( int i=0; i<mesh.nb_vertices(); i++)
            {
            	values.push_back(zero_fct(mesh.get_vertex(i)));
            }
            
            //F and Fe initialization
            std::vector< double > Fe(reference_functions.nb_functions(),0 );
            std::vector< double > F(mesh.nb_vertices(),0);
            
            for(int t=0; t < mesh.nb_triangles(); t++)
            {
            	ElementMapping elt_mapping = ElementMapping( mesh, false, t );
    		assemble_elementary_vector(elt_mapping, reference_functions, quadrature, unit_fct, Fe);

    		local_to_global_vector(mesh, false, t, Fe, F );
            }
            
            apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F );
            
            std::vector<double> U(mesh.nb_vertices());
            solve(K, F, U);
            std::string export_name = "source_dirichlet";
            mesh.save(export_name + ".mesh");
            save_solution(U, export_name + ".bb");
            
            
        }
        
        
        
        


	void sinus_bump_dirichlet_pb( const std::string& mesh_filename, bool verbose )
        {
            std::cout << "Solving a source term Dirichlet problem" << std::endl;           
            if ( verbose ) 
            {
                std::cout << " run simulation " << std::endl << "A null limit Dirichlet condition  is imposed to all the edges of the mesh with a source term equal to one." << std::endl;
            }
            
            Mesh mesh;
            mesh.load(mesh_filename);
            ShapeFunctions reference_functions  = ShapeFunctions( 2, 1 );
            Quadrature quadrature = Quadrature::get_quadrature(2, false);
            SparseMatrix K(mesh.nb_vertices());  
            
            
            //K initialization          
            for(int t=0; t < mesh.nb_triangles(); t++)
            {
            	ElementMapping elt_mapping = ElementMapping( mesh, false, t );
    		DenseMatrix Ke;
    		assemble_elementary_matrix(elt_mapping, reference_functions, quadrature, unit_fct, Ke);

    		local_to_global_matrix( mesh, t, Ke, K );
            }
            
            //attribut_is_dirichlet initialization
            std::vector< bool > attribute_is_dirichlet(1, true);
            /*std::vector< bool > attribute_is_dirichlet(2, false);
            att_is_dirichlet[1] = true;*/
            	
            mesh.set_attribute(unit_fct, 1, true );

            //Value initialization
            std::vector< double > values;
            for( int i=0; i<mesh.nb_vertices(); i++)
            {
            	values.push_back(zero_fct(mesh.get_vertex(i)));
            }
            
            //F and Fe initialization
            std::vector< double > Fe(reference_functions.nb_functions(),0 );
            std::vector< double > F(mesh.nb_vertices(),0);
            
            for(int t=0; t < mesh.nb_triangles(); t++)
            {
            	ElementMapping elt_mapping = ElementMapping( mesh, false, t );
    		assemble_elementary_vector(elt_mapping, reference_functions, quadrature, sinus_bump, Fe);

    		local_to_global_vector(mesh, false, t, Fe, F );
            }
            
            apply_dirichlet_boundary_conditions(mesh, attribute_is_dirichlet, values, K, F );
            
            std::vector<double> U(mesh.nb_vertices());
            solve(K, F, U);
            std::string export_name = "sin_bump_dirichlet";
            mesh.save(export_name + ".mesh");
            save_solution(U, export_name + ".bb");
            
            
        }
    }

}
