#pragma once

#include "mesh.h"
#include "fem.h"
#include "solver.h"

#include <assert.h>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <stdlib.h>

namespace FEM2A {
    namespace Tests {

        bool test_load_mesh()
        {
            Mesh mesh;
            mesh.load("data/square.mesh");

            std::cout << "Vertices <x> <y> <att>" << std::endl;
            for( int v = 0; v < mesh.nb_vertices(); v++ ) {
                std::cout << mesh.get_vertex(v).x << " " << mesh.get_vertex(v).y
                    << " " << mesh.get_vertex_attribute(v) << std::endl;
            }

            std::cout << "Edges <v0> <v1> <att>" << std::endl ;
            for( int ed = 0; ed < mesh.nb_edges(); ed++ ) {
                std::cout << mesh.get_edge_vertex_index(ed, 0) << " "
                    << mesh.get_edge_vertex_index(ed, 1) << " "
                    << mesh.get_edge_attribute(ed) << std::endl;
            }

            std::cout << "Triangles <v0> <v1> <v2> <att>" << std::endl ;
            for( int tr = 0; tr < mesh.nb_triangles(); tr++ ) {
                std::cout << mesh.get_triangle_vertex_index(tr, 0) << " "
                    << mesh.get_triangle_vertex_index(tr, 1) << " "
                    << mesh.get_triangle_vertex_index(tr, 2) << " "
                    << mesh.get_triangle_attribute(tr) << std::endl;
            }

            return true;
        }

        bool test_load_save_mesh()
        {
            Mesh mesh;
            mesh.load("data/geothermie_4.mesh");
            mesh.save("data/geothermie_4.mesh");
            return true;
        }
        
        bool test_quadrature(int order, bool border)
        {
        	Quadrature test;
        	test = test.get_quadrature(order, border);
        	double sum = 0;
        	std::cout << "Test Quadrature" << std::endl;
        	for (int i = 0; i<test.nb_points(); i++){
        		sum += test.weight(i);
        	}
        	std::cout<< sum << std::endl;
        	return true;
        }
        
        bool test_ElementMapping(std::string M, bool border, int i)
        {
        	Mesh mesh;
            	mesh.load(M);
            	
        	ElementMapping elmt = ElementMapping( mesh, border, i );
        	
        	return true;
    	}
    	
    	bool test_transf(double x, double y)
    	{
    		Mesh mesh;
            	mesh.load("data/square.mesh");
            	ElementMapping elmt = ElementMapping( mesh, false, 4 );
    		
    		vertex point;
    		point.x = x;
    		point.y = y;
    		
    		vertex real_point;
    		real_point = elmt.transform( point );
    		std::cout << real_point.x << " " << real_point.y << std::endl;
    		
    		return true;
    	}
    	
    	bool test_jacob(double x, double y)
    	{
    		Mesh mesh;
            	mesh.load("data/square.mesh");
            	ElementMapping elmt = ElementMapping( mesh, false, 4 );
    		
    		vertex point;
    		point.x = x;
    		point.y = y;
    		
    		//Test du calcul de la matrice jacobienne
    		std::cout << "Test de la matrice jacobienne"<< std::endl;
    		DenseMatrix jacob_point;
    		jacob_point = elmt.jacobian_matrix( point );
    		jacob_point.print();
    		
    		//Test du calcul du déterminant de la matrice jacobienne
    		std::cout << "Test determinant matrice jacobienne"<< std::endl;
    		double det;
    		det = elmt.jacobian( point );
    		std::cout << det << std::endl;
    		
    		
    		return true;
    	}
    	
    	bool test_shpfct(int dim, int order)
    	{
    		ShapeFunctions elmt = ShapeFunctions( dim, order );
    		
    		std::cout << elmt.nb_functions() << std::endl;
    		
    		return true;
    	}
    	
    	bool test_evaluate(int dim, int i, double x, double y)
    	{
    		ShapeFunctions elmt = ShapeFunctions( dim, 1 );
    		
    		
    		vertex point;
    		point.x = x;
    		point.y = y;
    		
    		//Test of the evaluation of the i-th shapefunction
    		std::cout << "Function evaluation" << elmt.evaluate(i, point) << std::endl;
    		
    		//Test of the evaluation of the gradiant of the i-th shapefunction
		std::cout << "Gradient evaluation" << elmt.evaluate_grad(i, point).x << elmt.evaluate_grad(i, point).y << std::endl;
    		
    		return true;
    	}
    	
    	double unit_fct( vertex v )
        {
            return 1.;
        }
    	
    	bool test_ass_el_matrix(std::string M, bool border, int order, int dim, int i)
    	{
    		Mesh mesh;
            	mesh.load(M);
            	
        	ElementMapping elt_mapping = ElementMapping( mesh, border, i );
    		ShapeFunctions reference_functions = ShapeFunctions( dim, order = 1 );
    		Quadrature quadrature;
    		quadrature = quadrature.get_quadrature(order, border);
    		
    		DenseMatrix Ke;
    		Ke.set_size(3,3);
    		
    		
    		assemble_elementary_matrix( elt_mapping, reference_functions, quadrature, unit_fct, Ke);
    		return true;
    	}
    }
}
