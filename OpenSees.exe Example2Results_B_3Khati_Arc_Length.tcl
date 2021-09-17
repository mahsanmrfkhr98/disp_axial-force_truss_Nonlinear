# Clear Previous Data
wipe
# units: lb, psi, in
# Shekl, safheye 1 az file "Session 02.pdf".
# Define Model Properties: model  BasicBuilder  -ndm  Number_of_Dimensions  -ndf  Nunber_of_DoFs_in_Each_Node
						   model  BasicBuilder  -ndm         2              -ndf              2
						   
# Define Nodes Information
	# Create nodes: node  node#  x[in]  y[in]
					node    1     0.0    0.0
					node    2     7.0    0.0
					node    3    12.0    4.0
					node    4     5.0    8.0
	
	# Assign nodes to a Group: region  Region#  -Type  Target_Nodes 
				   			   region    1      -node    1 2 3 4
	
	# Set Restrained DoFs: fix  Node#  X_Restrained?  Y_Restrained?
						   fix    1         1                1
					       fix    2         1                1
					       fix    3         1                1

# Define Elements Information
	#                                                Shekl, safheye 12 az file "Session 02.pdf".
    #                                                http://opensees.berkeley.edu/wiki/index.php/ElasticMultiLinear_Material
	# Create Multilinear Material: uniaxialMaterial  ElasticMultiLinear  Material#  -strain       Points_Strains[-]       -stress    Points_Stresses[psi]
	#                                                                                        ---------------------------           ------------------------
	#                                                                                         0    ey    ec=4ey  ef=6ey             0    Fy       Fc     0
                                   uniaxialMaterial  ElasticMultiLinear      1      -strain   0  0.0008  0.0032  0.0048   -stress   0  24000.0  26880.0  0
								   uniaxialMaterial  ElasticMultiLinear      2      -strain   0  0.0008  0.0032  0.0048   -stress   0  24000.0  26880.0  0
								   uniaxialMaterial  ElasticMultiLinear      3      -strain   0  0.0008  0.0032  0.0048   -stress   0  24000.0  26880.0  0
											  
	# Create Truss Elements:       element  truss  Element#  Start_Node  End_Node  A[in^2]  Assigned_Material#
							       element  Truss     1          1           4       6.0            1
							       element  Truss     2          2           4       4.0            2
							       element  Truss     3          3           4       4.0            3
	
	# Assign Elements to a Group:  region  Region#  -Type  Target_Elements 
				   			       region    2      -ele       1 2 3 
	
# Define Applied Loads
	# Create a Linear Loading TimeSeries: timeSeries  Linear  TimeSeries#
								          timeSeries  Linear      1
									 
	# Create a Plain load pattern associated with the TimeSeries: pattern  Plain  pattern#  Assigned_TimeSeries#  {
																  pattern  Plain     1              1             {
	# Create the nodal loads:																					     load  Node#  X_Value[lb]  Y_Value[lb]
																												     load    4        2.0         1.7
																												  }

# Define Analysis Parameters																							  
							system SparseGeneral
							numberer RCM
							constraints Transformation 						
							integrator ArcLength 0.02  0.01				
                            test EnergyIncr  1e-5  500
							algorithm Newton
							analysis Static

# Create Output Files
	# Output File for Nodal Info:    recorder  Node     -file  Output_File_Name          -time[?]  -node/Region  Region#  -dof  Target_DoFs  Type
								     recorder  Node     -file  OSNodalDisplacements.txt            -region          1     -dof      1 2      disp  
									                                                   
	# Output File for Elements Info: recorder  Element  -file  Output_File_Name          -time[?]  -ele/Region   Region#                     Type				   
 									 recorder  Element  -file  OSElementsForce.txt                 -region          2                        basicForces
									 
# Start the Analysis: analyze  Number_of_Steps									 
					  analyze       150001