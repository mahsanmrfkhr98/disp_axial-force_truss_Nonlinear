#clear previous data
wipe

#unit:lb,psi,in
#define model properties:model   BasicBuilder -ndm number_of_Dimensions -ndf  Nunber_of_DoFs_in_Each_Node
                          model  BasicBuilder -ndm  2                   -ndf        2
 #Define nodes information:
   #create nodes:  node      node#   x[in]   y[in] 
                   node       1       0.0    12.0
                   node       2       0.0     0.0 
				   node       3      17.0     0.0
				   node       4      24.0    12.0
				   node       5       9.0    17.0
    #Asign nodes to a Group:region  region#   -Type      Target_node	
 
    #set restrained DOFs: 	fix     node#  x_restrained	   y_restrained
                             fix      1          1                1
                             fix      2          1                1
                             fix      3          1                1
                             fix   	  4		     1                1	
# Define Elements Information                                            
	# Create Bilinear Material: uniaxialMaterial  Steel01  Material#  Fy[psi]    E[psi]     alpha[-]  
						        uniaxialMaterial  Steel01     1       25000.0  36000000.0    0.03
    #create Truss Elements:element  truss  Element#  Start_Node  End_Node  A[in^2]    Assigned_Material#
                           element  truss     1           1         5       5                 1
						   element  truss     2           2         5       4                 1
						   element  truss     3           3         5       5                 1
						   element  truss     4           4         5       6                 1
	# Assign Elements to a Group:  region  Region#  -Type  Target_Elements 
				   			       region    2      -ele       1 2 3 4
# Define Applied Loads
	# Create a Linear Loading TimeSeries: timeSeries   Linear  TimeSeries#					
	                                      timeSeries   Linear      1
	
	#Create a Plain Load Pattern associated with TimeSeries:Pattern  Plain   Pattern#
	 pattern  Plain     1              1             {
	# Create the nodal loads:																					    load  Node#  X_Value[lb]  Y_Value[lb]
																												    load    5        3.0         2.5
																												  }

# Define Analysis Parameters
							system SparseGeneral
							numberer RCM
							constraints Transformation 						
							integrator LoadControl 1.0					
                            test EnergyIncr 1e-5 500
							algorithm Newton
							analysis Static
	# Create Output Files
	# Output File for Nodal Info:    recorder  Node     -file   Output_File_Name          -time[?]  -node/Region  Node/Region#  -dof  Target_DoFs  Type
								     recorder  Node     -file   OSNodalDisplacements.txt            -node           1 2 3 4 5          -dof      1 2      disp   	
									 
# Start the Analysis: analyze  Number_of_Steps									 
					  analyze       150000						
							
																		
	