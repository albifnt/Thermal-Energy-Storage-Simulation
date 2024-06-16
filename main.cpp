#include <iostream>
#include <fstream>
#include <cmath>
#define MIC 1000000
using namespace std;

// 1.INITIALIZATION
// CALL OF THE FUNTIONS
// Input function used to implement all the functions of the code
void input_OVS(); // Input of the Order Verification Study function
void input_Storage_Design_Study(); // Input for the Storage Design Study
void input_Simulation();// Input for the Exact solution

void TS_discretization(); // Function which discretizes time and space
void Status_determination(); // Function which identifies the status of the system at each time-step
void Determination_constants(); // Function which calculates the constants and variables necessary for the implementation of the simulation
void Temperature(); // Main function which calculates the temperature at each grid point and at each time step
void Order_Verification_Study(); // Function which executes the GRS
void Simulation();
void Storage_Design_Study(); // Function which carries on the simulation to determine the efficiency and the capacity factor of the TES using different dimensions


// DEFINITION OF THE DIMENSIONS OF THE TES
double height, diameter;
double V; // It is introduced just to carry on the Storage Design Study

// TEMPORAL AND SPATIAL DISCRETIZATION
int n_cells, total_number_step, time_steps_visualization;
double dx, dt;

// CHARACTERISTIC OF THE CYCLE
double T1, T2, T3, T4, time_cycle, number_cycle, n_dt_cycle; // n_dt_cycle = number of step per cycle
double dt1, dt2, dt3, dt4; // Number of steps per each status
int time_gap_visualization, status; //time_gap_visualization = it is the period of time after which you want to see the temperatures
int X; // Variable for the printing of the X axis
// VARIABLES NECESSARY TO EVALUATE THE STATUS OF THE TES
double  G[MIC], S[MIC]; // G[]=vector which contains the status of the TES; S[]=vector which contains the time;
int state_T1, state_T2, state_T3, state_T4;

// TEMPERATURES
double Tf[MIC],Ts[MIC],Tf_new[MIC],Ts_new[MIC]; // Temperatures fluid and solid new and old
double T_left_boundary; // Initial temperature of the charging phase: at point 1/2
double T_right_boundary; // Initial temperature of the discharging phase: at point N+1/2
double T_0,Tf_0,Ts_0; //Initial Temperature (Solid and liquid phase) defined by the user

// PHYSICAL PROPERTIES OF THE FLUID AND THE SOLID
double u_f; //Fluid velocity during charging and discharging
double alpha_f, alpha_s; // Fluid and solid diffusivities

// CONSTANTS RESPONSIBLE FOR THE STABILITY OF THE NUMERICAL METHOD
double s,d,c3;	//Calculated constants

// VARIABLE REST FOR THE IMPLEMENTATION OF THE FUNCTION TEMPERATURE
double REST;

// VARIABLES FOR THE IMPLEMENTATION OF THE ORDER VERIFICATION STUDY
double flag; // This flag is also the variable which activates the different functions of this code
double Sf,Ss; // Fluid and solid source terms
double steady_f,steady_f_old, steady_s, steady_s_old; // Variable used to check when the solution has reached a steady state
double tol_OVS; // Tolerance of the Order Verification Study
double Norm_inf, Norm_1; // Norms used to carry the OVS
double n, k; // n= number of wavelengths; k=wavenumber

// VARIABLES NECESSARY TO COUPLE THE FLUID AND SOLID DISCRETIZED EQUATIONS
double Tf_Star, Ts_Star;
double hv_f, hv_s;
double c_hv;
double ks; // ks= wavenumber for the solid manufactured solution


// ADDITIONAL PHYSICAL FLUID AND SOLID QUATITITES (they are introduce to make the comparison with the exact solution - Assignment 5)
double eps, rho_f, Cpf, rho_s, Cs;
double k_s, k_f;
double mu_f;
double ds;
double M_f; // Mass flow
double Re, Pr, Nu_fs, hfs, h, hv;

// VARIABLES USED TO CARRY ON THE STORAGE DESIGN STUDY
double exergy_efficiency, exergy_efficiency_old, Ec_flux_in, Ec_flux_out, Ed_flux_in, Ed_flux_out ;
double DOS, tol_DOS; // DOS is the variable, together with the flag, that activates the different functions of the code
double Q_c, Q_d, Q_max, Capacity_Factor; // Q_c = Thermal energy stored before charging - Q_d= Thermal energy stored before discharging


// COMPUTATION OF THE TEMPERATURE DIFFERENCE
double T_out_c;
double DELTA_T;

// OPENING FILES
// FILES FOR THE DETERMINATION OF THE STATUS
ofstream System_Status;
ofstream Time_Vector;

// FILES FOR CHARGING, DISCGHARGING AND IDLE PHASES
ofstream TemperatureFluid;
ofstream TemperatureSolid;
ofstream Space_Vector;

// FILES FOR OVS
ofstream Fluid_Norm;
ofstream Solid_Norm;
ofstream Refinement;
ofstream Manufactured_Fluid;
ofstream Manufactured_Solid;

// FISE FOR THE STORAGE DESIGN STUDY
ofstream Efficiency;
ofstream CapacityFactor;
ofstream Diameter;
ofstream Temperature_difference;


// 2. CALL OF THE FUNCTIONS INSIDE THE MAIN
int main()
{

//Status_determination();
//Simulation();
//Order_Verification_Study();
//Storage_Design_Study();

return 0;
}

// 3. DEFINITION OF THE FUNCTION

void TS_discretization()
{
    // Cycle duration
	time_cycle=T1+T2+T3+T4;

	// Control to implement the Order_Verification_Study function
    if (flag==1 && DOS == 0) n_dt_cycle=time_cycle/dt;
    else dt=time_cycle/n_dt_cycle;

	dt1=n_dt_cycle*(T1/time_cycle);
		if (dt1-floor(dt1)>= 0.5) dt1=ceil(dt1);
		else dt1=floor(dt1);

	dt2=n_dt_cycle*(T2/time_cycle);
		if (dt2-floor(dt2)>= 0.5) dt2=ceil(dt2);
		else dt2=floor(dt2);

	dt3=n_dt_cycle*(T3/time_cycle);
		if (dt3-floor(dt3)>= 0.5) dt3=ceil(dt3);
		else dt3=floor(dt3);

	dt4=n_dt_cycle*(T4/time_cycle);
		if (dt4-floor(dt4)>= 0.5) dt4=ceil(dt4);
		else dt4=ceil(dt4);

    //Number of steps which lie between two consecutive visualizations
    time_steps_visualization=int(time_gap_visualization/dt);

    // Total number of steps for the duration of analysis
    total_number_step=number_cycle*n_dt_cycle;

    // Grid spacing
	dx=height/n_cells;

}


void Status_determination(){

     // CALL OF THE INPUTS
     input_Simulation();

     // CALL OF THE DISCRETIZATION FUNCTION
     TS_discretization();

     // STATE VALUES THAT YOU WANT TO SEE IN THE PLOT
     state_T1=1;
     state_T2=2;
     state_T3=3;
     state_T4=4;

     for (int j=0; j < dt1; j++){
         G[j]= state_T1;
          }
     for (int j=dt1; j < dt1+dt2; j++){
         G[j]= state_T2;
          }
     for (int j=dt1+dt2; j < dt1+dt2+dt3; j++){
         G[j]= state_T3;
          }
     for (int j=dt1+dt2+dt3; j < n_dt_cycle; j++){
         G[j]= state_T4;}

     // Definition of the vector status of the system; it is gathered into a unique vector and it should be displayed in a graph
     // as a function of time
     System_Status.open("System_Status.txt");
     for (int z=0; z < number_cycle; z++){
     for (int i=0; i < n_dt_cycle; i++){
          System_Status << G[i] << " ";
     }
     }
     System_Status.close();

     // Definition of the vector time to display the status of the system
     Time_Vector.open("Time_vector.txt");
     for (int j=0; j < total_number_step; j++){
          S[j]=j*dt;
          Time_Vector << S[j] << " ";
     }
     Time_Vector.close();
}


void Determination_constants()
{
    if ((flag==0 && DOS== 0)||(flag==0 && DOS==1)){ // Devi farlo anche per la soluzione esatta
    u_f=M_f/(M_PI*pow(diameter,2)*eps*rho_f/4);

    Re= (eps*rho_f*u_f*ds)/mu_f;
    Pr= (mu_f*Cpf)/k_f;
    Nu_fs=(0.255/eps)*(pow(Pr,1.0/3.0))*(pow(Re,2.0/3.0));
    hfs=Nu_fs*k_f/ds;
    h=1/(1/hfs + ds/(10*k_s));
    hv=(6*(1-eps)/ds)*h;
    hv_f=hv/(eps*rho_f*Cpf);
    hv_s=hv/((1-eps)*rho_s*Cs);
    alpha_f=k_f/(eps*rho_f*Cpf);
    alpha_s=k_s/((1-eps)*rho_s*Cs);

    }


    s=u_f*dt/dx;
    d=alpha_f*dt/(pow(dx,2));
    c3=alpha_s*dt/(pow(dx,2));
    c_hv=1/(1 + hv_s*dt + hv_f*dt);
}

void Temperature()
{
    Determination_constants();

	//Open file
	if (flag== 0 && DOS==0){
	TemperatureFluid.open ("Temperature_Fluid.txt");
	TemperatureSolid.open ("Temperature_Solid.txt");
    Space_Vector.open("Space_Vector.txt");
	}

	if(s*s>s+2*d || s+2*d>1) {cout<<"THE SOLUTION IS NOT STABLE! CHANGE VALUES OF THE CFL AND THE FLUID DIFFUSION NUMBERS"<<endl; return;}
    if(2*c3>1) {cout<< "THE SOLUTION IS NOT STABLE! CHANGE VALUE OF THE SOLID DIFFUSION NUMBER"<<endl; return;}

	//Initialization
    for(int i=0; i<=n_cells; i++)
    {

    // MMS Implementation
    if (flag==1 && DOS==0 ) {Tf[i]=0; Ts[i]=0;}

    // Solution without order verification study
    else if (flag==0 && DOS==0){Tf[i]=Tf_0;Ts[i]=Ts_0;}

    } // Close initialization

    // Definition of the energy stored before the charging
    if (flag == 0 && DOS == 1){
    for (int i=1; i<=n_cells; i++){Q_c+=(M_PI/4*pow(diameter,2))*(eps*rho_f*Cpf*(Tf[i]-T_right_boundary)*dx + (1-eps)*rho_s*Cs*(Ts[i]-T_right_boundary)*dx);}
    }

    // Loop Temperatures
    for(int j=0; j< total_number_step; j++)
    {
    REST = j%int(n_dt_cycle);

    // CHARGING PHASE
    if (REST < dt1)
    {

    // Initialization of variables which determine if the solution is or is not steady
    if(flag==1 && DOS==0){
    steady_f_old=0;
    steady_s_old=0;
    }


            for (int i=1; i<=n_cells;i++)
            {
                if (flag==1 && DOS==0){
                Sf=(dt/dx)*(u_f*(cos(k*dx*(i))-cos(k*dx*(i-1)))+(alpha_f*k + hv_f/k)*(sin(k*dx*(i))-sin(k*dx*(i-1))) - hv_f/ks*(sin(ks*dx*(i)) - sin(ks*dx*(i-1))));
                Ss=(dt/dx)*((alpha_s*ks + hv_s/ks)*(sin(ks*dx*(i))-sin(ks*dx*(i-1)))-hv_s/k*(sin(k*dx*(i))-sin(k*dx*(i-1))));
                }

                if(i==1)//Left Boundary Condition
                {
                Tf_Star=Tf[1]-s*(Tf[1]-T_left_boundary)+d*(Tf[2]-Tf[1])+flag*Sf;
                Ts_Star=Ts[1]+c3*(Ts[2]-Ts[1])+flag*Ss;
                }

                else if(i==n_cells)//Right Boundary Condition
                {
                Tf_Star=Tf[n_cells]-s*(Tf[n_cells]-Tf[n_cells-1])+d*(Tf[n_cells-1]-Tf[n_cells])+flag*Sf;
                Ts_Star=Ts[n_cells]+c3*(Ts[n_cells-1]-Ts[n_cells])+flag*Ss;
                }

                else
                {
                Tf_Star=Tf[i]-s*(Tf[i]-Tf[i-1])+d*(Tf[i+1]-2*Tf[i]+Tf[i-1])+flag*Sf;
                Ts_Star=Ts[i]+c3*(Ts[i+1]-2*Ts[i]+Ts[i-1])+flag*Ss;
                }

                Tf_new[i]= c_hv*((1 + hv_s*dt)*Tf_Star + hv_f*dt*Ts_Star);
                Ts_new[i]= c_hv*(hv_s*dt*Tf_Star + (1 + hv_f*dt)*Ts_Star);

            // Steady solution for the fluid
            if (flag == 1 && DOS==0){
            steady_f = fabs(Tf_new[i] - Tf[i]);
            if (steady_f > steady_f_old){steady_f_old=steady_f;}

            // Steady solution for the solid
            steady_s = fabs(Ts_new[i] - Ts[i]);
            if (steady_s > steady_s_old){steady_s_old=steady_s;}
            } // Close if steady solution

			} //Close Loop Space

            // Calculation of the energy fluxes of the Storage_Design_Study
            if (flag==0 && DOS ==1){
            Ec_flux_in+= M_f*Cpf*(Tf_new[1] - T_0 - T_0*log(Tf_new[1]/T_0))*dt;
            Ec_flux_out+= M_f*Cpf*(Tf_new[n_cells] - T_0 - T_0*log(Tf_new[n_cells]/T_0))*dt;
            }

            // Calculation of the errors
            // FLUID
			if (flag==1 && DOS==0 && (steady_f_old <= tol_OVS) && (steady_s_old <= tol_OVS))
        	{
        		//Fluid Error
        		Norm_inf=0;
        		Norm_1=0;
        		for(int i=1;i<=n_cells;i++)
        		{
                Norm_1+=fabs(Tf_new[i]-cos(k*dx*(i-0.5)));
                if(fabs(Tf_new[i]-cos(k*dx*(i-0.5)))>Norm_inf)
                {Norm_inf=fabs(Tf_new[i]-cos(k*dx*(i-0.5)));}
				}
				Norm_1=Norm_1/n_cells;
                // Data for the definition of the actual order of accuracy
		    	Fluid_Norm << log(Norm_1) << " " << log(Norm_inf) << endl;
		    	Refinement << log(dx) << endl;
				//File for the visualization of the manufactured
				for(int i=1; i<=n_cells; i++)
				{Manufactured_Fluid<<Tf_new[i]<<" "<<cos(k*dx*(i-0.5))<<" "<<i*dx<<endl;}

				//Solid Error
        		Norm_inf=0;
        		Norm_1=0;
        		for(int i=1;i<=n_cells;i++)
        		{
                Norm_1+=fabs(Ts_new[i]-cos(ks*dx*(i-0.5)));
                if(fabs(Ts_new[i]-cos(ks*dx*(i-0.5)))>Norm_inf)
                {Norm_inf=fabs(Ts_new[i]-cos(ks*dx*(i-0.5)));}
				}
				Norm_1=Norm_1/n_cells;
                // Data for the definition of the actual order of accuracy
                Solid_Norm << log(Norm_1) << " " << log(Norm_inf) << endl;
				// Filefor the visualization of the manufactured
				for (int i=1; i<=n_cells; i++)
				{Manufactured_Solid<<Ts_new[i]<<" "<<cos(ks*dx*(i-0.5))<<" "<<i*dx<<endl;}
                cout << "DONE" << endl;
				return; //
        	} // Close check steady solution

			// New temperatures
            for (int i=1; i<=n_cells; i++)
            {Tf[i]=Tf_new[i];	Ts[i]= Ts_new[i];

            // Print temperature
			if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0){TemperatureFluid << Tf[i] << " "; TemperatureSolid << Ts[i] << " ";} // Close print temperatures
            if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0 && X == 0){Space_Vector << i*dx << " ";}
            }

            // End to the Temperature Fluid an Solid file
            if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0){TemperatureFluid << endl; TemperatureSolid << endl;} // Close print temperatures
            if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0 && X == 0){X=1;}

    } // Close Rest <= dt1 (charging phase)
     // This return is written because assignment 5 aims to the calculation of the single charging

    // COMPUTATION OF THE DIFFERENCE
    if (flag == 0 && DOS == 1 && j==dt1-1){T_out_c=Tf[n_cells];}

    // IDLE PHASE (charging --> idle)
    else if (REST >= dt1 && REST < dt1+dt2)
    {
		for (int i=1; i<=n_cells;i++)
		{

			if (i==1)//Left Boundary Condition
			{
            Tf_Star=Tf[1]+d*(Tf[1+1]-Tf[1]);
			Ts_Star=Ts[1]+c3*(Ts[1+1]-Ts[1]);
			}

			else if (i==n_cells)//Right Boundary Condition
			{
            Tf_Star=Tf[n_cells]+d*(Tf[n_cells-1]-Tf[n_cells]);
			Ts_Star=Ts[n_cells]+c3*(Ts[n_cells-1]-Ts[n_cells]);
			}

			else
			{
			Tf_Star=Tf[i]+d*(Tf[i+1]-2*Tf[i]+Tf[i-1]);
			Ts_Star=Ts[i]+c3*(Ts[i+1]-2*Ts[i]+Ts[i-1]);
			}

            Tf_new[i]= c_hv*((1 + hv_s*dt)*Tf_Star + hv_f*dt*Ts_Star);
			Ts_new[i]= c_hv*(hv_s*dt*Tf_Star + (1 + hv_f*dt)*Ts_Star);

            } // Close loop Space

            // New temperatures
            for (int i=1; i<=n_cells; i++)
            {Tf[i]=Tf_new[i];	Ts[i]= Ts_new[i];

            // Print temperature
			if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0){TemperatureFluid << Tf[i] << " "; TemperatureSolid << Ts[i] << " ";} // Close print temperatures
            }

            // End to the Temperature Fluid an Solid file
            if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0){TemperatureFluid << endl; TemperatureSolid << endl;} // Close print temperatures
    } // Close Rest <= dt1+dt2 (idle phase)

    // Definition of the energy stored before the discharging
    if (flag==0 && DOS== 1 && REST==dt1+dt2-1){
    for (int i=1; i<=n_cells; i++){Q_d+=(M_PI/4*pow(diameter,2))*(eps*rho_f*Cpf*(Tf[i]-T_right_boundary)*dx + (1-eps)*rho_s*Cs*(Ts[i]-T_right_boundary)*dx);}
    }

    // DISCHARGING PHASE
    else if (REST >= dt1+dt2 && REST < dt1+dt2+dt3)
    {

		for (int i=1; i<=n_cells;i++)
		{

			if (i==1)//Left Boundary condition
			{
            Tf_Star=Tf[1]+s*(Tf[1+1]-Tf[1])+d*(Tf[1+1]-Tf[1]);
			Ts_Star=Ts[1]+c3*(Ts[1+1]-Ts[1]);
			}

			else if (i==n_cells)//Right Boundary Condition
			{
            Tf_Star=Tf[n_cells]+s*(T_right_boundary-Tf[n_cells])+d*(Tf[n_cells-1]-Tf[n_cells]);
			Ts_Star=Ts[n_cells]+c3*(Ts[n_cells-1]-Ts[n_cells]);
			}

			else
			{
			Tf_Star=Tf[i]+s*(Tf[i+1]-Tf[i])+d*(Tf[i+1]-2*Tf[i]+Tf[i-1]);
			Ts_Star=Ts[i]+c3*(Ts[i+1]-2*Ts[i]+Ts[i-1]);
			}

            Tf_new[i]= c_hv*((1 + hv_s*dt)*Tf_Star + hv_f*dt*Ts_Star);
			Ts_new[i]= c_hv*(hv_s*dt*Tf_Star + (1 + hv_f*dt)*Ts_Star);

		} // Close loop Space

            // Calculation for the evaluation of the Storage_Design_Study
            if (flag==0 && DOS==1){
            Ed_flux_out+= M_f*Cpf*(Tf_new[1] - T_0 - T_0*log(Tf_new[1]/T_0))*dt;
            Ed_flux_in+= M_f*Cpf*(Tf_new[n_cells] - T_0 - T_0*log(Tf_new[n_cells]/T_0))*dt;
            }

            // New temperatures
            for (int i=1; i<=n_cells; i++)
            {Tf[i]=Tf_new[i];	Ts[i]= Ts_new[i];

            // Print temperature
			if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0){TemperatureFluid << Tf[i] << " "; TemperatureSolid << Ts[i] << " ";} // Close print temperatures
            }

            // End to the Temperature Fluid an Solid file
            if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0){TemperatureFluid << endl; TemperatureSolid << endl;} // Close print temperatures

    } // Close REST <= dt1+dt2+dt3 (discharging phase)

    // IDLE PHASE (discharging --> idle)
    else if (REST >= dt1+dt2+dt3 && REST < n_dt_cycle)
    {

		for (int i=1; i<=n_cells;i++)
		{
			if (i==1) //Left Boundary Condition
			{
            Tf_Star=Tf[1]+d*(Tf[1+1]-Tf[1]);
			Ts_Star=Ts[1]+c3*(Ts[1+1]-Ts[1]);
			}

			else if (i==n_cells)//Right Boundary Condition
			{
            Tf_Star=Tf[n_cells]+d*(Tf[n_cells-1]-Tf[n_cells]);
			Ts_Star=Ts[n_cells]+c3*(Ts[n_cells-1]-Ts[n_cells]);
			}//Central cells

			else
			{
            Tf_Star=Tf[i]+d*(Tf[i+1]-2*Tf[i]+Tf[i-1]);
			Ts_Star=Ts[i]+c3*(Ts[i+1]-2*Ts[i]+Ts[i-1]);
			}

            Tf_new[i]= c_hv*((1 + hv_s*dt)*Tf_Star + hv_f*dt*Ts_Star);
			Ts_new[i]= c_hv*(hv_s*dt*Tf_Star + (1 + hv_f*dt)*Ts_Star);

		} // Close loop space

            // New temperatures
            for(int i=1; i<=n_cells; i++)
            {Tf[i]=Tf_new[i];	Ts[i]= Ts_new[i];

            // Print temperature
			if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0){TemperatureFluid << Tf[i] << " "; TemperatureSolid << Ts[i] << " ";} // Close print temperatures
            }

            // End to the Temperature Fluid an Solid file
            if (j%time_steps_visualization == 0 && flag == 0 && DOS == 0){TemperatureFluid << endl; TemperatureSolid << endl;} // Close print temperatures

    } // Close REST <= n_dt_cycle (idle phase)

    } // Close total number of step j

    if (flag==0 && DOS==0){
    TemperatureFluid.close();
	TemperatureSolid.close();
	Space_Vector.close();
    }

} // Close temperature function



void Order_Verification_Study()
{
	Manufactured_Fluid.open ("Manufactured_Fluid.txt");
	Manufactured_Solid.open("Manufactured_Solid.txt");
	Fluid_Norm.open("Fluid_Norm.txt");
	Solid_Norm.open("Solid_Norm.txt");
	Refinement.open("Refinement.txt");
	input_OVS();
    TS_discretization();
	for (int grid_refinement=0; grid_refinement<8; grid_refinement++)
	{
		cout << grid_refinement << " ";

        dx=height/n_cells;

		Temperature();

		n_cells=2*n_cells;
	}
	Manufactured_Fluid.close();
	Manufactured_Solid.close();
    Fluid_Norm.close();
	Solid_Norm.close();
	Refinement.close();
}


void input_OVS()
{
    // DETERMINATION OF THE FLAG
    flag=1;
    DOS=0;

    // STATE DURATION
	T1=1000;
	T2=0;
	T3=0;
	T4=0;

	// 0SPATIAL DISCRETIZATION
	n_cells=10;

    // TEMPORAL DISCRETIZATION
    dt=0.005;
	number_cycle=1;
	time_gap_visualization=10;

    // TEMPERATURES
    // Initial temperature of fluid and solid
	T_0=1;
    Tf_0=T_0;
	Ts_0=T_0;
    // Boundary temperatures
	T_left_boundary=T_0;
	T_right_boundary=T_0;

	// DIMENSIONS OF THE THERMAL ENERGY STORAGE
    height=1;
	diameter=1;

	// FLUID AND SOLID PROPERTIES (Assignment 5)
    // Assignment 5
    eps=0.4;
    ds=0.03;
    rho_s=2600.0;
	rho_f=1835.6;
	Cs=900.0;
	Cpf=1511.8;
    k_s=2.0;
    k_f=0.52;
    mu_f=2.63;
    M_f=0.1;

    // Fluid velocity and fluid and solid diffusivity
    u_f=0.1;
	alpha_f=2e-7;
	alpha_s=9e-7;
	// dEFINISCILE
    hv_f=1000;
	hv_s=1000;

	// VARIABLES FOR THE IMPLEMENTATION OF THE ORDER VERIFICATION STUDY
	n=1;
	k=(2*M_PI*n)/height;
	ks=2*k; // Kappa solid
    tol_OVS=1e-20;
}


void Storage_Design_Study()
{
    Efficiency.open("Efficiency.txt");
    CapacityFactor.open("CapacityFactor.txt");
    Diameter.open("Diameter.txt");
    Temperature_difference.open("Temperature_difference.txt");

    // CALL OF THE IMPUTS (flag=0 && DOS=1)
    input_Storage_Design_Study();

    // Definition of the tolerance
    tol_DOS=1e-6;

    // Iteration
    int k;

    // Cycle for the variation of diameter
    for (int z=0; z<=4; z++)
    {
    diameter+=1; // increase of the diameter (for the storage design study the initial diameter is 3)

    // Iteration
    k=0;

    // Calculation of the height of the TES
    height=V/(M_PI*pow(diameter,2)/4);

    // Call of TS_discretization
    TS_discretization();

    // Initialization of the energy_efficiencies
    exergy_efficiency_old=0;
    exergy_efficiency=100000;

    // Initialization temperatures
    for(int i=0; i<=n_cells; i++){Tf[i]=Tf_0;Ts[i]=Ts_0;}

    // Definition of the maximum capacity factor
    Q_max=(eps*rho_f*Cpf+(1-eps)*rho_s*Cs)*(M_PI/4*pow(diameter,2))*height*(T_left_boundary - T_right_boundary);

    while (abs(exergy_efficiency - exergy_efficiency_old) > tol_DOS){
    exergy_efficiency_old=exergy_efficiency;

    // Iteration
    k+=1;

    // Call of the function Temperature
    Temperature();

    // Evaluation of the energy efficiency
    exergy_efficiency=(Ed_flux_out-Ed_flux_in)/(Ec_flux_in-Ec_flux_out); // Hai invertito perchè out charging = in discharging
    // File for the energy efficiency
    Efficiency << exergy_efficiency << " " << k <<endl;

    // Re-Initialization of all the energy fluxes
    Ec_flux_in=0;    // Charging in
    Ec_flux_out=0;    // Charging out
    Ed_flux_in=0;    // Discharging in
    Ed_flux_out=0;    // Discharging out

    // Computation of the temperature difference
    DELTA_T=T_out_c - T_right_boundary;
    Temperature_difference << DELTA_T << " " << k << endl;

    // Definition of the capacity factor
    Capacity_Factor=(Q_d - Q_c)/Q_max;
    // File for the Capacity Factor
    CapacityFactor << Capacity_Factor << " " << k << endl;

    // Re-Initialization of all the capacity factors
    Q_c=0;    // Capacity factor before charging
    Q_d=0;    // Capacity factor before discharging
    }
    Diameter << diameter << endl;

    } // Close for of the diameters
}


// DECISION OF WHAT YOU WANT TO DO
    // (flag == 0 && DOS == 0) ---> Temperature Function (the function status also)
    // (flag == 1 && DOS == 0) ---> Order Verification Study
    // (flag == 0 && DOS == 1) ---> Storage Design Study


void input_Storage_Design_Study()
{
    // DETERMINATION OF THE FLAG
    flag=0;
    DOS=1;

    // STATE DURATION
	T1=3600*6;	// CHARGING
	T2=3600*6;	    // IDLE (charging-idle)
	T3=3600*6;	// DISCHARGING
	T4=3600*6; 	// IDLE (discharging - idle)

	// SPATIAL DISCRETIZATION
	n_cells=2000;

    // TEMPORAL DISCRETIZATION
    time_gap_visualization=10000; // The time_gap_visualization is high in order not to print too many temperatures
	number_cycle=1;
	n_dt_cycle=24*3600;

    // TEMPERATURES
    // Initial temperature of fluid and solid
	T_0=288.15;
    Tf_0=293;
	Ts_0=293;
    // Boundary temperatures
	T_left_boundary=873;
	T_right_boundary=293;

	// DIMENSIONS OF THE THERMAL ENERGY STORAGE
	diameter=3;
    V=300;

	// FLUID AND SOLID PROPERTIES (Assignment 5)
    // Assignment 5
    eps=0.4;
    ds=0.03;
    rho_s=2600.0;
	rho_f=1835.6;
	Cs=900.0;
	Cpf=1511.8;
    k_s=2.0;
    k_f=0.52;
    mu_f=2.63;
    M_f=10;

}


void Simulation()
{
    // VARIABLE FOR THE PRINTING OF THE X AXIS
    X=0;

    // CALL OF THE INPUTS
    input_Simulation();

    // CALL OF TS_discretization
    TS_discretization();

    // CALL OF THE TEMPERATURE FUNTION
    Temperature();
}


void input_Simulation()
{
    // DETERMINATION OF THE FLAG
    flag=0;
    DOS=0;

    // STATE DURATION
    T1=5000;
	T2=0;
	T3=0;
	T4=0;

	// SPATIAL DISCRETIZATION
    n_cells=2000;

    // TEMPORAL DISCRETIZATION
    time_gap_visualization=100;
    number_cycle=1;
    n_dt_cycle=500000;

    // TEMPERATURES
    // Initial temperature of fluid and solid
    T_0=293;
    Tf_0=T_0;
	Ts_0=T_0;
    // Boundary temperatures
    T_left_boundary=873;
	T_right_boundary=293;

	// DIMENSIONS OF THE THERMAL ENERGY STORAGE
	//V=300;
	diameter=1;
    // Calculation of the height of the TES
    height=1;
    //height=V/(M_PI*pow(diameter,2)/4);

	// FLUID AND SOLID PROPERTIES (Assignment 5)
    eps=0.4;
    ds=0.03;
    rho_s=2600.0;
	rho_f=1835.6;
	Cs=900.0;
	Cpf=1511.8;
    k_s=2.0;
    k_f=0.52;
    mu_f=2.63;

	// Mass flow rate
    M_f=0.1;
}
