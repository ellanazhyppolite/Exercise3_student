#include <iostream>       // basic input output streams
#include <fstream>        // input output file stream class
#include <cmath>          // librerie mathematique de base
#include <iomanip>        // input output manipulators
#include <valarray>
#include "/Users/ella/Documents/EPFL/BA4/PHYSNUM/physnum_ex3/common/ConfigFile.h" // Il contient les methodes pour lire inputs et ecrire outputs 
//#include "/Users/eldidi/Desktop/physnum/physnum_ex3/physnum_ex3/common/ConfigFile.h" 
#include <numeric>

using namespace std; // ouvrir un namespace avec la librerie c++ de base


class GravitationalSimulation
{
private:
    // =========================
    // Physical and numerical parameters
    // =========================
    std::size_t numBodies = 0;
    int timeScheme = 0;              // 0 = fixed RK4, 1 = adaptive RK4
    int sampling = 1;
    double t = 0.0;
    double tEnd = 0.0;
    double dt = 0.0;
    double tolerance = 1e-6;
    double G = 0.0;
    // Optional atmosphere model
    bool useAtmosphere = false;
    double rho0 = 0.0;               // Density at reference radius
    double atmosphereScale = 1.0;    // Exponential scale height
    double dragArea = 0.0;           // Cross-sectional area
    double dragCoefficient = 0.0;    // Cx
    std::size_t dragBody = 0;        // Which moving body feels drag
    std::size_t dragCenterBody = 0;  // Which body has the atmosphere
    // Moving-body masses and state vector
    std::valarray<double> masses;    // size = numBodies
    std::valarray<double> state;     // size = 4*numBodies
    // Moving-body radii
    std::valarray<double> radii;   // size = numBodies
    // Work arrays for RK4
    std::valarray<double> k1, k2, k3, k4; // preallocate valarrays for the RK4 step parts
    std::valarray<double> trial1, halfStep, trial2, derivativeBuffer; // preallocate valarrays for the different steps and dy/dt vector
    // Force array
    std::valarray<double> Force;
    // autre
    double d;
    // Output
    std::ofstream outputFile;
    int stepsSinceLastOutput = 0;

    // Pour récupérer indices
    std::size_t ix(std::size_t i) const { return 2 * i; }
    std::size_t iy(std::size_t i) const { return 2 * i + 1; }
    std::size_t ivx(std::size_t i) const { return 2 * numBodies + 2 * i; }
    std::size_t ivy(std::size_t i) const { return 2 * numBodies + 2 * i + 1; }
    std::size_t Fijx(std::size_t i, std::size_t j) const { return 2 * i * numBodies + 2 * j; }
    std::size_t Fijy(std::size_t i, std::size_t j) const { return 2 * i * numBodies + 2 * j + 1; }
    



    // verifie que la distance entre chaque corps soit au moins la somme de leur rayon --> renvoie true si collision
    bool checkCollision(){
        for(std::size_t i = 0; i < numBodies ; ++i){
            for(std::size_t j = 0; j < numBodies ; ++j){
                if(i != j){
                    double norm_rij(sqrt( pow(state[ix(j)]-state[ix(i)],2) + pow(state[iy(j)]- state[iy(i)],2))); 
                    double  Rtot(radii[i] + radii[j]);
                    if(norm_rij <= Rtot){
                        return true;
                    }
                }
            }
        }
        return false;
    }
    //////ATTENTION A UTILISER DANS LA BOUCLE RUN 

    void printOut(bool force){
        if (force || stepsSinceLastOutput >= sampling) {
            outputFile << t << " ";

            for (size_t i = 0; i < numBodies; ++i) {
                outputFile << state[ix(i)] << " " << state[iy(i)] << " " << state[ivx(i)] << " " << state[ivy(i)] << " ";   //verifier si ajouter acceleration
        }

        outputFile << endl;

        stepsSinceLastOutput = 0;
        }
        else {
        stepsSinceLastOutput++;
        }
        //cout << "t = " << t << endl;
        //cout << "dt = " << dt << endl;
    }
/*
    // Force appliquée par j sur i
    std::valarray<double> Force(const std::valarray<double>& state_, std::size_t i, std::size_t j){
        std::valarray<double> Force;
        Force.resize(2);
        double dx = state_[ix(i)] - state_[ix(j)];
        double dy = state_[iy(i)] - state_[iy(j)];
        double norm_rij(sqrt( pow(dx,2) + pow(dy,2)));
        if(i != j){
            Force[0] = -G * masses[i] * masses [j] / pow(norm_rij,3) * dx;
            Force[1] = -G * masses[i] * masses [j] / pow(norm_rij,3) * dy;
         }else{
            Force[0] = 0;
            Force[1] = 0;
        }
        if((i == dragBody) and (j == dragCenterBody) and useAtmosphere){      //force de l'athmosphere de dragCenterBody sur dragBody
            double dvx = state_[ivx(dragBody)] - state_[ivx(dragCenterBody)];
            double dvy = state_[ivy(dragBody)] - state_[ivy(dragCenterBody)];
            double norm_vij(sqrt( pow(dvx,2) + pow(dvy,2)));
            double Force_atm_norm(-0.5 * rho0 * exp(-(norm_rij - radii[dragCenterBody])/ atmosphereScale) * dragArea * dragCoefficient * norm_vij) ;
            Force[0] += Force_atm_norm * dvx;
            Force[1] += Force_atm_norm * dvy;
        }
        return Force;
    }*/


        // Force appliquée par j sur i
    std::valarray<double> Force_calcul(const std::valarray<double>& state_){
        for(std::size_t i = 0 ; i < numBodies; ++i){
            for(std::size_t j = 0 ; j < numBodies ; ++j){
                double dx = state_[ix(i)] - state_[ix(j)];
                double dy = state_[iy(i)] - state_[iy(j)];
                double norm_rij(sqrt( pow(dx,2) + pow(dy,2)));
                if(i != j){
                    Force[Fijx(i,j)] = -G * masses[i] * masses [j] / pow(norm_rij,3) * dx;
                    Force[Fijy(i,j)] = -G * masses[i] * masses [j] / pow(norm_rij,3) * dy;
                }else{
                    Force[Fijx(i,j)] = 0;
                    Force[Fijy(i,j)] = 0;
                }
                if((i == dragBody) and (j == dragCenterBody) and useAtmosphere){      //force de l'athmosphere de dragCenterBody sur dragBody
                    double dvx = state_[ivx(dragBody)] - state_[ivx(dragCenterBody)];
                    double dvy = state_[ivy(dragBody)] - state_[ivy(dragCenterBody)];
                    double norm_vij(sqrt( pow(dvx,2) + pow(dvy,2)));
                    double Force_atm_norm(-0.5 * rho0 * exp(-(norm_rij - radii[dragCenterBody])/ atmosphereScale) * dragArea * dragCoefficient * norm_vij) ;
                    Force[Fijx(dragBody,dragCenterBody)] += Force_atm_norm * dvx;
                    Force[Fijy(dragBody,dragCenterBody)] += Force_atm_norm * dvy;
                }
            }
        }
        
        return Force;
    }



    // donne à derivativeBuffer la dérivée de state_ ( f(y) dans la notice )
    std::valarray<double> dydt(const std::valarray<double>& state_){
        for(std::size_t i = 0 ; i < numBodies ; ++i){
            derivativeBuffer[ix(i)] = state_[ivx(i)] ;
            derivativeBuffer[iy(i)] = state_[ivy(i)] ;
        }
        //derivativeBuffer[std::slice(0,1,numBodies*2)] = state_[std::slice(ivx(0),1,numBodies*2)];
        Force_calcul(state_);
        for(std::size_t i = 0 ; i < numBodies ; ++i){
                derivativeBuffer[ivx(i)] = 0;
                derivativeBuffer[ivy(i)] = 0;
            for(std::size_t j = 0; j < numBodies ; ++j){
                derivativeBuffer[ivx(i)] +=  1/masses[i] * Force[Fijx(i,j)];
                derivativeBuffer[ivy(i)] +=  1/masses[i] * Force[Fijy(i,j)];
            }
        }
        return derivativeBuffer;
    }


    // Norme d'un vecteur
    double norm(const std::valarray<double>& v){
        double n(0.0);
        for(size_t i = 0; i < v.size() ; ++i){
            n += pow(v[i],2);
        }
        return pow(n,0.5);
    }

    // tout les ki ont la bonne valeur après
    void RungeKutta(const std::valarray<double>& state_){
        k1 = dt*dydt(state_);
        k2 = dt*dydt(state_ + 0.5 * k1);
        k3 = dt*dydt(state_ + 0.5 * k2);
        k4 = dt*dydt(state_ + k3);
    }

    // passe de y_t à y_t+1 avec le schéma adaptatif
    void step(){
        if(timeScheme == 1){
            double dt_saved(0.0);
            do{
                RungeKutta(state);
                trial1 = state + 1.0/6.0 *  (k1 + 2*k2 + 2*k3 + k4);    //eq 2.141 notes de cours

                dt_saved = dt;
                dt = dt/2.0;
                RungeKutta(state);
                halfStep = state + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);
                RungeKutta(halfStep);
                trial2 = halfStep + 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);

                dt = dt_saved;
                d = norm(trial1 - trial2);

                cout << "d = " << d << endl;
                cout << "dt = " << dt << endl;

                if(d > tolerance){
                    dt = 0.9 * dt * pow(tolerance/d, 1.0/5.0); 
                }else{
                    dt = dt * pow(tolerance/d, 1.0/5.0);
                }
            }while(d > tolerance)

            cout << "%%%%%%" << endl;

            state = trial2;
            t += dt;

        }
        if(timeScheme == 0){
            RungeKutta(state);
            state += 1.0/6.0 * (k1 + 2*k2 + 2*k3 + k4);
            t += dt;
        }
    }


public:
// Initialisation de la classe
    GravitationalSimulation(int argc, char* argv[])
    {
        std::string inputPath = "configuration.in";
        if (argc > 1) {
            inputPath = argv[1];
        }
        ConfigFile configFile(inputPath);
        for (int i = 2; i < argc; ++i) {
            configFile.process(argv[i]);
        }
        // =========================
        // Required parameters
        // =========================
        numBodies  = static_cast<std::size_t>(configFile.get<int>("numBodies"));
        timeScheme = configFile.get<int>("timeScheme");
        sampling   = configFile.get<int>("sampling");
        tEnd       = configFile.get<double>("tEnd");
        dt         = configFile.get<double>("dt");
        tolerance  = configFile.get<double>("tolerance");
        G          = configFile.get<double>("G");
        // Optional switches
        useAtmosphere       = configFile.get<bool>("useAtmosphere", false);
 
        // Atmosphere / drag
        if (useAtmosphere) {
            rho0             = configFile.get<double>("rho0");
            atmosphereScale  = configFile.get<double>("atmosphereScale");
            dragArea         = configFile.get<double>("dragArea");
            dragCoefficient  = configFile.get<double>("dragCoefficient");
            dragBody         = static_cast<std::size_t>(configFile.get<int>("dragBody"));
            dragCenterBody   = static_cast<std::size_t>(configFile.get<int>("dragCenterBody", 0));
        }
        // Allocate arrays
        masses.resize(numBodies);
        radii.resize(numBodies);
        state.resize(4*numBodies);
        k1.resize(4 * numBodies);
        k2.resize(4 * numBodies);
        k3.resize(4 * numBodies);
        k4.resize(4 * numBodies);
        trial1.resize(4 * numBodies);
        halfStep.resize(4 * numBodies);
        trial2.resize(4 * numBodies);
        derivativeBuffer.resize(4*numBodies);
        Force.resize(2*numBodies*numBodies);
        // Read masses and initial conditions
        for (std::size_t i = 0; i < numBodies; ++i) {
            const std::string id = std::to_string(i + 1);
            masses[i]      = configFile.get<double>("m" + id);
            radii[i] = configFile.get<double>("r" + id, 0.0); // default 0 if not present
            state[ix(i)]   = configFile.get<double>("x" + id);
            state[iy(i)]   = configFile.get<double>("y" + id);
            state[ivx(i)]  = configFile.get<double>("vx" + id);
            state[ivy(i)]  = configFile.get<double>("vy" + id);
           
        }
        // Output file
        outputFile.open(configFile.get<std::string>("output").c_str());
        outputFile.precision(15);
        
    }

// run est une méthode de la classe GravitationalSimulation
    void run()
    {

        if(timeScheme == 0){cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;}
        else{cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;}

        stepsSinceLastOutput = 0;
        printOut(true);
        t = 0;

        while ((t < tEnd) and (not checkCollision())) {

            step();
            printOut(false);

        }
        
        if(checkCollision()){
            cerr << "Collision!" << endl;
        }

         printOut(true); // ensure final state is written
        if(timeScheme == 0){cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!" << endl;}
        else{cout << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$" << endl;}
    }
};
// Comprendre printOut !!!!!!!!!!!!!!!!!!!!!




int main(int argc, char* argv[])
{
    try {
        GravitationalSimulation simulation(argc, argv);
        simulation.run();
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
};

