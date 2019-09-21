#include <iostream>
#include <math.h>
#include <string>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <time.h>
#include <mpi.h>
#include <iomanip>


 
/* Résolution de Laplacien(u) = 0 par un algorithme de Monte-Carlo */

using namespace std;


int main(int argc, char ** argv)
{
	//ajout de l'identifiant et du nombre de processeurs	
	MPI_Init(&argc, &argv);
	int rank;
	int nb_proc;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nb_proc);

	double t1 = MPI_Wtime(); //Valeur du temps au début

	// Lecture des données d'entrée dans le fichier passé en paramètre
	// argc = nombre d'arguments passés en ligne de commande. Vaut 1 au minimum (nom de l'exécutable)

	if ( argc < 2 )
	{
		if ( rank == 0 )
		  cout << "Il manque le nom du fichier de données !" << endl;
		exit(1);
        }

	int nx, ny, nb_tirages;
	double temps_init, temps_final;

	// Lecture du fichier de données par les processeur de rang 0
	// argv est un tableau de chaînes de caractères passées en ligne de commande. argv[0] contient le nom de l'éxecutable. 
	
	
	
	if ( rank == 0 )
	{
		ifstream fic_donnees;
		fic_donnees.open(argv[1], ifstream::in);

		fic_donnees >> nx >> ny >> nb_tirages;
		fic_donnees.close();

		// Affichage des données lues : permet de vérifier que tout va bien.

		cout << "************************" << endl;
		cout << "Calcul exécuté sur " << nb_proc << " coeurs" << endl;
		cout << "Dimension selon x : " << nx << endl;
		cout << "Dimension selon y : " << ny << endl;
		cout << "Nombre de lancés par case : " << nb_tirages << endl;
		for (int j=1; j<nb_proc; j++)
			{
			MPI_Send(&nx, 1, MPI_INT, j, 1, MPI_COMM_WORLD);
			MPI_Send(&ny, 1, MPI_INT, j, 2, MPI_COMM_WORLD);
			MPI_Send(&nb_tirages, 1, MPI_INT, j, 3, MPI_COMM_WORLD);
			}
	}

	else{
		
		MPI_Status status;
		MPI_Recv(&nx, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		MPI_Recv(&ny, 1, MPI_INT, 0, 2, MPI_COMM_WORLD, &status);
		MPI_Recv(&nb_tirages, 1, MPI_INT, 0, 3, MPI_COMM_WORLD, &status);
	}

	MPI_Barrier(MPI_COMM_WORLD);

	//Partition de j
	int div = (ny - 2) / nb_proc;
	int reste = (ny - 2) % nb_proc;
	int debut_tab = div * rank + 1;
	int fin_tab = div * (rank + 1);
	
	int r;
	for ( r=0; r<reste; r++){
		if (rank==r){
			fin_tab = fin_tab + 1;
		}
		if (rank>r){
			fin_tab = fin_tab + 1;
			debut_tab = debut_tab + 1;
		}
	}

	cout << "Proc : " << rank << ", debut = " << debut_tab << ", fin = " << fin_tab << endl;

	// Initialisation et allocation
	
	double* grille = new double[nx * ny];
	
	int i, j;

	double valeurs_aux_bords[4] = { 0,1,2,3 }; // valeurs de la solution sur les bords y=0, x=nx, y=ny et x=0
	
	for (i=0; i<nx*ny; i++)
		grille[i] = 0;  // Initialisation de la solution
			
	// Initialisation de la solution sur les bords
	
	for (i = 0; i < nx; i++) // y=0
		grille[i] = valeurs_aux_bords[0];
	
	for (i = 1; i < ny - 1; i++) // x=nx
		grille[i * nx + nx - 1] = valeurs_aux_bords[1];
	
	for (i = 0; i < nx; i++) // y=ny
		grille[(ny - 1) * nx + i] = valeurs_aux_bords[2];

	for (i = 1; i < ny - 1; i++) // x=0
		grille[i * nx] = valeurs_aux_bords[3];

	
	srand(time(NULL)); // Initialisation du générateur de nombres aléatoire

	// 1. Boucles sur j et i
	
	for (j = debut_tab; j < fin_tab + 1; j++)
	{
	  for (i = 1; i < nx - 1; i++) 
	  {
	    for (int n = 0; n < nb_tirages; n++) // on effectue nb_tirages tirages de particules
	    {
			int position[2] = { i,j }; // Initialisation de la position de la particule
			bool stop = 0; // variable booléenne qui passe à 1 quand on rencontre un bord
			double valeur = 0;
			while ( stop == 0 ) // tant que stop est nul, faire ...
			{
				int decision = rand() % 2; // nombre aléatoire [0,1]
		  
				if (decision == 0)
					position[0]--;
				else
					position[0]++;

				decision = rand() % 2;
				if (decision == 1)
					position[1]--;
				else
					position[1]++;

				// Test des bords : est-on sur un bord, si oui lequel ?

				if ( ( position[0] == 0 ) || ( position[0] == nx-1 ) || ( position[1] == 0 ) || ( position[1] == ny-1 ) )
				{
					valeur = grille[position[1] * nx + position[0]]; // On récupère la valeur au bord atteint
					stop = 1;
				}
			}
	      
		  grille[j * nx + i] += valeur;
	    }
	    
	    // Valeur finale de la solution en (i,j) : c'est la moyenne des tirages
	    
		grille[j * nx + i] /= nb_tirages;
	  }
	}

	MPI_Barrier(MPI_COMM_WORLD);
	//reconstruction de la grille

	if (rank != 0)// Tous les autres processeurs envoie leur partie du tableau 
	{
		MPI_Send(&grille[int(debut_tab * nx)], nx* (fin_tab - debut_tab + 1), MPI_DOUBLE, 0, 1, MPI_COMM_WORLD);
	}
	else
	{	
		// Le processeur 0 récupère les parties de chaque grille qu'il associe à sa grille à la bonne place
		for (int k = 1; k < nb_proc; k++)
		{
			int debut = div * k + 1;
			int fin = div * (k + 1);
			for ( r = 0; r < reste; r++){
				if (k == r){
					fin = fin + 1;
				}
				if (k>r){
					fin = fin + 1;
					debut = debut + 1;
				}
			}
			MPI_Status status;
			MPI_Recv(&grille[debut * nx], nx* (fin - debut + 1), MPI_DOUBLE, k, 1, MPI_COMM_WORLD, &status);
		}
	
	}

if (rank == 0){
for (i = 0; i < nx * ny; i++)
		{
		cout << "proc" << rank << "rang : " << i << " valeur : " << grille[i] << endl;
		}
}

	
MPI_Barrier(MPI_COMM_WORLD);
double t2 = MPI_Wtime();
cout << "temps : " << t2 - t1 << endl;
	
	// 2. Ecriture au format vtk
	if (rank == 0){

		cout << "Écriture du fichier vtk" << endl;
		stringstream nom;
		nom << "stochastique";
		int Nbnoe = nx*ny;
		double dx = 1.0 / (nx - 1);
		double dy = 1.0 / (ny - 1);

		nom << ".vtk" << '\0';
		ofstream fic(nom.str().c_str());

		fic << "# vtk DataFile Version 2.0" << endl;
		fic << "Laplacien stochastique" << endl;
		fic << "ASCII" << endl;
		fic << "DATASET STRUCTURED_POINTS" << endl;
		fic << "DIMENSIONS " << nx << "  " << ny << "  1 " << endl;
		fic << "ORIGIN 0 0 0" << endl;
		fic << "SPACING " << dx << "  " << dy << "  1" << endl;
		fic << "POINT_DATA " << Nbnoe << endl;
		fic << "SCALARS Concentration float" << endl;
		fic << "LOOKUP_TABLE default" << endl;
		for (i = 0; i < Nbnoe; i++)
		{
			fic << grille[i] << endl;
		}
		fic.close();				
	}

	// Fin

	delete [] grille;
	MPI_Finalize();
}


