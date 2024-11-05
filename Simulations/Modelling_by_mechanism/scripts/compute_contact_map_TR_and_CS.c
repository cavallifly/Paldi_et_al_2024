#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>
#include<unistd.h>
#include"/home/marco.di-stefano/TOOLS/My_Memory.c"
#include"/home/marco.di-stefano/TOOLS/My_IO.c"

#define MAXLINE  200

/* Fuctions */
void xyz(double **coordinates, char *filename_input, int n_particles);
void COMs(double ** COM, double **coordinates, int n_bins, int resolution, int dim);
int contact_between_points(double *P1, double *P2, int dim, double particle_contact_distance);
int contact_between_particles(double **coordinates, int start1, int stop1, int start2, int stop2, int dim, double particle_contact_distance); 
double distance_between_points(double *P1, double *P2, int dim, double particle_contact_distance);
void print_help(); 

int main(int argc,char** argv)
{
  clock_t begin, end;                            /* Self explanatory variables */
  int c, i, j, l, m, k, particle, contact, file_number;
  int chrom, chrom1, chrom2;                     /* Variables for computing the Trans-contact ratio */
  int start1, start2, end1, end2;                /* Variables for computing the Trans-contact ratio */
  int *Ccis, *Ctrans, *Ctot;                     /* Variables for computing the Trans-contact ratio */
  double *cis, *trans, *tot;                     /* Variables for computing the Trans-contact ratio */
  char compName;                                 /* Variables for computing the Trans-contact ratio */
  int comp, comp1, comp2;                        /* Variables for computing the Compartment-Strength */
  int nchains, genDist, indi, indj;              /* Variables for computing the Compartment-Strength */
  double *expected;                              /* Variables for computing the Compartment-Strength */
  int *Cexpected;                                /* Variables for computing the Compartment-Strength */
  double obsOverExp;                             /* Variables for computing the Compartment-Strength */
  int resolution, n_particles, ploidy, n_bins;   /* Self explanatory variables */
  int *mapping;                                  /* Vector to store the mapping of the particles with the reference haploid genome */
  int *flag;                                     /* Vector to store the flag to exclude the centromeres */
  int *chromosome;                               /* Vector to store the chromosome per bin */
  int *compartment;                              /* Vector to store the compartment per bin */  
  int **contact_matrix;                          /* Matrix to store the number of contacts */
  int inter, j_min, j_max;                       /* User-defined variable to consider (0) or not (1) inter copies contacts */
  double time_spent, distance;                   /* Self explanatory variables */
  double particle_contact_distance;              /* Contact distance between particles */
  double b;                                      /* Monomer diameter  */
  double M, Rg, lk, bin_contact_distance;        /* Contact distance between bins  */
  double **COM, **coordinates;
  char line[MAXLINE], chr[10], filename_input[400];  
  FILE *fp_input_files, *fp_input, *fp_output_TR, *fp_output_CS, *fp_error, *fp_contacts, *fp_chromosomes, *fp_compartments;

  fp_output_TR   = open_w("TR_per_bin.txt");
  fp_contacts = open_w("contacts.tab");
  fp_error    = open_w("output.log");  
  inter       = 0; /* Default value of the inter variable */

  /* Get options executing the file */
  if(argc < 2)
    {
      print_help();
      exit(1);
    }
  while ((c = getopt (argc, argv, ":h:r:p:k:d:b:i:")) != -1)
    switch (c)
      {
      case 'r':
        resolution  = atoi(optarg);
        fprintf(fp_error, "Resolution in particles (In Rosa et al. model one particle is about 3kb): %d\n", resolution);
        break;
      case 'p':
        n_particles = atoi(optarg);
        fprintf(fp_error, "Number of particles in the system: %d\n", n_particles);
        break;
      case 'k':
        ploidy = atoi(optarg);
        fprintf(fp_error, "Ploidy of the system (Number of chromosome copies): %d\n", ploidy);
        break;	
      case 'd':
        particle_contact_distance = atof(optarg);
        fprintf(fp_error, "Contact distance between particles in nm: %lf\n", particle_contact_distance);
        break;
      case 'b':
        b = atof(optarg);
        fprintf(fp_error, "Monomer diameter (b): %lf\n", b);
        break;
      case 'i':
        inter = atoi(optarg);
        fprintf(fp_error, "Consider (0) or not (1) inter copies contacts: %d\n", inter);
        break;	
      case '?':
        if (optopt == 'c')
          fprintf (stderr, "Option -%c requires an argument.\n", optopt);
        else if (isprint (optopt))
          fprintf (stderr, "Unknown option `-%c'.\n", optopt);
        else
          fprintf (stderr, "Unknown option character `\\x%x'.\n", optopt);
        return 1;
      default:
        abort ();
      }

  /* From the input (number of particles and the resolution in particles) 
     we can compute the number of bins */  
  n_bins = (int) (n_particles/ploidy)/resolution + 1;
  fprintf(fp_error, "Number of bins in the system: %d\n", n_bins);
  /* From the input (contact distance between particles)
     we can compute the square to ease the contact assessment. We avoid to compute
     the sqrt of the distance */  
  particle_contact_distance = pow((particle_contact_distance/b), 2.0);    
  fflush(fp_error);
  
  /* Allocating memory for the distances and check arrays */
  fprintf(fp_error, "Allocating memory: ");
  begin = clock();
  coordinates    = matrix2_d(n_particles, 3);
  mapping        = vector1_i(n_particles);
  flag           = vector1_i(n_particles);
  contact_matrix = matrix2_i(n_bins, n_bins);
  chromosome     = vector1_i(n_bins);
  compartment    = vector1_i(n_bins);
  cis            = vector1_d(n_bins);
  trans          = vector1_d(n_bins);
  tot            = vector1_d(n_bins);
  Ccis           = vector1_i(n_bins);
  Ctrans         = vector1_i(n_bins);
  Ctot           = vector1_i(n_bins);
  expected       = vector1_d(2000);
  Cexpected      = vector1_i(2000);
  /*  COM            = matrix2_d(n_bins, 3); */
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  fprintf(fp_error, "DONE in %lf\n", time_spent);
  fflush(fp_error);
  
  /* Mapping particles to the reference haploid genome */
  //fp_input = open_r("mapping_particles.txt");

  //fprintf(fp_error, "Mapping particles to the reference haploid genome\n");
  //while( fgets(line, MAXLINE, fp_input) != NULL )
  for(i=0; i<n_particles; i++)
    {
      //sscanf(line, "%s %d %d %d", chr, &i, &k);
      mapping[i] = (round) ( i % (n_particles/ploidy)) / resolution ;
      //mapping[i] = i % (n_particles/ploidy);
      flag[i]  = 1;
      //fprintf(fp_error, "Particle %d mapped on bin %d, with module %d\n", i, mapping[i], i % (n_particles/ploidy)); 
    }
  fflush(fp_error);

  /* Getting the chromosome per bin */
  fp_input = open_r("_chromosomes");

  fprintf(fp_error, "Getting the chromosome per bin\n");
  nchains = 0;
  while( fgets(line, MAXLINE, fp_input) != NULL )
    {
      sscanf(line, "%d %d", &particle, &chrom);
      chromosome[particle] = chrom;      
      //fprintf(fp_error,"%d %d\n",particle,chromosome[particle]);
      if(chrom > nchains)
	{
	  nchains = chrom;
	}
    }
  fprintf(fp_error, "Number of chromosomes %d\n", nchains);
  fclose(fp_input);
  fflush(fp_error);
  
  /* Getting the compartment per bin */
  fp_input = open_r("_compartments");

  fprintf(fp_error, "Getting the compartment per bin\n");
  while( fgets(line, MAXLINE, fp_input) != NULL )
    {
      sscanf(line, "%d %d", &particle, &comp);
      compartment[particle] = comp;
      //fprintf(fp_error,"%d %d %d\n",particle,chromosome[particle],compartment[particle]);
    }
  fclose(fp_input);
  fflush(fp_error);
  
  /* Opening the file with the names of the snapshot files */ 
  fp_input_files = open_r("DATA_FILES_INPUT.txt");

  fprintf(fp_error, "Computing the number of contacts \n");
  fflush(fp_error);

  /* Computing the number of contacts */
  /* Inizializing the counter of analysed files */  
  file_number = 0;
  /* Reading the file with input files name row by row  */
  while( (fscanf(fp_input_files, "%s", &filename_input)) != EOF )
    { 
      //fprintf(fp_output,"%s\n",filename_input);      

      /* Getting the coordinates */ 
      begin = clock();
      fprintf(fp_error, "Reading coordinates in %s: ", filename_input);
      xyz(coordinates, filename_input, n_particles);
      end = clock();
      time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      fprintf(fp_error, "DONE in %lf s\n", time_spent);
      fflush(fp_error);
      /* Calculating the centres of mass for each bin 
      fprintf(fp_error, "Computing the COMs of the bins: ");
      COMs(COM, coordinates, n_bins, resolution, 3);
      fprintf(fp_error, "DONE\n"); */

      fprintf(fp_error, "Computing the contacts between pairs of particles in file number %d: ", file_number);
      begin = clock();
      j_min = 0;
      j_max = n_particles;
      //fprintf(fp_output,"#Bin1 Bin2 dist_in_nm\n");
      for(i = 0; i < n_particles; ++i)
	{
	  if(flag[i] == 0) continue;

	  if(inter != 0)
	    {
	      j_min =  i ;
	      j_max = ((int)(i/(n_bins-1)) + 1) * (n_bins-1) ;
	      /*fprintf(fp_output,"For particle %d, I am considering j from %d in copy %d to %d in copy %d\n",i,j_min,(int)(j_min/(n_bins-1)),j_max-1,(int)((j_max-1)/(n_bins-1)));*/
	    }

	  for(j = j_min; j < j_max; ++j)
	    {
	      if(flag[j] == 0) continue;
	      /* if(inter != 0 && ((int)(i/n_bins) != (int)(j/n_bins))){fprintf(fp_output,"Not considered %d %d in copy %d vs %d in copy %d\n",inter,i,(int)(i/n_bins),j,(int)(j/n_bins)); exit(1); continue;}; 
	      if(inter != 0 && ((int)(i/(n_bins-1)) != (int)(j/(n_bins-1))))
		{
		   if(i==0) fprintf(fp_output,"Avoid %d %d in copy %d vs %d in copy %d\n",inter,i,(int)(i/(n_bins-1)),j,(int)(j/(n_bins-1)));
		  continue;
	        } */
	      
	      distance = distance_between_points(coordinates[i], coordinates[j], 3, particle_contact_distance);
	      if(distance!=-1.)
		{
		  /*fprintf(fp_output,"%d %d %lf\n",mapping[i],mapping[j],sqrt(distance)*b); */
		  contact_matrix[mapping[i]][mapping[j]] += 1;
		  if(mapping[i] != mapping[j]) contact_matrix[mapping[j]][mapping[i]] += 1;
		  //fflush(fp_output_TR);
		}
	    }
	}
      end = clock();
      time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
      fprintf(fp_error, "DONE in %lf s\n", time_spent);
      fflush(fp_error);      
      
      /* Increasing the number of analysed files */	  
      file_number = file_number + 1;

    }
  fclose(fp_input_files);

  /* Writing the contact matrix 
  fprintf(fp_error, "Writing the contact matrix: ");
  begin = clock();
  for(i = 0; i < n_bins; i++)
    {
      for(j = i; j < n_bins; j++)
	{
	  fprintf(fp_contacts, "%-7d %-7d %d\n", i, j, contact_matrix[i][j]); 
	}
    }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  fprintf(fp_error, "DONE in %lf\n", time_spent); */

  /* Compute the TR per bin */
  fprintf(fp_error, "Computing the Trans-contact ratio per bin: ");
  fflush(fp_error);  
  begin = clock();
  for(i = 0; i < n_bins-1; i++)
    {
      for(j = i; j < n_bins-1; j++)
	{
	  if(contact_matrix[i][j] == 0){continue;}
	  
	  chrom1 = chromosome[i];	  
	  chrom2 = chromosome[j];
	  
	  if(chrom1==chrom2)
	    {
	      cis[i] += contact_matrix[i][j];
	      Ccis[i]++;
	      if(i != j)
		{
		  cis[j] += contact_matrix[i][j];
		  Ccis[j]++;		  
		}	      
	    }
	  else
	    {	    
	      trans[i] += contact_matrix[i][j];
	      Ctrans[i]++;
	      if(i != j)
		{
		  trans[j] += contact_matrix[i][j];
		  Ctrans[j]++;		  
		}
	    }

	  tot[i] += contact_matrix[i][j];
	  Ctot[i]++;	  
	  if(i != j)
	    {
	      tot[j] += contact_matrix[i][j];
	      Ctot[j]++;	  	      
	    }

	  //if(contact_matrix[i][j] != 0)
	  //{
	  //fprintf(fp_error, "%d %d chr%d chr%d %lf %lf %lf %lf %lf %lf %d\n", i, j, chrom1, chrom2, cis[i], cis[j], trans[i], trans[j], tot[i], tot[j], contact_matrix[i][j]);
	  //}
	  fflush(fp_error);
	  //exit(1);
	}
    }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  fprintf(fp_error, "DONE in %lf\n", time_spent);
  fflush(fp_error);

  /* Write the TR per bin */
  fprintf(fp_error, "Writing the Trans-contact ratio per bin: ");
  fflush(fp_error);  
  begin = clock();
  for(i = 0; i < n_bins-1; i++)
    {
      chrom1 = chromosome[i];
      start1 = (i % 2000) * 10000;
      end1   = start1+10000;	  
      comp1  = compartment[i];
      //fprintf(fp_error, "%d %d %d %d %d %d %d %d\n",i,chrom1,start1,end1,comp1,cis[i],trans[i],tot[i]);
      
      if(comp1 == 0)
	{
	  compName = 'T';
	}
      if(comp1 == 1)
	{
	  compName = 'A';
	}
      if(comp1 == 2)
	{
	  compName = 'B';
	}
      if(tot[i]==0)
	{
	  fprintf(fp_output_TR, "%d_%d_%d_%d_%c %d %d %d %d %d %d %lf %lf\n", i, chrom1, start1, end1, compName, (int) cis[i], Ccis[i], (int) trans[i], Ctrans[i], (int) tot[i], Ctot[i], cis[i], trans[i]);
	}
      else
	{
	  fprintf(fp_output_TR, "%d_%d_%d_%d_%c %d %d %d %d %d %d %lf %lf\n", i, chrom1, start1, end1, compName, (int) cis[i], Ccis[i], (int) trans[i], Ctrans[i], (int) tot[i], Ctot[i], cis[i]/tot[i], trans[i]/tot[i]);	
	  fflush(fp_output_TR);  	  
	}

      //exit(1);
    }
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  fprintf(fp_error, "DONE in %lf\n", time_spent);
  /* Compute Trans-contacts per bin */

  /* Compute Compartment Strength intra-chain */
  for(chrom=0; chrom<=nchains; chrom++)
    {
      /* Re-initializing the cis and trans vectors for the Compartment-Strength calculation */
      for(i = 0; i < n_bins-1; i++)
	{
	  cis[i]    = 0.;
	  trans[i]  = 0.;
	  tot[i]    = 0.;
	  Ccis[i]   = 0;
	  Ctrans[i] = 0;
	  Ctot[i]   = 0;
	}
      for(i = 0; i < 2000; i++)
	{
	  expected[i] = 0.;
	  Cexpected[i] = 0;
	}
        
      fprintf(fp_error, "BEGIN: Computing the expected number of counts for chain %d\n", chrom);      
      for(i = 0; i < n_bins-1; i++)
	{
	  for(j = i; j < n_bins-1; j++)
	    {
	      if(chromosome[i] == chrom && chromosome[j] == chrom)
		{
		  genDist = abs(i-j);
		  if(genDist<=2){continue;}
		  /*if(contact_matrix[i][j] != 0)
		    {
		      fprintf(fp_error,"%d %d %d\n", i, j, contact_matrix[i][j]);
		      }*/

		  expected[genDist] += (double) contact_matrix[i][j];
		  Cexpected[genDist]++;
		}
	    }
	}
      for(genDist = 0; genDist < 2000; genDist++)
	{
	  if(Cexpected[genDist] == 0){continue;}
	  expected[genDist] /= Cexpected[genDist];
	  if(expected[genDist] != 0)
	    {
	      fprintf(fp_error,"%d %lf %d\n",genDist, expected[genDist], Cexpected[genDist]);	      
	    }
	}
      fprintf(fp_error, "END: Computing the expected number of counts for chain %d\n", chrom);
      fflush(fp_error);		      
      
      fprintf(fp_error, "BEGIN: Computing the Compartment-Strength per bin for chain %d\n", chrom);      
      for(i = 0; i < n_bins-1; i++)
	{
	  for(j = i; j < n_bins-1; j++)
	    {
	      if(chromosome[i] == chrom && chromosome[j] == chrom)
		{		  
		  genDist = abs(i-j);
		  /* Ignore contribution for genDist <= 2 bins */
		  if(genDist<=2){continue;}
		  if(compartment[i] == 0 || compartment[j] == 0){continue;}
		  
		  indi = (i % 2000);
		  indj = (j % 2000);
		  if(expected[genDist] != 0.0)
		    {
		      obsOverExp = (double) contact_matrix[i][j] / (double) expected[genDist];
		      //fprintf(fp_error, "obsOverExp: %d %d %lf %d %d\n", i, j, obsOverExp, contact_matrix[i][j], expected[genDist]);
		      fflush(fp_error);		      
		    }
		  else
		    {
		      obsOverExp = 0.0;
		    }
		  if(compartment[i] == compartment[j])
		    {
		      cis[indi] += obsOverExp;
		      Ccis[indi]++;

		      cis[indj] += obsOverExp;
		      Ccis[indj]++;
		    }
		  tot[indi] += obsOverExp;
		  Ctot[indi]++;
		  
		  tot[indj] += obsOverExp;
		  Ctot[indj]++;

		  //fprintf(fp_error, "%d %d %d %d %lf %d %lf %d %lf %d %lf %d\n", i, j, indi, indj, cis[indi], Ccis[indi], cis[indj], Ccis[indj], tot[indi], Ctot[indi], tot[indj], Ctot[indj]);
		}
	    }
	}
      fprintf(fp_error, "END: Computing the Compartment-Strength per bin for chain %d\n", chrom);
  
      /* Write Compartment Strength intra-chain */
      sprintf(line,"CS_per_bin_%d.txt", chrom);
      fp_output_CS   = open_w(line);
      for(i = 0; i < 2000; i++)
	{
	  comp1 = compartment[i];
	  
	  if(comp1 == 0){continue;};
	  if(comp1 == 1)
	    {
	      compName = 'A';
	    }
	  if(comp1 == 2)
	    {
	      compName = 'B';
	    }
	  if(tot[i] == 0)
	    {
	      fprintf(fp_output_CS, "%d %c %lf %d %lf %d %lf\n", i, compName, cis[i], Ccis[i],tot[i],Ctot[i],0.);
	    }
	  else
	    {	  
	      fprintf(fp_output_CS, "%d %c %lf %d %lf %d %lf\n", i, compName, cis[i], Ccis[i],tot[i],Ctot[i],(cis[i]/Ccis[i])/(tot[i]/Ctot[i]));
	    }
	  
	}
      //exit(1);
      fclose(fp_output_CS);
    }
  
  /* Free Memory */
  fprintf(fp_error, "Freeing memory: ");
  fflush(fp_error);
  begin = clock();
  free1_i(mapping);
  free1_i(flag); 
  free2_d(coordinates, n_particles); 
  free2_i(contact_matrix, n_bins);
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  fprintf(fp_error, "DONE in %lf\n", time_spent);

  /* Closing Files */
  fprintf(fp_error, "Closing Files: ");
  begin = clock();
  fclose(fp_output_TR);
  fclose(fp_contacts);
  end = clock();
  time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
  fprintf(fp_error, "DONE in %lf\n", time_spent);
  fflush(fp_error);
  fclose(fp_error);
   
  return(0);
}

/*************************/

void print_help()
{
  fprintf(stderr," \n");
  fprintf(stderr,"OPTIONS: \n");
  fprintf(stderr,"\t-h:\t\t Print this help and exit.\n");
  fprintf(stderr,"\t-r:\t\t Resolution of the map in monomers\n");
  fprintf(stderr,"\t-p:\t\t Number of particles in the system\n");
  fprintf(stderr,"\t-d:\t\t Contact distance cutoff between particles (nm)\n");
  fprintf(stderr,"\t-b:\t\t Monomer diameter (nm)\n");
  fprintf(stderr,"\t-k:\t\t Ploidy of the system\n");
  fprintf(stderr,"\t-i:\t\t Consider (0) or not (1) inter copies contacts\n");
  fprintf(stderr," \n");
  fprintf(stderr,"INPUT:   \n");
  fprintf(stderr,"\t File DATA_FILES_INPUT.TXT contains the absolute paths of the coordinate files\n");
  fprintf(stderr," \n");
  fprintf(stderr,"OUTPUT:   \n");
  fprintf(stderr,"\t File output.txt contains the computed contact map\n");
  fprintf(stderr,"\t File output.log contains some details of the computation\n");
  fprintf(stderr," \n");
  exit(1);
}

/*************************/

void xyz(double **coordinates, char *filename_input, int n_particles)
{
  int i, j, type, n_atoms, particle; 
  double *x;
  char line[200];
  FILE *fp_input;

  /* Allocating Memory*/  
  x           = vector1_d(3);
  
  /* Apertura del file di input */
  fp_input = open_r(filename_input);
  
  /* Skip the first 9 lines */
  for(i = 0; i < 9 ; ++i) fgets(line, MAXLINE, fp_input);

  while( fgets(line, MAXLINE, fp_input) != NULL)
    {
      //sscanf(line, "%i %i %lf %lf %lf", &i, &type, &x[0], &x[1], &x[2]);
      sscanf(line, "%i %lf %lf %lf", &i, &x[0], &x[1], &x[2]);
      for(j = 0; j < 3; ++j) coordinates[i-1][j] = x[j];
      //      fprintf(stderr,"%d %lf %lf %lf %lf %lf %lf\n", i, x[0], x[1], x[2], coordinates[i-1][0], coordinates[i-1][1], coordinates[i-1][2]);
    }
  fclose(fp_input);
  
  /* Free Memory */
  free1_d(x);

}

/*************************/

void COMs(double **COM, double **coordinates, int n_bins, int resolution, int dim)
{
  int n; /* Variables for the cycles on the number of bins between 0 and n_bins-1     */
  int i; /* Variables for the cycles within the bin between 0 and resolution-1        */
  int k; /* Variables for the cycles on the Cartesian coordinates between 0 and dim-1 */
  
  for(n = 1; n < n_bins+1; ++n)
    {
      /*      fprintf(stderr,"%d %d\n", (n-1)*resolution, n*resolution-1);*/
      for(i = (n-1)*resolution; i < n*resolution-1; ++i)
	{
	  for(k = 0; k < 3; ++k)
	    {
	      COM[n-1][k] += coordinates[i][k];
	    }
	}
      for(k = 0; k < 3; ++k)
	{
	  COM[n-1][k] = COM[n-1][k] / resolution;
	}
      /*      fprintf(stderr,"%d %lf %lf %lf\n", n, COM[n-1][0], COM[n-1][1], COM[n-1][2]); */
    }
}

/*************************/

int contact_between_points(double *P1, double *P2, int dim, double particle_contact_distance)
{
  int    k;
  double distance;

  /* Initializing the distance to 0.0 */
  distance = 0.0;
  for(k = 0; k < dim; ++k)
    {
      distance += (P1[k]-P2[k])*(P1[k]-P2[k]);
      if(distance > particle_contact_distance) return(0);
    }

  return(1);
}

/*************************/

int contact_between_particles(double **coordinates, int start1, int stop1, int start2, int stop2, int dim, double particle_contact_distance)
{
  int i, j;

  for(i = start1; i < stop1; ++i)
    {
      for(j = start2; j < stop2; ++j)
	{
	  if(contact_between_points(coordinates[i], coordinates[j], 3, particle_contact_distance) == 1) return(1);
	}
    }

  return(0);
}

/*************************/

double distance_between_points(double *P1, double *P2, int dim, double particle_contact_distance)
{
  int    k;
  double distance;

  /* Initializing the distance to 0.0 */
  distance = 0.0;
  for(k = 0; k < dim; ++k)
    {
      distance += (P1[k]-P2[k])*(P1[k]-P2[k]);
      if(distance > particle_contact_distance) return(-1.);
    }
  if(distance > particle_contact_distance) return(-1.);

  return(distance);
}
